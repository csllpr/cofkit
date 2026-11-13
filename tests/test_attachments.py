"""Attachment metadata preserves the chosen atomic periodic connectivity."""

from collections import Counter
from dataclasses import replace
import json

import gemmi
import pytest

from cofkit import (
    AttachmentRecipe,
    BatchGenerationConfig,
    BatchStructureGenerator,
    CIFWriter,
    COFEngine,
    COFProject,
    build_from_attachment,
    build_rdkit_monomer,
    decompose_cif_to_cofid,
    load_attachment_recipe,
    verify_attachment_cif,
)
from cofkit.attachment import AttachmentBond, compare_attachment_graphs
from cofkit.build_workflows.ring_forming import RingFormingStructureGenerator
from cofkit.cli import main
from cofkit.decompose import _build_bonded_mol
from cofkit.decompose_cif import read_periodic_cif_atoms
from cofkit.decompose_events import (
    _generate_and_validate_hypotheses,
    _select_event_result,
    detect_linkage_events,
)
from cofkit.periodic_geometry import format_p1_shift, p1_shift
from test_decompose import ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES


FAMILIES = (
    "imine",
    "hydrazone",
    "azine",
    "bken",
    "boest",
    "vinylene",
    "boroxine",
    "triazine",
)


@pytest.fixture(scope="module", params=FAMILIES)
def source(request, tmp_path_factory):
    family = request.param
    root = tmp_path_factory.mktemp(f"attachment_{family}")
    path = root / "source.cif"
    if family in {"boroxine", "triazine"}:
        kind, smiles, template = (
            ("boronic_acid", "OB(O)c1ccc(B(O)O)cc1", "boroxine_trimerization")
            if family == "boroxine"
            else ("nitrile", "N#Cc1ccc(C#N)cc1", "triazine_trimerization")
        )
        monomers = (build_rdkit_monomer("p", "p", smiles, kind, num_conformers=1),)
        candidate = RingFormingStructureGenerator().generate(monomers[0], template)
    else:
        template, _, first, second, _ = next(
            c for c in ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES if c[1] == family
        )
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=(template,),
                topology_ids=("hcb",),
                rdkit_num_conformers=1,
                max_workers=1,
            )
        )
        _, candidate = generator.generate_pair_candidate(first, second, write_cif=False)
        assert candidate is not None
        monomers = tuple(generator.build_monomer(r).monomer for r in (first, second))
    CIFWriter().write_candidate(path, candidate, monomers)
    result = decompose_cif_to_cofid(path, topology="hcb", linkage=family)
    assert result.ok, result.reason
    assert result.metadata["attachment"]["status"] == "available", result.metadata[
        "attachment"
    ]
    return family, path, result, monomers


def test_all_supported_linkages_replay_source_heavy_graph(source, tmp_path):
    family, _, result, monomers = source
    recipe = result.attachment_recipe
    assert isinstance(recipe, AttachmentRecipe)
    assert recipe.cofid == result.cofid
    assert recipe.provenance["bond_source"] == "explicit_cif"
    assert recipe.provenance["experimental_status"] == "not_asserted"
    assert len(recipe.instances) == len(
        result.metadata["atom_ledger"]["recovered_fragments"]
    )
    assert all(
        len(e.participants) == (3 if family in {"boroxine", "triazine"} else 2)
        for e in recipe.events
    )
    grouped = Counter(e.source_event_id for e in recipe.events)
    assert set(grouped.values()) == ({2} if family == "azine" else {1})
    # Pass the original, sometimes noncanonical-SMILES atom order. The recipe
    # must bind a full molecular correspondence rather than matching indices.
    replay = build_from_attachment(recipe, monomers)
    assert (
        replay.candidate.metadata["attachment"]["verification"]["status"] == "verified"
    )
    assert replay.candidate.state.cell == recipe.cell
    assert len(replay.candidate.state.monomer_poses) == len(recipe.instances)
    path = tmp_path / "replayed.cif"
    CIFWriter().write_candidate(
        path, replay.candidate, replay.monomers, cofid=result.cofid
    )
    assert verify_attachment_cif(recipe, path, replay.monomers)["status"] == "verified"
    assert (
        decompose_cif_to_cofid(path, topology="hcb", linkage=family).cofid
        == result.cofid
    )


def test_recipe_json_round_trip_and_all_loader_surfaces(source, tmp_path):
    _, _, result, _ = source
    recipe = result.attachment_recipe
    for payload in (recipe.to_dict(), result.metadata["attachment"], result.to_dict()):
        path = tmp_path / "recipe.json"
        path.write_text(json.dumps(payload))
        restored = load_attachment_recipe(path)
        assert restored.to_dict() == recipe.to_dict()
        assert restored.fingerprint == recipe.fingerprint


def test_cli_decomposition_to_constrained_build(source, tmp_path, capsys):
    family, path, result, _ = source
    sidecar = tmp_path / "source.attachment.json"
    main(
        [
            "analyze",
            "decompose",
            str(path),
            "--topology",
            "hcb",
            "--linkage",
            family,
            "--attachment-output",
            str(sidecar),
            "--json",
        ]
    )
    report = json.loads(capsys.readouterr().out)
    assert report["cofid"] == result.cofid
    command = "ring-forming" if family in {"boroxine", "triazine"} else "single-pair"
    output = tmp_path / "build"
    main(
        [
            "build",
            command,
            "--attachment-file",
            str(sidecar),
            "--output-dir",
            str(output),
            "--num-conformers",
            "1",
        ]
    )
    summary = json.loads((output / "summary.json").read_text())
    assert summary["requested_cofid"] == result.cofid
    assert load_attachment_recipe(output / "attachment.json").cofid == result.cofid
    row = summary["results"][0]
    attachment = row.get("attachment") or row["metadata"]["attachment"]
    assert attachment["status"] == "applied"
    assert attachment["verification"]["status"] == "verified"
    # Geometry can independently withhold the artifact; it cannot replace the
    # requested attachment by an unconstrained topology-generated candidate.
    assert row["cif_export"]["status"] in {"written", "blocked"}


def test_engine_and_batch_apis_honor_recipe(source, tmp_path):
    family, _, result, monomers = source
    recipe = result.attachment_recipe
    template = recipe.events[0].template_id
    ensemble = COFEngine().run(
        COFProject(
            monomers, (template,), target_topologies=("hcb",), attachment_recipe=recipe
        )
    )
    assert len(ensemble.candidates) == 1
    assert (
        ensemble.candidates[0].metadata["attachment"]["verification"]["status"]
        == "verified"
    )
    if family in {"boroxine", "triazine"}:
        generated = RingFormingStructureGenerator().generate(
            monomers[0], template, attachment_recipe=recipe
        )
        assert generated.metadata["attachment"]["verification"]["status"] == "verified"
    else:
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                allowed_reactions=(template,),
                topology_ids=("hcb",),
                attachment_recipe=recipe,
                max_workers=1,
            )
        )
        summary, candidate = generator.generate_monomer_pair_candidate(
            *monomers, out_dir=tmp_path
        )
        assert candidate is not None
        assert summary.metadata["attachment"]["verification"]["status"] == "verified"
        if summary.cif_path:
            assert (
                summary.metadata["validation"]["attachment_verification"]["status"]
                == "verified"
            )


def test_seed_positions_can_be_replaced_by_an_explicit_initial_state(source):
    _, _, result, monomers = source
    original = build_from_attachment(result.attachment_recipe, monomers)
    recipe = replace(
        result.attachment_recipe,
        instances=tuple(
            replace(i, seed_positions=()) for i in result.attachment_recipe.instances
        ),
    )
    replay = build_from_attachment(
        recipe, monomers, initial_state=original.candidate.state
    )
    assert (
        replay.candidate.metadata["attachment"]["verification"]["status"] == "verified"
    )


def _copy_cif(source, destination, *, repeat=1, rebase=False):
    """Independent literal graph-preserving cell expansion / atom rebasing."""
    block = gemmi.cif.read_file(str(source)).sole_block()
    small = gemmi.make_small_structure_from_block(block)
    sites = list(small.sites)
    shifts = {
        s.label: ((i % 3) - 1, (i % 2), 0) if rebase else (0, 0, 0)
        for i, s in enumerate(sites)
    }
    rows = []
    for tile in range(repeat):
        for site in reversed(sites):
            offset = shifts[site.label]
            rows.append(
                f"t{tile}_{site.label} {site.element.name} {(site.fract.x + tile) / repeat + offset[0]} {site.fract.y + offset[1]} {site.fract.z + offset[2]} 1"
            )
    bonds = []
    for tile in range(repeat):
        for a, b, s1, s2, d, kind in block.find(
            [
                "_geom_bond_atom_site_label_1",
                "_geom_bond_atom_site_label_2",
                "_geom_bond_site_symmetry_1",
                "_geom_bond_site_symmetry_2",
                "_geom_bond_distance",
                "_ccdc_geom_bond_type",
            ]
        ):
            a, b = str(a), str(b)
            g = tuple(y - x for x, y in zip(p1_shift(str(s1)), p1_shift(str(s2))))
            target = (tile + g[0]) % repeat
            image = ((tile + g[0]) // repeat, g[1], g[2])
            image = tuple(v + x - y for v, x, y in zip(image, shifts[a], shifts[b]))
            bonds.append(
                f"t{tile}_{a} t{target}_{b} . {format_p1_shift(image)} {d} {kind}"
            )
    destination.write_text(
        f"""data_representatives
_cell_length_a {small.cell.a * repeat}
_cell_length_b {small.cell.b}
_cell_length_c {small.cell.c}
_cell_angle_alpha {small.cell.alpha}
_cell_angle_beta {small.cell.beta}
_cell_angle_gamma {small.cell.gamma}
_space_group_name_H-M_alt 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
"""
        + "\n".join(rows)
        + """
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_site_symmetry_1
_geom_bond_site_symmetry_2
_geom_bond_distance
_ccdc_geom_bond_type
"""
        + "\n".join(reversed(bonds))
        + "\n"
    )


def test_atom_order_wrapping_and_supercells_preserve_attachment(source, tmp_path):
    family, path, result, _ = source
    modified = tmp_path / "expanded.cif"
    _copy_cif(path, modified, repeat=2, rebase=True)
    recovered = decompose_cif_to_cofid(modified, topology="hcb", linkage=family)
    assert recovered.ok, recovered.reason
    assert recovered.cofid == result.cofid
    recipe = recovered.attachment_recipe
    assert recipe is not None, recovered.metadata["attachment"]
    assert len(recipe.instances) == 2 * len(result.attachment_recipe.instances)
    replay = build_from_attachment(recipe)
    assert (
        replay.candidate.metadata["attachment"]["verification"]["status"] == "verified"
    )


def test_attachment_replay_under_unimodular_cell_change(source):
    _, _, result, monomers = source
    recipe = result.attachment_recipe

    def transform(image):
        return image[0], image[1] - image[0], image[2]

    # a' = a + b and b' = b: transform every periodic event/reference gain
    # and fractional ring center, keeping Cartesian placement hints fixed.
    events = []
    for event in recipe.events:
        metadata = dict(event.metadata)
        if "ring_center_fractional" in metadata:
            metadata["ring_center_fractional"] = transform(
                metadata["ring_center_fractional"]
            )
        events.append(
            replace(
                event,
                participants=tuple(
                    replace(s, periodic_image=transform(s.periodic_image))
                    for s in event.participants
                ),
                metadata=metadata,
            )
        )
    changed = replace(
        recipe,
        cell=(
            tuple(a + b for a, b in zip(recipe.cell[0], recipe.cell[1])),
            recipe.cell[1],
            recipe.cell[2],
        ),
        events=tuple(events),
        reference_bonds=tuple(
            replace(b, image=transform(b.image)) for b in recipe.reference_bonds
        ),
    )
    replay = build_from_attachment(changed, monomers)
    assert (
        replay.candidate.metadata["attachment"]["verification"]["status"] == "verified"
    )


def test_missing_recipe_mapping_is_nonfatal(source, monkeypatch):
    family, path, result, _ = source

    def fail(*args, **kwargs):
        raise ValueError("synthetic canonical correspondence failure")

    monkeypatch.setattr("cofkit.attachment_decompose.extract_attachment_recipe", fail)
    recovered = decompose_cif_to_cofid(path, topology="hcb", linkage=family)
    assert recovered.ok
    assert recovered.cofid == result.cofid
    assert recovered.attachment_recipe is None
    assert (
        "synthetic canonical correspondence failure"
        in recovered.metadata["attachment"]["error"]
    )
    with pytest.raises(ValueError, match="available recipe"):
        load_attachment_recipe(recovered.to_dict())


def test_same_cofid_alternatives_need_attachment_equivalence(source):
    family, path, _, _ = source
    atoms = read_periodic_cif_atoms(path)
    built = _build_bonded_mol(atoms)
    detection = detect_linkage_events(built)
    hypotheses, generation = _generate_and_validate_hypotheses(
        path, atoms, built, detection.events, topology="hcb", bond_mode="auto"
    )
    selected = next(h for h in hypotheses if h.complete and h.family == family)
    report = json.loads(json.dumps(selected.metadata["attachment"]))
    report["recipe"]["provenance"]["source_graph_sha256"] = "different source graph"
    different = replace(
        selected,
        hypothesis_id="different-attachment",
        metadata={**selected.metadata, "attachment": report},
    )
    result = _select_event_result(
        path,
        requested_family=family,
        topology="hcb",
        detection=detection,
        hypotheses=(selected, different),
        generation_metadata=generation,
    )
    assert result.ok  # compact precursor identity remains available
    assert result.metadata["attachment"]["status"] == "ambiguous"
    assert result.attachment_recipe is None
    assert len(result.metadata["attachment"]["alternatives"]) == 2
    with pytest.raises(ValueError, match="available recipe"):
        load_attachment_recipe(result.to_dict())


def test_periodic_graph_check_keeps_multiedges_and_self_images():
    atoms = {("a", 0): "C", ("b", 0): "N"}
    edges = (
        AttachmentBond(("a", 0), ("b", 0)),
        AttachmentBond(("a", 0), ("b", 0), (1, 0, 0)),
        AttachmentBond(("a", 0), ("a", 0), (0, 1, 0)),
    )
    rebased = (
        AttachmentBond(("a", 0), ("b", 0), (3, -2, 0)),
        AttachmentBond(("a", 0), ("b", 0), (4, -2, 0)),
        edges[2],
    )
    assert (
        compare_attachment_graphs(atoms, edges, atoms, rebased)["status"] == "verified"
    )
    assert (
        compare_attachment_graphs(atoms, edges, atoms, rebased[:2])["status"]
        == "mismatch"
    )
    altered = (*rebased[:2], replace(edges[2], image=(0, 0, 1)))
    assert (
        compare_attachment_graphs(atoms, edges, atoms, altered)["status"] == "mismatch"
    )
