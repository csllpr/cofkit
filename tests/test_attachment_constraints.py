"""Connectivity alternatives, invalid recipes, and full-cell replay contracts."""

from collections import Counter, defaultdict
from dataclasses import replace
import json

import gemmi
import pytest
from rdkit import Chem

from cofkit import (
    BatchGenerationConfig,
    BatchStructureGenerator,
    CIFWriter,
    build_from_attachment,
    build_rdkit_monomer,
    decompose_cif_to_cofid,
    load_attachment_recipe,
)
from cofkit.cli import main
from cofkit.cofid import cofid_to_build_request
from test_attachments import _copy_cif


@pytest.fixture(scope="module")
def asymmetric_source(tmp_path_factory):
    root = tmp_path_factory.mktemp("asymmetric_attachments")
    amine = build_rdkit_monomer(
        "a",
        "a",
        "Nc1ccc(C(=C(c2ccc(N)cc2)c2ccc(N)cc2)c2ccc(N)cc2)cc1",
        "amine",
        num_conformers=1,
    )
    aldehyde = build_rdkit_monomer(
        "b", "b", "O=Cc1ccc(C=O)cc1", "aldehyde", num_conformers=1
    )
    generator = BatchStructureGenerator(
        BatchGenerationConfig(topology_ids=("kgm",), max_workers=1)
    )
    summary, candidate = generator.generate_monomer_pair_candidate(
        amine, aldehyde, write_cif=False
    )
    path = root / "source.cif"
    CIFWriter().write_candidate(
        path, candidate, (amine, aldehyde), cofid=summary.metadata["cofid"]
    )
    result = decompose_cif_to_cofid(path, topology="kgm", linkage="imine")
    assert result.attachment_recipe is not None, result.metadata.get("attachment")
    return path, result, (amine, aldehyde)


def _inequivalent_swap(recipe):
    request = cofid_to_build_request(recipe.cofid)
    instance = next(
        i
        for i in recipe.instances
        if request.monomers[i.monomer_index].connectivity == 4
    )
    mol = Chem.MolFromSmiles(request.monomers[instance.monomer_index].canonical_smiles)
    sites = sorted(
        {
            a
            for e in recipe.events
            for s in e.participants
            if s.instance_id == instance.id
            for a in s.atom_indices
        }
    )
    for a in sites:
        for b in sites:
            if a >= b:
                continue
            swap = {a: b, b: a}
            if any(
                len(Chem.GetShortestPath(mol, x, y))
                != len(Chem.GetShortestPath(mol, swap.get(x, x), swap.get(y, y)))
                for x in sites
                for y in sites
                if x != y
            ):
                return replace(
                    recipe,
                    events=tuple(
                        replace(
                            e,
                            participants=tuple(
                                replace(
                                    s,
                                    atom_indices=tuple(
                                        swap.get(x, x) for x in s.atom_indices
                                    ),
                                )
                                if s.instance_id == instance.id
                                else s
                                for s in e.participants
                            ),
                        )
                        for e in recipe.events
                    ),
                )
    raise AssertionError("fixture lacks an inequivalent attachment permutation")


def _triangle_corner_distances(recipe):
    """Independent chemical invariant of the decorated kgm net."""
    request = cofid_to_build_request(recipe.cofid)
    instances = {i.id: i for i in recipe.instances}
    linkers = defaultdict(list)
    for event in recipe.events:
        a, b = event.participants
        if request.monomers[instances[a.instance_id].monomer_index].connectivity == 2:
            a, b = b, a
        linkers[b.instance_id].append(
            (
                a.instance_id,
                a.atom_indices[0],
                tuple(y - x for x, y in zip(a.periodic_image, b.periodic_image)),
            )
        )
    adjacency = defaultdict(list)
    for ends in linkers.values():
        (a, pa, ga), (b, pb, gb) = ends
        gain = tuple(x - y for x, y in zip(ga, gb))
        adjacency[a].append((b, gain, pa))
        adjacency[b].append((a, tuple(-g for g in gain), pb))
    distances = []
    for node, edges in adjacency.items():
        mol = Chem.MolFromSmiles(
            request.monomers[instances[node].monomer_index].canonical_smiles
        )
        for i, (a, ga, pa) in enumerate(edges):
            for b, gb, pb in edges[i + 1 :]:
                if any(
                    other == b and tuple(x + y for x, y in zip(ga, g)) == gb
                    for other, g, _ in adjacency[a]
                ):
                    distances.append(len(Chem.GetShortestPath(mol, pa, pb)) - 1)
    return Counter(distances)


def test_distinct_attachment_patterns_with_one_cofid_replay_independently(
    asymmetric_source, tmp_path
):
    _, original, monomers = asymmetric_source
    recipe = original.attachment_recipe
    alternative = replace(
        _inequivalent_swap(recipe),
        reference_atoms=(),
        reference_bonds=(),
        provenance={"source": "modeling_choice", "experimental_status": "not_asserted"},
    )
    assert _triangle_corner_distances(recipe) != _triangle_corner_distances(alternative)
    built = build_from_attachment(alternative, monomers)
    assert (
        built.candidate.metadata["attachment"]["verification"]["status"]
        == "not_checked"
    )
    source = tmp_path / "alternative.cif"
    CIFWriter().write_candidate(source, built.candidate, built.monomers)
    recovered = decompose_cif_to_cofid(source, topology="kgm", linkage="imine")
    assert recovered.cofid == original.cofid
    assert recovered.attachment_recipe is not None, recovered.metadata["attachment"]
    assert _triangle_corner_distances(
        recovered.attachment_recipe
    ) == _triangle_corner_distances(alternative)
    replay = build_from_attachment(recovered.attachment_recipe, monomers)
    assert (
        replay.candidate.metadata["attachment"]["verification"]["status"] == "verified"
    )


def test_distance_bond_provenance_does_not_claim_experimental_resolution(
    asymmetric_source,
):
    path, original, monomers = asymmetric_source
    inferred = decompose_cif_to_cofid(
        path, topology="kgm", linkage="imine", bond_mode="distance"
    )
    assert inferred.cofid == original.cofid
    recipe = inferred.attachment_recipe
    assert recipe.provenance["bond_source"] == "distance_inferred"
    assert recipe.provenance["experimental_status"] == "not_asserted"
    assert (
        build_from_attachment(recipe, monomers).candidate.metadata["attachment"][
            "verification"
        ]["status"]
        == "verified"
    )


def test_legacy_sidecar_is_explicitly_unavailable(asymmetric_source, tmp_path, capsys):
    path, _, _ = asymmetric_source
    output = tmp_path / "legacy.attachment.json"
    main(
        [
            "analyze",
            "decompose",
            str(path),
            "--topology",
            "kgm",
            "--linkage",
            "imine",
            "--decomposition-mode",
            "legacy",
            "--attachment-output",
            str(output),
        ]
    )
    assert json.loads(output.read_text())["status"] == "unavailable"
    assert "status=unavailable" in capsys.readouterr().err
    with pytest.raises(ValueError, match="available recipe"):
        load_attachment_recipe(output)


@pytest.mark.parametrize(
    "data",
    [
        {"metadata": None},
        {"metadata": {"attachment": None}},
        [],
        {"status": "ambiguous", "recipe": None},
    ],
)
def test_malformed_or_ambiguous_envelopes_are_not_build_inputs(data):
    with pytest.raises(ValueError):
        load_attachment_recipe(data)


def test_modified_sites_cannot_silently_override_source_reference(asymmetric_source):
    _, result, monomers = asymmetric_source
    with pytest.raises(ValueError, match="reference graph"):
        build_from_attachment(_inequivalent_swap(result.attachment_recipe), monomers)


@pytest.mark.parametrize(
    "mutation",
    [
        "version",
        "fractional_image",
        "duplicate_instance",
        "missing_instance",
        "wrong_template",
        "negative_atom",
        "nonfinite_cell",
        "boolean_index",
        "wrong_element",
    ],
)
def test_invalid_recipe_rejected(asymmetric_source, mutation):
    _, result, _ = asymmetric_source
    data = result.attachment_recipe.to_dict()
    if mutation == "version":
        data["schema_version"] = 999
    elif mutation == "fractional_image":
        data["events"][0]["participants"][0]["periodic_image"][0] = 0.5
    elif mutation == "duplicate_instance":
        data["instances"].append(data["instances"][0])
    elif mutation == "missing_instance":
        data["events"][0]["participants"][0]["instance_id"] = "missing"
    elif mutation == "wrong_template":
        data["events"][0]["template_id"] = "vinylene_bridge"
    elif mutation == "negative_atom":
        data["events"][0]["participants"][0]["atom_indices"] = [-1]
    elif mutation == "nonfinite_cell":
        data["cell"][0][0] = float("nan")
    elif mutation == "boolean_index":
        data["instances"][0]["monomer_index"] = True
    else:
        data["reference_atoms"][0]["element"] = "Xe"
    with pytest.raises(ValueError):
        load_attachment_recipe(data)


def test_duplicate_reactive_site_rejected(asymmetric_source):
    _, result, monomers = asymmetric_source
    recipe = result.attachment_recipe
    first = recipe.events[0]
    same_instance = next(
        e
        for e in recipe.events[1:]
        if any(
            s.instance_id == first.participants[0].instance_id for s in e.participants
        )
    )
    replaced = replace(
        same_instance,
        participants=tuple(
            first.participants[0]
            if s.instance_id == first.participants[0].instance_id
            else s
            for s in same_instance.participants
        ),
    )
    modified = replace(
        recipe,
        events=tuple(
            replaced if e.id == same_instance.id else e for e in recipe.events
        ),
    )
    with pytest.raises(ValueError, match="already consumed"):
        build_from_attachment(modified, monomers)


def test_cli_rejects_cofid_or_stacking_conflicts(asymmetric_source, tmp_path):
    _, result, _ = asymmetric_source
    path = tmp_path / "recipe.json"
    path.write_text(json.dumps(result.attachment_recipe.to_dict()))
    for flags in (
        ["--stacking", "AA"],
        ["--cofid", result.cofid.replace("kgm", "sql")],
    ):
        with pytest.raises(SystemExit, match="error:"):
            main(
                [
                    "build",
                    "single-pair",
                    "--attachment-file",
                    str(path),
                    "--output-dir",
                    str(tmp_path / "build"),
                    *flags,
                ]
            )


def test_actual_written_connectivity_is_verified(
    asymmetric_source, monkeypatch, tmp_path
):
    _, result, monomers = asymmetric_source
    recipe = result.attachment_recipe
    original = CIFWriter.write_candidate

    def corrupt(self, path, *args, **kwargs):
        output = original(self, path, *args, **kwargs)
        block = gemmi.cif.read_file(str(path)).sole_block()
        labels = block.find_loop("_geom_bond_atom_site_label_1")
        symmetries = block.find_loop("_geom_bond_site_symmetry_2")
        for i, label in enumerate(labels):
            if "_H" not in str(label):
                symmetries[i] = "1_655" if str(symmetries[i]) != "1_655" else "."
                break
        block.write_file(str(path))
        return output

    monkeypatch.setattr(CIFWriter, "write_candidate", corrupt)
    generator = BatchStructureGenerator(
        BatchGenerationConfig(
            topology_ids=("kgm",), attachment_recipe=recipe, max_workers=1
        )
    )
    summary, _ = generator.generate_monomer_pair_candidate(*monomers, out_dir=tmp_path)
    assert summary.status == "export-failed"
    assert summary.cif_export.status == "failed"
    assert "attachment reference graph" in summary.cif_export.error
    assert summary.cif_path is None


def test_attachment_mismatch_does_not_abort_batch(asymmetric_source, tmp_path):
    _, result, monomers = asymmetric_source
    library = tmp_path / "library"
    library.mkdir()
    (library / "amines_count_4.txt").write_text(
        "smiles\n"
        + monomers[0].metadata["source_smiles"]
        + "\nNc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1\n"
    )
    (library / "aldehydes_count_2.txt").write_text("smiles\nO=Cc1ccc(C=O)cc1\n")
    generator = BatchStructureGenerator(
        BatchGenerationConfig(
            topology_ids=("kgm",),
            attachment_recipe=result.attachment_recipe,
            rdkit_num_conformers=1,
            max_workers=1,
        )
    )
    summary = generator.run_binary_bridge_batch(library, tmp_path / "build")
    rows = [
        json.loads(line)
        for line in (tmp_path / "build/manifest.jsonl").read_text().splitlines()
    ]
    assert len(rows) == 2
    assert sorted(r["status"] for r in rows) == ["generation-failed", "ok"]
    assert len(summary.record_failures) == 1
    assert "ValueError:" in next(iter(summary.record_failures.values()))


def test_three_dimensional_cell_multiplicity_is_preserved(tmp_path):
    a = build_rdkit_monomer(
        "a",
        "a",
        "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1",
        "amine",
        num_conformers=1,
    )
    b = build_rdkit_monomer("b", "b", "O=Cc1ccc(C=O)cc1", "aldehyde", num_conformers=1)
    generator = BatchStructureGenerator(
        BatchGenerationConfig(
            target_dimensionality="3D", topology_ids=("dia",), max_workers=1
        )
    )
    _, candidate = generator.generate_monomer_pair_candidate(a, b, write_cif=False)
    original, expanded = tmp_path / "original.cif", tmp_path / "expanded.cif"
    CIFWriter().write_candidate(original, candidate, (a, b))
    before = decompose_cif_to_cofid(original, topology="dia", linkage="imine")
    _copy_cif(original, expanded, repeat=2, rebase=True)
    after = decompose_cif_to_cofid(expanded, topology="dia", linkage="imine")
    assert before.ok and after.ok
    assert before.cofid == after.cofid
    assert len(after.attachment_recipe.instances) == 2 * len(
        before.attachment_recipe.instances
    )
    replay = build_from_attachment(after.attachment_recipe, (a, b))
    assert (
        replay.candidate.metadata["attachment"]["verification"]["status"] == "verified"
    )


@pytest.mark.parametrize(
    "family,smiles,kind",
    [
        ("triazine", "N#Cc1cc(C#N)cc(C#N)c1", "nitrile"),
        ("boroxine", "OB(O)c1cc(B(O)O)cc(B(O)O)c1", "boronic_acid"),
    ],
)
def test_ring_uses_distinct_periodic_copies_of_one_higher_connected_precursor(
    tmp_path, family, smiles, kind
):
    from cofkit import RingFormingStructureGenerator

    monomer = build_rdkit_monomer("p", "p", smiles, kind, num_conformers=1)
    candidate = RingFormingStructureGenerator().generate(
        monomer, f"{family}_trimerization"
    )
    path = tmp_path / "ring.cif"
    CIFWriter().write_candidate(path, candidate, (monomer,))
    recovered = decompose_cif_to_cofid(path, topology="hcb", linkage=family)
    recipe = recovered.attachment_recipe
    assert recipe is not None, recovered.metadata["attachment"]
    assert len(recipe.instances) == 1
    for event in recipe.events:
        assert len({(s.instance_id, s.periodic_image) for s in event.participants}) == 3
    replay = build_from_attachment(recipe, (monomer,))
    assert (
        replay.candidate.metadata["attachment"]["verification"]["status"] == "verified"
    )
