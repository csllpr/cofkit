"""Reproducible common-COF precursor seeds for pair-search comparisons.

These models exercise all six binary families and both ring builders. They
are construction seeds, not experimental coordinates or physical minima.
Only installed cofkit APIs and sibling fixture paths are used.
"""

from dataclasses import asdict, replace
import json
from pathlib import Path

import gemmi

from cofkit.batch import BatchGenerationConfig, BatchStructureGenerator
from cofkit.build_workflows.ring_forming import RingFormationConfig, RingFormingStructureGenerator
from cofkit.chem.rdkit import build_rdkit_monomer
from cofkit.cif import CIFWriter
from cofkit.cofid import cofid_to_build_request
from cofkit.seed_placement import SeedPlacementConfig


COF_CASES = json.loads((Path(__file__).parent / "fixtures/soft_relax/cofs.json").read_text())
# Heuristic — pending calibration: inexpensive fixed conformers for this
# search-equivalence cohort; geometry accuracy is not a benchmark label.
COF_TEST_CONFORMERS = 2
# Heuristic — pending calibration: deliberately overcrowded layers for a
# non-equilibrium regression input, not a recommended stacking clearance.
CROWDED_LAYER_GAP_ANGSTROM = 1.0
CROWDED_LAYER_CASE_IDS = ("cof-5", "tppa-1")


def write_cof_seed(case, path):
    request = cofid_to_build_request(case["cofid"])
    monomers = tuple(
        build_rdkit_monomer(
            f"reactant_{i}", f"reactant_{i}", record.canonical_smiles,
            record.motif_kind, num_conformers=COF_TEST_CONFORMERS,
        )
        for i, record in enumerate(request.monomers, start=1)
    )
    if len(monomers) == 1:
        generator = RingFormingStructureGenerator(RingFormationConfig(topology_id=request.topology))
        candidate = generator.generate(monomers[0], request.template_id)
    else:
        generator = BatchStructureGenerator(BatchGenerationConfig(
            allowed_reactions=(request.template_id,), topology_ids=(request.topology,),
            target_dimensionality=request.target_dimensionality, write_cif=False,
            # Test only the pair search: avoid an additional independent
            # atomistic seed-screening pass before the frozen comparison seed.
            seed_placement_config=SeedPlacementConfig(enabled=False),
        ))
        summary, candidate = generator.generate_monomer_pair_candidate(*monomers, write_cif=False)
        if candidate is None:
            raise ValueError(f"{case['id']}: {summary.status}: {summary.metadata}")
    export = CIFWriter().write_candidate(path, candidate, monomers, cofid=case["cofid"])
    if export.mode != "atomistic":
        raise ValueError(f"{case['id']}: expected atomistic seed, got {export.mode!r}")
    return export


def normalized_report(report):
    return asdict(replace(report, input_path="", output_path=""))


def replicate_cif(source, target, repeats):
    """Exact diagonal P1 supercell retaining every explicit periodic edge."""
    from itertools import product
    from cofkit.periodic_geometry import format_p1_shift, p1_shift
    from cofkit.soft_relax import SoftRelaxConfig, _parse_system
    from cofkit.validation import _instance_id

    block = gemmi.cif.read_file(str(source)).sole_block()
    small = gemmi.make_small_structure_from_block(block)
    system, _ = _parse_system(block, SoftRelaxConfig())
    # Recover each finite precursor's atom-image offsets. A wrapped precursor
    # can straddle the original cell: replica instance labels must follow its
    # bonded atoms, rather than assigning an instance per geometric tile.
    adjacency = {label: [] for label in system.labels}
    for i, j, gain in system.bonds:
        a, b = system.labels[i], system.labels[j]
        if _instance_id(a) == _instance_id(b):
            adjacency[a].append((b, gain))
            adjacency[b].append((a, tuple(-v for v in gain)))
    offsets = {}
    for root in system.labels:
        if root in offsets:
            continue
        offsets[root] = (0, 0, 0)
        pending = [root]
        while pending:
            a = pending.pop()
            for b, gain in adjacency[a]:
                offset = tuple(p + g for p, g in zip(offsets[a], gain))
                if b in offsets:
                    if offsets[b] != offset:
                        raise ValueError("replica instance labels require finite precursor graphs")
                else:
                    offsets[b] = offset
                    pending.append(b)
    cell = small.cell
    nx, ny, nz = repeats
    lines = [
        "data_search_supercell",
        f"_cell_length_a {cell.a * nx:.12g}", f"_cell_length_b {cell.b * ny:.12g}",
        f"_cell_length_c {cell.c * nz:.12g}",
        f"_cell_angle_alpha {cell.alpha:.12g}", f"_cell_angle_beta {cell.beta:.12g}",
        f"_cell_angle_gamma {cell.gamma:.12g}", "loop_",
        "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
        "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy",
    ]
    images = list(product(range(nx), range(ny), range(nz)))

    def label(base, image):
        # Preserve CIFWriter's final _ElementN atom tag so bridge validation
        # continues to group atoms by precursor instance after replication.
        instance, separator, atom = base.rpartition("_")
        prefix = instance if separator else "site"
        tag = atom if separator else base
        precursor_image = tuple((t - offset) % n for t, offset, n in zip(image, offsets[base], repeats))
        return prefix + "_r" + "_".join(str(v) for v in precursor_image) + "_" + tag

    for image in images:
        for site in small.sites:
            f = [(v + t) / n for v, t, n in zip(site.fract, image, repeats)]
            lines.append(f"{label(site.label, image)} {site.element.name} " + " ".join(f"{v:.12g}" for v in f) + " 1")
    tags = (
        "_geom_bond_atom_site_label_1", "_geom_bond_atom_site_label_2",
        "_geom_bond_distance", "_geom_bond_site_symmetry_1",
        "_geom_bond_site_symmetry_2", "_ccdc_geom_bond_type",
    )
    columns = [list(block.find_values(tag)) for tag in tags]
    lines.extend(("loop_", *tags))
    for image in images:
        for a, b, distance, first_sym, second_sym, order in zip(*columns):
            first_shift, second_shift = p1_shift(first_sym), p1_shift(second_sym)
            relative = tuple(t + q - p for t, p, q in zip(image, first_shift, second_shift))
            endpoint = tuple(t % n for t, n in zip(relative, repeats))
            gain = tuple(t // n for t, n in zip(relative, repeats))
            lines.append(f"{label(a, image)} {label(b, endpoint)} {distance} . {format_p1_shift(gain)} {order}")
    Path(target).parent.mkdir(parents=True, exist_ok=True)
    Path(target).write_text("\n".join(lines) + "\n")


def write_crowded_layers(source, target):
    """Two intact copies of a layer separated by an intentional clash gap."""
    original = gemmi.make_small_structure_from_block(gemmi.cif.read_file(str(source)).sole_block())
    replicate_cif(source, target, (1, 1, 2))
    document = gemmi.cif.read_file(str(target))
    block = document.sole_block()
    labels = block.find_values("_atom_site_label")
    z = block.find_values("_atom_site_fract_z")
    fractional_shift = (original.cell.c - CROWDED_LAYER_GAP_ANGSTROM) / (2 * original.cell.c)
    for index, label in enumerate(labels):
        if str(label).rsplit("_", 1)[0].endswith("_r0_0_1"):
            z[index] = f"{float(z[index]) - fractional_shift:.12g}"
    document.write_file(str(target))
