"""Explicit periodic edges must survive changes of atom representatives."""

from collections import deque
from dataclasses import replace

import gemmi
import networkx as nx
import pytest

from cofkit.cif import CIFWriter
from cofkit.batch import BatchGenerationConfig, BatchStructureGenerator
from cofkit.chem.rdkit import build_rdkit_monomer
from cofkit.model import AssemblyState, Candidate, MonomerSpec, Pose
from cofkit.periodic_geometry import format_p1_shift, p1_shift
from cofkit.reaction_realization import (
    ReactionRealizationResult,
    RealizedAtom,
    RealizedBond,
)
from cofkit.soft_relax import (
    SoftRelaxConfig,
    _rewrite_cif_coordinates,
    relax_cif_clashes,
)


def _edges(block):
    rows = block.find(
        [
            "_geom_bond_atom_site_label_1",
            "_geom_bond_atom_site_label_2",
            "_geom_bond_site_symmetry_1",
            "_geom_bond_site_symmetry_2",
            "_geom_bond_distance",
        ]
    )
    return [
        (
            str(a),
            str(b),
            tuple(y - x for x, y in zip(p1_shift(str(s1)), p1_shift(str(s2)))),
            float(d),
        )
        for a, b, s1, s2, d in rows
    ]


def _ball(edges, radius=6):
    """Independent exact finite neighborhood in the infinite periodic lift."""
    adjacency = {}
    for a, b, gain, _ in edges:
        adjacency.setdefault(a, []).append((b, gain))
        adjacency.setdefault(b, []).append((a, tuple(-v for v in gain)))
    root = ("A", (0, 0, 0))
    depths = {root: 0}
    queue = deque([root])
    while queue:
        atom, image = queue.popleft()
        if depths[atom, image] == radius:
            continue
        for neighbor, gain in adjacency[atom]:
            state = (neighbor, tuple(x + y for x, y in zip(image, gain)))
            if state not in depths:
                depths[state] = depths[atom, image] + 1
                queue.append(state)
    graph = nx.MultiGraph()
    for state in depths:
        graph.add_node(state, color=(state[0], state == root))
    for atom, image in depths:
        for neighbor, gain in adjacency[atom]:
            state = (neighbor, tuple(x + y for x, y in zip(image, gain)))
            if state in depths and (atom, image) < state:
                graph.add_edge((atom, image), state)
    return graph


def _chain_cif(path, *, gamma=90, replicas=1):
    atoms = []
    bonds = []
    for i in range(replicas):
        a = "A" if i == 0 else f"A{i}"
        b = "B" if i == 0 else f"B{i}"
        next_a = f"A{i + 1}" if i + 1 < replicas else "A"
        atoms.extend(
            [
                f"{a} C {(i + 0.25) / replicas} .5 .5",
                f"{b} C {(i + 0.75) / replicas} .5 .5",
            ]
        )
        bonds.extend(
            [
                f"{a} {b} . . 1.4 D",
                f"{b} {next_a} 1_655 1_655 1.4 D"
                if i + 1 < replicas
                else f"{b} {next_a} . 1_655 1.4 D",
            ]
        )
    path.write_text(
        f"""data_chain
_cell_length_a {2.8 * replicas}
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma {gamma}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
"""
        + "\n".join(atoms)
        + """
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_site_symmetry_1
_geom_bond_site_symmetry_2
_geom_bond_distance
_ccdc_geom_bond_type
"""
        + "\n".join(bonds)
        + "\n"
    )


@pytest.mark.parametrize("gamma", [90, 25])
@pytest.mark.parametrize("replicas", [1, 3])
@pytest.mark.parametrize("move", [False, True])
def test_coordinate_rewrite_preserves_lifted_chain(tmp_path, gamma, replicas, move):
    source, target = tmp_path / "chain.cif", tmp_path / "rewritten.cif"
    _chain_cif(source, gamma=gamma, replicas=replicas)
    block = gemmi.cif.read_file(str(source)).sole_block()
    small = gemmi.make_small_structure_from_block(block)
    labels = [s.label for s in small.sites]
    positions = [[s.fract.x, s.fract.y, s.fract.z] for s in small.sites]
    if move:
        # Different endpoints traverse different cell boundaries; motions
        # need not be chemically plausible to test connectivity bookkeeping.
        for i, frac in enumerate(positions):
            frac[0] += 1.1 if i % 2 else -0.8
            frac[1] += -1.2 if i % 2 else 0.6
    _rewrite_cif_coordinates(source, target, positions, labels, small.cell)
    after = gemmi.cif.read_file(str(target)).sole_block()
    before_edges, after_edges = _edges(block), _edges(after)
    assert len(after_edges) == len(before_edges)
    assert nx.is_isomorphic(
        _ball(before_edges),
        _ball(after_edges),
        node_match=lambda a, b: a["color"] == b["color"],
    )
    if not move:
        assert [d for *_, d in after_edges] == pytest.approx([1.4] * (2 * replicas))


def test_soft_relax_keeps_periodic_multiedges(tmp_path):
    source, target = tmp_path / "chain.cif", tmp_path / "relaxed.cif"
    _chain_cif(source)
    report = relax_cif_clashes(
        source, target, config=SoftRelaxConfig(repulsion_coefficients=(0.0,))
    )
    assert report.converged
    edges = _edges(gemmi.cif.read_file(str(target)).sole_block())
    assert len(_ball(edges)) == 13  # six bonds in both directions from the root
    assert report.max_bond_drift == pytest.approx(0)


def test_coordinate_rewrite_materializes_implicit_identity_images(tmp_path):
    source, target = tmp_path / "molecule.cif", tmp_path / "wrapped.cif"
    _chain_cif(source)
    # Replace the periodic chain by one explicit identity edge, omitting
    # optional symmetry columns as ordinary molecular CIFs may do.
    source.write_text(
        source.read_text().rpartition("loop_")[0]
        + """loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_ccdc_geom_bond_type
A B 1.4 D
"""
    )
    _rewrite_cif_coordinates(
        source,
        target,
        [[0.9, 0.5, 0.5], [1.4, 0.5, 0.5]],
        ["A", "B"],
        gemmi.UnitCell(2.8, 10, 10, 90, 90, 90),
    )
    edges = _edges(gemmi.cif.read_file(str(target)).sole_block())
    assert len(edges) == 1
    assert edges[0][2] == (1, 0, 0)
    assert edges[0][3] == pytest.approx(1.4)


def _candidate(cell, poses):
    return Candidate(
        id="periodic-test",
        score=0,
        state=AssemblyState(cell=cell, monomer_poses=poses, stacking_state="disabled"),
        events=(),
        metadata={"instance_to_monomer": {name: name for name in poses}},
    )


def test_cif_writer_keeps_intended_intramonomer_bond_in_small_cell():
    monomer = MonomerSpec(
        id="m",
        name="m",
        motifs=(),
        atom_symbols=("C", "C"),
        atom_positions=((0, 0, 0), (1.4, 0, 0)),
        bonds=((0, 1, 1.0),),
    )
    candidate = _candidate(((2, 0, 0), (0, 10, 0), (0, 0, 10)), {"m": Pose()})
    result = CIFWriter().export_candidate(candidate, {"m": monomer})
    edges = _edges(gemmi.cif.read_string(result.text).sole_block())
    assert len(edges) == 1
    assert edges[0][2] == (0, 0, 0)
    assert edges[0][3] == pytest.approx(1.4)  # nearest image would be only 0.6 A


def test_cif_writer_keeps_parallel_reaction_images_and_reverse_override():
    specs = {
        name: MonomerSpec(
            id=name,
            name=name,
            motifs=(),
            atom_symbols=("C",),
            atom_positions=((0, 0, 0),),
        )
        for name in ["a", "b"]
    }
    candidate = _candidate(
        ((2.8, 0, 0), (0, 10, 0), (0, 0, 10)),
        {"a": Pose(translation=(0.7, 0, 0)), "b": Pose(translation=(2.1, 0, 0))},
    )
    result = ReactionRealizationResult(
        atoms_by_instance={
            name: (RealizedAtom(0, name.upper(), "C", (0, 0, 0)),) for name in specs
        },
        bonds=(
            RealizedBond("A", "B", 1.4),
            RealizedBond("A", "B", 1.4, symmetry_2="1_455"),
            RealizedBond("B", "A", 1.4, symmetry_2="1_655", bond_order=2),
        ),
    )

    class Realizer:
        def realize(self, *args):
            return result

    writer = CIFWriter(reaction_realizer=Realizer())
    exported = writer.export_candidate(candidate, specs)
    block = gemmi.cif.read_string(exported.text).sole_block()
    edges = _edges(block)
    assert len(edges) == 2
    assert len(_ball(edges)) == 13
    assert sorted(block.find_loop("_ccdc_geom_bond_type")) == ["D", "S"]

    # Atom representatives and cell basis are changed together; the two
    # physical reaction bonds and their lengths must stay intact.
    # a' = a + b, b' = b transforms each old gain (x,y,z) to
    # (x,y-x,z). Atom Cartesian coordinates are unchanged.
    sheared = replace(
        candidate,
        state=replace(candidate.state, cell=((2.8, 10, 0), (0, 10, 0), (0, 0, 10))),
    )
    result = replace(
        result,
        bonds=(
            result.bonds[0],
            replace(result.bonds[1], symmetry_2="1_465"),
            replace(result.bonds[2], symmetry_2="1_645"),
        ),
    )
    sheared_edges = _edges(
        gemmi.cif.read_string(writer.export_candidate(sheared, specs).text).sole_block()
    )
    assert nx.is_isomorphic(
        _ball(edges),
        _ball(sheared_edges),
        node_match=lambda a, b: a["color"] == b["color"],
    )
    assert [e[3] for e in sheared_edges] == pytest.approx([1.4, 1.4])


def test_p1_formatter_rejects_unrepresentable_gain():
    assert p1_shift(format_p1_shift((-5, 4, 0))) == (-5, 4, 0)
    with pytest.raises(ValueError, match="formatter range"):
        format_p1_shift((5, 0, 0))


def test_cif_writer_preserves_bond_to_self_image():
    spec = MonomerSpec(
        id="a", name="a", motifs=(), atom_symbols=("C",), atom_positions=((0, 0, 0),)
    )
    candidate = _candidate(((2.8, 0, 0), (0, 10, 0), (0, 0, 10)), {"a": Pose()})

    class Realizer:
        def realize(self, *args):
            return ReactionRealizationResult(
                atoms_by_instance={"a": (RealizedAtom(0, "A", "C", (0, 0, 0)),)},
                bonds=(RealizedBond("A", "A", 2.8, symmetry_2="1_655"),),
            )

    text = (
        CIFWriter(reaction_realizer=Realizer())
        .export_candidate(candidate, {"a": spec})
        .text
    )
    assert _edges(gemmi.cif.read_string(text).sole_block()) == [
        ("A", "A", (1, 0, 0), 2.8)
    ]


def test_common_event_translation_does_not_change_realized_connectivity():
    specs = (
        build_rdkit_monomer("a", "a", "Nc1cc(N)cc(N)c1", "amine", num_conformers=1),
        build_rdkit_monomer(
            "b", "b", "O=Cc1cc(C=O)cc(C=O)c1", "aldehyde", num_conformers=1
        ),
    )
    generator = BatchStructureGenerator(
        BatchGenerationConfig(topology_ids=("hcb",), max_workers=1)
    )
    _, candidate = generator.generate_monomer_pair_candidate(*specs, write_cif=False)
    assert candidate is not None
    translated = replace(
        candidate,
        events=tuple(
            replace(
                event,
                participants=tuple(
                    replace(
                        ref,
                        periodic_image=tuple(
                            x + y for x, y in zip(ref.periodic_image, (8, -7, 6))
                        ),
                    )
                    for ref in event.participants
                ),
            )
            for event in candidate.events
        ),
    )
    writer = CIFWriter()
    original = _edges(
        gemmi.cif.read_string(
            writer.export_candidate(candidate, specs).text
        ).sole_block()
    )
    shifted = _edges(
        gemmi.cif.read_string(
            writer.export_candidate(translated, specs).text
        ).sole_block()
    )
    assert [row[:3] for row in original] == [row[:3] for row in shifted]
    assert [row[3] for row in original] == pytest.approx(
        [row[3] for row in shifted], abs=1e-6
    )
