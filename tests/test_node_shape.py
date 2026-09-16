import unittest

from cofkit import (
    BatchGenerationConfig,
    BatchStructureGenerator,
    COFEngine,
    COFProject,
    Frame,
    MonomerSpec,
    ReactiveMotif,
)
from cofkit.node_shape import (
    SHAPE_RECTANGULAR,
    SHAPE_SQUARE,
    SHAPE_TETRAHEDRAL,
    SHAPE_UNKNOWN,
    classify_monomer_node_shape,
    classify_topology_node_shape,
    shapes_compatible,
)
from cofkit.planner import NetPlanner, TopologyHint
from cofkit import build_rdkit_monomer


# Curated labels for the four-connecting imine precursors shipped in the
# default example library.  These labels describe the symmetry family used
# for topology selection (C4/Td versus D2h), not a claim about the relaxed
# 3-D conformation of the isolated molecule.
DEFAULT_TETRATOPIC_IMINE_SHAPE_CASES = (
    ("tetraphenylmethane_tetraamine", "amine", "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1", "c4_td"),
    ("tetraphenylmethane_tetrabenzaldehyde", "aldehyde", "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1", "c4_td"),
    ("1,2,4,5_tetraaminobenzene", "amine", "Nc1cc(N)c(N)cc1N", "d2h"),
    ("1,2,4,5_tetrabenzaldehyde", "aldehyde", "O=Cc1cc(C=O)c(C=O)cc1C=O", "d2h"),
    # Inert methyl/fluoro substitution must not change the reactive-site
    # symmetry family used by the classifier.
    ("substituted_tetraphenylmethane_tetraamine", "amine", "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1C", "c4_td"),
    ("substituted_tetraphenylmethane_tetrabenzaldehyde", "aldehyde", "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1F", "c4_td"),
    ("substituted_d2h_tetraamine", "amine", "Nc1cc(N)c(N)c(F)c1N", "d2h"),
    ("substituted_d2h_tetrabenzaldehyde", "aldehyde", "O=Cc1cc(C=O)c(C=O)c(F)c1C=O", "d2h"),
)


class DefaultMonomerShapeLabelTests(unittest.TestCase):
    def test_labeled_default_tetratopic_imine_set_is_parseable(self):
        for monomer_id, kind, smiles, expected_family in DEFAULT_TETRATOPIC_IMINE_SHAPE_CASES:
            with self.subTest(monomer_id=monomer_id):
                spec = build_rdkit_monomer(monomer_id, monomer_id, smiles, kind)
                self.assertEqual(len(spec.motifs), 4)
                self.assertIn(expected_family, {"c4_td", "d2h"})
                from cofkit.node_shape import classify_monomer_graph_shape
                self.assertEqual(classify_monomer_graph_shape(spec), expected_family)


def _motif(motif_id, origin):
    return ReactiveMotif(
        id=motif_id,
        kind="amine",
        atom_ids=(0,),
        frame=Frame(origin=origin, primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
    )


def _monomer(monomer_id, origins):
    return MonomerSpec(
        id=monomer_id,
        name=monomer_id,
        motifs=tuple(_motif(f"m{index}", origin) for index, origin in enumerate(origins)),
    )


SQUARE_ORIGINS = ((2.0, 0.0, 0.0), (0.0, 2.0, 0.0), (-2.0, 0.0, 0.0), (0.0, -2.0, 0.0))
RECTANGULAR_ORIGINS = ((1.732, 1.0, 0.0), (-1.732, 1.0, 0.0), (-1.732, -1.0, 0.0), (1.732, -1.0, 0.0))
TETRAHEDRAL_ORIGINS = ((1.0, 1.0, 1.0), (1.0, -1.0, -1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0))
IRREGULAR_ORIGINS = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (-1.2, 0.1, 0.0), (0.1, -1.1, 0.0))
LINKER_ORIGINS = ((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0))


class TopologyShapeTests(unittest.TestCase):
    def test_supported_layouts(self):
        self.assertEqual(classify_topology_node_shape("sql").label, SHAPE_SQUARE)
        self.assertEqual(classify_topology_node_shape("kgm").label, SHAPE_RECTANGULAR)
        self.assertEqual(classify_topology_node_shape("dia").label, SHAPE_TETRAHEDRAL)

    def test_unresolvable_layouts_are_unknown(self):
        for topology_id in ("nbo", "srs", "she", "bex", "pcu", "hxl", "htb"):
            self.assertEqual(
                classify_topology_node_shape(topology_id).label,
                SHAPE_UNKNOWN,
                msg=topology_id,
            )

    def test_unknown_topology_id_is_unknown(self):
        self.assertEqual(classify_topology_node_shape("does-not-exist").label, SHAPE_UNKNOWN)


D2H_TETRA_AMINE = "Nc1cc(N)c(N)cc1N"
D2H_TETRA_ALDEHYDE = "O=Cc1cc(C=O)c(C=O)cc1C=O"


def _d2h_tetratopic_monomers():
    return (
        build_rdkit_monomer("d2h_amine", "d2h_amine", D2H_TETRA_AMINE, "amine"),
        build_rdkit_monomer("d2h_aldehyde", "d2h_aldehyde", D2H_TETRA_ALDEHYDE, "aldehyde"),
    )


class PlannerShapeFilterTests(unittest.TestCase):
    """Explicit topology requests survive a node-shape disagreement as warnings."""

    @classmethod
    def setUpClass(cls):
        cls.monomers = _d2h_tetratopic_monomers()

    def test_explicit_request_records_shape_conflict_instead_of_dropping_it(self):
        plans = NetPlanner().propose(self.monomers, (), "2D", target_topologies=("sql",))

        self.assertEqual([plan.topology.id for plan in plans], ["sql"])
        warnings = plans[0].metadata["shape_warnings"]
        self.assertEqual(len(warnings), 2)
        for warning in warnings:
            self.assertIn("requested topology 'sql'", warning)
            self.assertIn("rectangular", warning)

    def test_shape_filter_toggle_keeps_conflicting_topology_without_warning(self):
        plans = NetPlanner(shape_aware_topology_filter=False).propose(
            self.monomers,
            (),
            "2D",
            target_topologies=("sql",),
        )

        self.assertEqual([plan.topology.id for plan in plans], ["sql"])
        self.assertNotIn("shape_warnings", plans[0].metadata)


class EngineShapeFilterTests(unittest.TestCase):
    """COFEngine keeps explicit topology requests and reports the conflict."""

    def test_explicit_topology_request_reports_shape_conflict_on_candidate(self):
        project = COFProject(
            monomers=_d2h_tetratopic_monomers(),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
            target_topologies=("sql",),
        )

        candidate = COFEngine().run(project).top(1)[0]

        self.assertEqual(candidate.metadata["net_plan"]["topology"], "sql")
        warnings = candidate.metadata["shape_warnings"]
        self.assertEqual(len(warnings), 1)
        self.assertIn("rectangular", warnings[0])


if __name__ == "__main__":
    unittest.main()
