import math
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


# Curated graph-family labels for the four-connecting imine precursors
# shipped in the default example library.  These labels come from
# classify_monomer_graph_shape, a graph-automorphism diagnostic; spatial
# topology filtering uses the geometric classifier instead.
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


class GeometricNodeShapeTests(unittest.TestCase):
    """Connector-shape classification from conformer geometry."""

    def test_square_planar_connectors_classify_square(self):
        self.assertEqual(classify_monomer_node_shape(_monomer("m", SQUARE_ORIGINS)).label, SHAPE_SQUARE)

    def test_rectangular_connectors_classify_rectangular(self):
        self.assertEqual(classify_monomer_node_shape(_monomer("m", RECTANGULAR_ORIGINS)).label, SHAPE_RECTANGULAR)

    def test_tetrahedral_connectors_classify_tetrahedral(self):
        self.assertEqual(classify_monomer_node_shape(_monomer("m", TETRAHEDRAL_ORIGINS)).label, SHAPE_TETRAHEDRAL)

    def test_classification_is_rigid_motion_invariant(self):
        angle = 0.7
        ca, sa = math.cos(angle), math.sin(angle)

        def rotate(point):
            x, y, z = point
            # arbitrary proper rotation: Rz(0.7) then Rx-like tilt
            x, y = ca * x - sa * y, sa * x + ca * y
            y, z = 0.6 * y - 0.8 * z, 0.8 * y + 0.6 * z
            return (x + 5.0, y - 3.0, z + 1.0)

        for origins, expected in (
            (SQUARE_ORIGINS, SHAPE_SQUARE),
            (RECTANGULAR_ORIGINS, SHAPE_RECTANGULAR),
            (TETRAHEDRAL_ORIGINS, SHAPE_TETRAHEDRAL),
        ):
            with self.subTest(expected=expected):
                moved = tuple(rotate(origin) for origin in origins)
                self.assertEqual(classify_monomer_node_shape(_monomer("m", moved)).label, expected)

    def test_irregular_connectors_are_unknown(self):
        self.assertEqual(classify_monomer_node_shape(_monomer("m", IRREGULAR_ORIGINS)).label, SHAPE_UNKNOWN)

    def test_ditopic_monomer_is_unknown(self):
        self.assertEqual(classify_monomer_node_shape(_monomer("m", LINKER_ORIGINS)).label, SHAPE_UNKNOWN)

    def test_missing_or_degenerate_evidence_is_unknown(self):
        # No atom coordinates at all (bare motif frames at the origin).
        degenerate = _monomer("m", ((0.0, 0.0, 0.0),) * 4)
        self.assertEqual(classify_monomer_node_shape(degenerate).label, SHAPE_UNKNOWN)
        # Non-finite coordinates carry no evidence.
        nan_monomer = _monomer("m", ((float("nan"), 0.0, 0.0),) + SQUARE_ORIGINS[1:])
        self.assertEqual(classify_monomer_node_shape(nan_monomer).label, SHAPE_UNKNOWN)

    def test_real_d2h_and_td_monomers_classify_spatially(self):
        self.assertEqual(
            classify_monomer_node_shape(build_rdkit_monomer("a", "a", D2H_TETRA_AMINE, "amine")).label,
            SHAPE_RECTANGULAR,
        )
        self.assertEqual(
            classify_monomer_node_shape(
                build_rdkit_monomer(
                    "b",
                    "b",
                    "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1",
                    "amine",
                )
            ).label,
            SHAPE_TETRAHEDRAL,
        )


class ShapeFilteredEnumerationTests(unittest.TestCase):
    """Enumerated topology pools react to the geometric node shape."""

    @staticmethod
    def _generator(**overrides):
        config = BatchGenerationConfig(rdkit_num_conformers=1, **overrides)
        return BatchStructureGenerator(config)

    def test_square_node_enumerates_sql(self):
        shape = classify_monomer_node_shape(_monomer("sq", SQUARE_ORIGINS))
        self.assertEqual(shape.label, SHAPE_SQUARE)
        ids = self._generator()._topology_ids_for_pair(
            connectivities=(4, 2), pair_mode="node-linker", node_shapes=(shape,)
        )
        self.assertIn("sql", ids)
        self.assertNotIn("kgm", ids)

    def test_rectangular_node_enumerates_kgm_not_sql(self):
        shape = classify_monomer_node_shape(_monomer("re", RECTANGULAR_ORIGINS))
        ids = self._generator()._topology_ids_for_pair(
            connectivities=(4, 2), pair_mode="node-linker", node_shapes=(shape,)
        )
        self.assertIn("kgm", ids)
        self.assertNotIn("sql", ids)

    def test_tetrahedral_node_enumerates_dia_for_3d_target(self):
        shape = classify_monomer_node_shape(_monomer("td", TETRAHEDRAL_ORIGINS))
        ids_3d = self._generator(target_dimensionality="3D")._topology_ids_for_pair(
            connectivities=(4, 4), pair_mode="node-node", node_shapes=(shape, shape)
        )
        self.assertIn("dia", ids_3d)
        # A 2D target offers no confident-shape 2D family to a tetrahedral node.
        ids_2d = self._generator()._topology_ids_for_pair(
            connectivities=(4, 4), pair_mode="node-node", node_shapes=(shape, shape)
        )
        self.assertNotIn("sql", ids_2d)

    def test_unknown_shape_excludes_nothing(self):
        shape = classify_monomer_node_shape(_monomer("un", IRREGULAR_ORIGINS))
        self.assertEqual(shape.label, SHAPE_UNKNOWN)
        ids = self._generator()._topology_ids_for_pair(
            connectivities=(4, 2), pair_mode="node-linker", node_shapes=(shape,)
        )
        self.assertIn("sql", ids)
        self.assertIn("kgm", ids)
        ids_3d = self._generator(target_dimensionality="3D")._topology_ids_for_pair(
            connectivities=(4, 4), pair_mode="node-node", node_shapes=(shape, shape)
        )
        self.assertIn("dia", ids_3d)


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


class ExplicitDimensionalityConflictTests(unittest.TestCase):
    """An explicit topology of the wrong dimensionality is unsatisfiable.

    Unlike heuristic shape conflicts (warning + build anyway), a
    dimensionality conflict is exact metadata: the request cannot be
    satisfied, so planner, engine, and batch generator reject it with a
    clear error instead of silently building another dimensionality or
    falling back to a topology-free plan.
    """

    def test_planner_raises_on_dimensionality_conflict(self):
        with self.assertRaises(ValueError) as raised:
            NetPlanner().propose(
                _d2h_tetratopic_monomers(),
                (),
                "2D",
                target_topologies=("dia",),
            )
        self.assertIn("'dia' is 3D", str(raised.exception))

    def test_planner_accepts_consistent_explicit_request(self):
        plans = NetPlanner().propose(
            _d2h_tetratopic_monomers(),
            (),
            "2D",
            target_topologies=("sql",),
        )
        self.assertEqual([plan.topology.id for plan in plans], ["sql"])

    def test_engine_raises_on_dimensionality_conflict(self):
        project = COFProject(
            monomers=_d2h_tetratopic_monomers(),
            allowed_reactions=("imine_bridge",),
            target_dimensionality="2D",
            target_topologies=("dia",),
        )
        with self.assertRaises(ValueError) as raised:
            COFEngine().run(project)
        self.assertIn("'dia' is 3D", str(raised.exception))

    def test_batch_generator_raises_on_dimensionality_conflict(self):
        with self.assertRaises(ValueError) as raised:
            BatchStructureGenerator(
                BatchGenerationConfig(target_dimensionality="2D", topology_ids=("dia",))
            )
        self.assertIn("'dia' is 3D", str(raised.exception))
        with self.assertRaises(ValueError):
            BatchStructureGenerator(
                BatchGenerationConfig(target_dimensionality="3D", single_node_topology_ids=("sql",))
            )

    def test_batch_config_rejects_unknown_dimensionality(self):
        with self.assertRaises(ValueError):
            BatchGenerationConfig(target_dimensionality="4D")

    def test_unknown_explicit_topology_id_is_not_rejected_at_construction(self):
        # Unresolvable ids keep the existing per-pair KeyError reporting.
        BatchStructureGenerator(
            BatchGenerationConfig(target_dimensionality="2D", topology_ids=("no-such-net",))
        )


if __name__ == "__main__":
    unittest.main()
