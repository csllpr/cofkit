import sys
import unittest
from math import cos, pi, sin
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import (
    AssemblyState,
    CIFWriter,
    COFEngine,
    COFProject,
    Frame,
    MonomerSpec,
    MotifRef,
    Pose,
    ReactiveMotif,
    Candidate,
    ReactionEvent,
)


def trigonal_motifs(prefix: str, kind: str, radius: float) -> tuple[ReactiveMotif, ...]:
    motifs = []
    for idx, angle in enumerate((0.0, 2.0 * pi / 3.0, -2.0 * pi / 3.0), start=1):
        origin = (radius * cos(angle), radius * sin(angle), 0.0)
        primary = (cos(angle), sin(angle), 0.0)
        motifs.append(
            ReactiveMotif(
                id=f"{prefix}{idx}",
                kind=kind,
                atom_ids=(idx,),
                frame=Frame(origin=origin, primary=primary, normal=(0.0, 0.0, 1.0)),
            )
        )
    return tuple(motifs)


def build_imine_project() -> tuple[COFProject, tuple[MonomerSpec, MonomerSpec]]:
    tri_amine = MonomerSpec(
        id="tapb",
        name="TAPB-like triamine",
        motifs=trigonal_motifs("n", "amine", radius=4.5),
    )
    tri_aldehyde = MonomerSpec(
        id="tfp",
        name="TFP-like trialdehyde",
        motifs=trigonal_motifs("c", "aldehyde", radius=2.4),
    )
    project = COFProject(
        monomers=(tri_amine, tri_aldehyde),
        allowed_reactions=("imine_bridge",),
        target_dimensionality="2D",
        target_topologies=("hcb",),
    )
    return project, (tri_amine, tri_aldehyde)


class CIFWriterTests(unittest.TestCase):
    def test_engine_candidate_exports_legal_coarse_cif(self):
        project, monomers = build_imine_project()
        candidate = COFEngine().run(project).top(1)[0]

        result = CIFWriter().export_candidate(candidate, monomers, data_name="mock imine")

        self.assertEqual(candidate.metadata["net_plan"]["topology"], "hcb")
        self.assertEqual(candidate.metadata["embedding"]["placement_mode"], "single-node-bipartite")
        self.assertEqual(result.mode, "coarse")
        self.assertEqual(result.data_name, "mock_imine")
        self.assertEqual(result.n_sites, 8)
        self.assertIn("data_mock_imine", result.text)
        self.assertIn("_space_group_name_H-M_alt 'P 1'", result.text)
        self.assertIn("_atom_site_label", result.text)

        atom_lines = [
            line for line in result.text.splitlines()
            if line.startswith("m") and len(line.split()) == 6
        ]
        self.assertEqual(len(atom_lines), 8)
        for line in atom_lines:
            parts = line.split()
            for value in parts[2:5]:
                coord = float(value)
                self.assertGreaterEqual(coord, 0.0)
                self.assertLess(coord, 1.0)

    def test_atomistic_export_uses_monomer_atom_positions(self):
        monomer = MonomerSpec(
            id="mono",
            name="atomistic monomer",
            motifs=(
                ReactiveMotif(
                    id="n1",
                    kind="amine",
                    atom_ids=(1,),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                ),
            ),
            atom_symbols=("C", "N"),
            atom_positions=((0.0, 0.0, 0.0), (1.5, 0.0, 0.0)),
            bonds=((0, 1, 1.0),),
        )
        candidate = Candidate(
            id="atomistic-demo",
            score=0.0,
            state=AssemblyState(
                cell=((10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)),
                monomer_poses={"m1": Pose(translation=(11.0, -1.0, 5.0))},
                stacking_state="disabled",
            ),
            events=(),
            metadata={"instance_to_monomer": {"m1": "mono"}},
        )

        result = CIFWriter().export_candidate(candidate, {"mono": monomer})

        self.assertEqual(result.mode, "atomistic")
        self.assertEqual(result.n_sites, 2)
        self.assertIn("m1_C1 C 0.100000 0.900000 0.500000 1.00", result.text)
        self.assertIn("m1_N2 N 0.250000 0.900000 0.500000 1.00", result.text)
        self.assertIn("_geom_bond_atom_site_label_1", result.text)
        self.assertIn("m1_C1 m1_N2 . . 1.500000", result.text)

    def test_keto_enamine_export_shortens_tautomerized_carbonyl_bond_and_removes_phenolic_hydrogen(self):
        amine = MonomerSpec(
            id="amine",
            name="minimal amine",
            motifs=(
                ReactiveMotif(
                    id="ami1",
                    kind="amine",
                    atom_ids=(0, 1, 2, 3),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("keto_enamine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 1},
                ),
            ),
            atom_symbols=("N", "C", "H", "H"),
            atom_positions=((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.1, 1.0, 0.0), (0.1, -1.0, 0.0)),
            bonds=((0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0)),
        )
        keto_aldehyde = MonomerSpec(
            id="keto_aldehyde",
            name="minimal keto aldehyde",
            motifs=(
                ReactiveMotif(
                    id="kal1",
                    kind="keto_aldehyde",
                    atom_ids=(0, 1, 2, 3, 4, 5),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("keto_enamine_bridge",),
                    metadata={
                        "reactive_atom_id": 0,
                        "anchor_atom_id": 2,
                        "aldehyde_oxygen_atom_id": 1,
                        "ortho_hydroxyl_oxygen_atom_id": 3,
                        "ortho_hydroxyl_hydrogen_atom_id": 4,
                        "ortho_hydroxyl_anchor_atom_id": 2,
                    },
                ),
            ),
            atom_symbols=("C", "O", "C", "O", "H", "H"),
            atom_positions=(
                (0.0, 0.0, 0.0),
                (0.0, 1.2, 0.0),
                (1.2, 0.0, 0.0),
                (1.2, 1.5, 0.0),
                (1.2, 2.4, 0.0),
                (-0.8, 0.0, 0.0),
            ),
            bonds=((0, 1, 2.0), (0, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0), (0, 5, 1.0)),
        )
        candidate = Candidate(
            id="keto-enamine-demo",
            score=0.0,
            state=AssemblyState(
                cell=((20.0, 0.0, 0.0), (0.0, 20.0, 0.0), (0.0, 0.0, 10.0)),
                monomer_poses={
                    "m1": Pose(translation=(0.0, 0.0, 0.0)),
                    "m2": Pose(translation=(1.36, 0.0, 0.0)),
                },
                stacking_state="disabled",
            ),
            events=(
                ReactionEvent(
                    id="rxn1",
                    template_id="keto_enamine_bridge",
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id="amine", motif_id="ami1"),
                        MotifRef(monomer_instance_id="m2", monomer_id="keto_aldehyde", motif_id="kal1"),
                    ),
                ),
            ),
            metadata={"instance_to_monomer": {"m1": "amine", "m2": "keto_aldehyde"}},
        )

        result = CIFWriter().export_candidate(candidate, {"amine": amine, "keto_aldehyde": keto_aldehyde})

        realization = result.metadata["reaction_realization"]
        self.assertEqual(realization["removed_atom_symbols"], {"H": 2, "O": 1})
        self.assertNotIn("m2_H5", result.text)
        bond_lines = [line for line in result.text.splitlines() if line.startswith("m") and len(line.split()) == 5]
        carbonyl_bonds = [line for line in bond_lines if line.startswith("m2_C3 m2_O4") or line.startswith("m2_O4 m2_C3")]
        self.assertEqual(len(carbonyl_bonds), 1)
        self.assertAlmostEqual(float(carbonyl_bonds[0].split()[-1]), 1.24, places=2)

    def test_hydrazone_export_removes_both_terminal_hydrazide_hydrogens(self):
        hydrazide = MonomerSpec(
            id="hydrazide",
            name="minimal hydrazide",
            motifs=(
                ReactiveMotif(
                    id="hdz1",
                    kind="hydrazide",
                    atom_ids=(0, 1, 2, 3, 4, 5, 6),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("hydrazone_bridge",),
                    metadata={
                        "reactive_atom_id": 0,
                        "anchor_atom_id": 1,
                        "hydrogen_atom_ids": (5, 6),
                        "internal_nitrogen_atom_id": 1,
                        "carbonyl_carbon_atom_id": 2,
                        "carbonyl_oxygen_atom_id": 3,
                    },
                ),
            ),
            atom_symbols=("N", "N", "C", "O", "C", "H", "H"),
            atom_positions=(
                (0.0, 0.0, 0.0),
                (-1.1, 0.0, 0.0),
                (-2.3, 0.0, 0.0),
                (-3.4, 0.0, 0.0),
                (-2.3, 1.2, 0.0),
                (0.2, 1.0, 0.0),
                (0.2, -1.0, 0.0),
            ),
            bonds=((0, 1, 1.0), (0, 5, 1.0), (0, 6, 1.0), (1, 2, 1.0), (2, 3, 2.0), (2, 4, 1.0)),
        )
        aldehyde = MonomerSpec(
            id="aldehyde",
            name="minimal aldehyde",
            motifs=(
                ReactiveMotif(
                    id="ald1",
                    kind="aldehyde",
                    atom_ids=(0, 1, 2, 3),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("hydrazone_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 2},
                ),
            ),
            atom_symbols=("C", "O", "C", "H"),
            atom_positions=((0.0, 0.0, 0.0), (0.0, 1.2, 0.0), (1.2, 0.0, 0.0), (-0.8, 0.0, 0.0)),
            bonds=((0, 1, 2.0), (0, 2, 1.0), (0, 3, 1.0)),
        )
        candidate = Candidate(
            id="hydrazone-demo",
            score=0.0,
            state=AssemblyState(
                cell=((20.0, 0.0, 0.0), (0.0, 20.0, 0.0), (0.0, 0.0, 10.0)),
                monomer_poses={
                    "m1": Pose(translation=(0.0, 0.0, 0.0)),
                    "m2": Pose(translation=(1.3, 0.0, 0.0)),
                },
                stacking_state="disabled",
            ),
            events=(
                ReactionEvent(
                    id="rxn1",
                    template_id="hydrazone_bridge",
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id="hydrazide", motif_id="hdz1"),
                        MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
                    ),
                ),
            ),
            metadata={"instance_to_monomer": {"m1": "hydrazide", "m2": "aldehyde"}},
        )

        result = CIFWriter().export_candidate(candidate, {"hydrazide": hydrazide, "aldehyde": aldehyde})

        realization = result.metadata["reaction_realization"]
        self.assertEqual(realization["removed_atom_symbols"], {"H": 2, "O": 1})
        self.assertNotIn("m1_H6", result.text)
        self.assertNotIn("m1_H7", result.text)
        self.assertIn("m2_C1 m1_N1 . . 1.300000", result.text)


if __name__ == "__main__":
    unittest.main()
