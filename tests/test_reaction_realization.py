import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import AssemblyState, Candidate, Frame, MonomerSpec, Pose, ReactiveMotif, ReactionEvent, MotifRef
from cofkit.reaction_realization import EventRealization, ReactionEventRealizationRegistry, ReactionRealizer


def _single_event_candidate(
    template_id: str,
    first_ref: MotifRef,
    second_ref: MotifRef,
    *,
    distance: float = 1.3,
) -> Candidate:
    return Candidate(
        id=f"{template_id}-demo",
        score=0.0,
        state=AssemblyState(
            cell=((20.0, 0.0, 0.0), (0.0, 20.0, 0.0), (0.0, 0.0, 10.0)),
            monomer_poses={
                first_ref.monomer_instance_id: Pose(translation=(0.0, 0.0, 0.0)),
                second_ref.monomer_instance_id: Pose(translation=(distance, 0.0, 0.0)),
            },
            stacking_state="disabled",
        ),
        events=(
            ReactionEvent(
                id="rxn1",
                template_id=template_id,
                participants=(first_ref, second_ref),
            ),
        ),
    )


class ReactionRealizationTests(unittest.TestCase):
    def _world_position(self, pose: Pose, local_position):
        return ReactionRealizer()._world_position(pose, local_position)

    def _assert_retained_hydrogen_reoriented(
        self,
        *,
        result,
        monomer: MonomerSpec,
        instance_id: str,
        parent_atom_id: int,
        hydrogen_atom_id: int,
    ) -> None:
        realizer = ReactionRealizer()
        realized_atoms = {atom.atom_id: atom for atom in result.atoms_by_instance[instance_id]}
        parent_local = realized_atoms[parent_atom_id].local_position
        hydrogen_local = realized_atoms[hydrogen_atom_id].local_position
        hydrogen_vector = (
            hydrogen_local[0] - parent_local[0],
            hydrogen_local[1] - parent_local[1],
            hydrogen_local[2] - parent_local[2],
        )
        self.assertGreater(realizer._distance(hydrogen_local, monomer.atom_positions[hydrogen_atom_id]), 0.25)
        self.assertAlmostEqual(hydrogen_vector[0], 0.0, places=6)
        self.assertGreater(abs(hydrogen_vector[1]), 0.7)

    def test_imine_realization_removes_oxygen_and_both_amine_hydrogens(self):
        amine = MonomerSpec(
            id="amine",
            name="minimal amine",
            motifs=(
                ReactiveMotif(
                    id="ami1",
                    kind="amine",
                    atom_ids=(0, 1, 2, 3),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("imine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 1},
                ),
            ),
            atom_symbols=("N", "C", "H", "H"),
            atom_positions=((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.1, 1.0, 0.0), (0.1, -1.0, 0.0)),
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
                    allowed_reaction_templates=("imine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 2},
                ),
            ),
            atom_symbols=("C", "O", "C", "H"),
            atom_positions=((0.0, 0.0, 0.0), (0.0, 1.2, 0.0), (1.2, 0.0, 0.0), (-0.8, 0.0, 0.0)),
        )
        candidate = Candidate(
            id="imine-demo",
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
                    template_id="imine_bridge",
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id="amine", motif_id="ami1"),
                        MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
                    ),
                ),
            ),
            metadata={"instance_to_monomer": {"m1": "amine", "m2": "aldehyde"}},
        )

        result = ReactionRealizer().realize(candidate, {"amine": amine, "aldehyde": aldehyde}, {"m1": "amine", "m2": "aldehyde"})

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["applied_event_count"], 1)
        self.assertEqual(result.metadata["removed_atom_symbols"], {"H": 2, "O": 1})
        self.assertEqual([atom.symbol for atom in result.atoms_by_instance["m1"]], ["N", "C"])
        self.assertEqual([atom.symbol for atom in result.atoms_by_instance["m2"]], ["C", "C", "H"])
        self.assertEqual(len(result.bonds), 1)
        self.assertEqual(result.bonds[0].label_1, "m2_C1")
        self.assertEqual(result.bonds[0].label_2, "m1_N1")
        self.assertAlmostEqual(result.bonds[0].distance, 1.3, places=6)

    def test_hydrazone_realization_removes_oxygen_and_both_terminal_hydrazide_hydrogens(self):
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
        candidate = _single_event_candidate(
            "hydrazone_bridge",
            MotifRef(monomer_instance_id="m1", monomer_id="hydrazide", motif_id="hdz1"),
            MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
        )

        result = ReactionRealizer().realize(candidate, {"hydrazide": hydrazide, "aldehyde": aldehyde}, {"m1": "hydrazide", "m2": "aldehyde"})

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["applied_templates"], {"hydrazone_bridge": 1})
        self.assertEqual(result.metadata["removed_atom_symbols"], {"H": 2, "O": 1})
        self.assertEqual([atom.symbol for atom in result.atoms_by_instance["m1"]], ["N", "N", "C", "O", "C"])
        self.assertEqual(len(result.bonds), 1)
        self.assertEqual(result.bonds[0].label_1, "m2_C1")
        self.assertEqual(result.bonds[0].label_2, "m1_N1")
        self.assertAlmostEqual(result.bonds[0].distance, 1.3, places=6)
        self.assertEqual(result.metadata["hydrogen_cleanup"]["atom_labels"], ("m2_H4",))
        self._assert_retained_hydrogen_reoriented(
            result=result,
            monomer=aldehyde,
            instance_id="m2",
            parent_atom_id=0,
            hydrogen_atom_id=3,
        )

    def test_azine_realization_removes_oxygen_and_both_hydrazine_hydrogens(self):
        hydrazine = MonomerSpec(
            id="hydrazine",
            name="minimal hydrazine",
            motifs=(
                ReactiveMotif(
                    id="hyd1",
                    kind="hydrazine",
                    atom_ids=(0, 1, 2, 3),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("azine_bridge",),
                    metadata={
                        "reactive_atom_id": 0,
                        "anchor_atom_id": 1,
                        "hydrogen_atom_ids": (2, 3),
                        "internal_nitrogen_atom_id": 1,
                    },
                ),
            ),
            atom_symbols=("N", "N", "H", "H", "H", "H"),
            atom_positions=(
                (0.0, 0.0, 0.0),
                (-1.1, 0.0, 0.0),
                (0.2, 1.0, 0.0),
                (0.2, -1.0, 0.0),
                (-1.3, 1.0, 0.0),
                (-1.3, -1.0, 0.0),
            ),
            bonds=((0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (1, 4, 1.0), (1, 5, 1.0)),
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
                    allowed_reaction_templates=("azine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 2},
                ),
            ),
            atom_symbols=("C", "O", "C", "H"),
            atom_positions=((0.0, 0.0, 0.0), (0.0, 1.2, 0.0), (1.2, 0.0, 0.0), (-0.8, 0.0, 0.0)),
            bonds=((0, 1, 2.0), (0, 2, 1.0), (0, 3, 1.0)),
        )
        candidate = _single_event_candidate(
            "azine_bridge",
            MotifRef(monomer_instance_id="m1", monomer_id="hydrazine", motif_id="hyd1"),
            MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
        )

        result = ReactionRealizer().realize(candidate, {"hydrazine": hydrazine, "aldehyde": aldehyde}, {"m1": "hydrazine", "m2": "aldehyde"})

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["applied_templates"], {"azine_bridge": 1})
        self.assertEqual(result.metadata["removed_atom_symbols"], {"H": 2, "O": 1})
        self.assertEqual([atom.symbol for atom in result.atoms_by_instance["m1"]], ["N", "N", "H", "H"])
        self.assertEqual(len(result.bonds), 1)
        self.assertEqual(result.bonds[0].label_1, "m2_C1")
        self.assertEqual(result.bonds[0].label_2, "m1_N1")
        self.assertAlmostEqual(result.bonds[0].distance, 1.3, places=6)
        self.assertEqual(result.metadata["hydrogen_cleanup"]["atom_labels"], ("m2_H4",))
        self._assert_retained_hydrogen_reoriented(
            result=result,
            monomer=aldehyde,
            instance_id="m2",
            parent_atom_id=0,
            hydrogen_atom_id=3,
        )

    def test_azine_realization_bends_shared_hydrazine_bridge_toward_120_degrees(self):
        hydrazine = MonomerSpec(
            id="hydrazine",
            name="linked hydrazine",
            motifs=(
                ReactiveMotif(
                    id="hyd_left",
                    kind="hydrazine",
                    atom_ids=(0, 1, 2, 3),
                    frame=Frame(origin=(-0.71735, 0.0, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("azine_bridge",),
                    metadata={
                        "reactive_atom_id": 0,
                        "anchor_atom_id": 1,
                        "hydrogen_atom_ids": (2, 3),
                        "internal_nitrogen_atom_id": 1,
                    },
                ),
                ReactiveMotif(
                    id="hyd_right",
                    kind="hydrazine",
                    atom_ids=(0, 1, 4, 5),
                    frame=Frame(origin=(0.71735, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("azine_bridge",),
                    metadata={
                        "reactive_atom_id": 1,
                        "anchor_atom_id": 0,
                        "hydrogen_atom_ids": (4, 5),
                        "internal_nitrogen_atom_id": 0,
                    },
                ),
            ),
            atom_symbols=("N", "N", "H", "H", "H", "H"),
            atom_positions=(
                (-0.71735, 0.0, 0.0),
                (0.71735, 0.0, 0.0),
                (-0.9, 0.95, 0.0),
                (-0.9, -0.95, 0.0),
                (0.9, 0.95, 0.0),
                (0.9, -0.95, 0.0),
            ),
            bonds=((0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (1, 4, 1.0), (1, 5, 1.0)),
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
                    allowed_reaction_templates=("azine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 2},
                ),
            ),
            atom_symbols=("C", "O", "C", "H"),
            atom_positions=((0.0, 0.0, 0.0), (0.0, 1.2, 0.0), (1.2, 0.0, 0.0), (-0.8, 0.0, 0.0)),
            bonds=((0, 1, 2.0), (0, 2, 1.0), (0, 3, 1.0)),
        )
        candidate = Candidate(
            id="azine-double-demo",
            score=0.0,
            state=AssemblyState(
                cell=((20.0, 0.0, 0.0), (0.0, 20.0, 0.0), (0.0, 0.0, 10.0)),
                monomer_poses={
                    "m1": Pose(translation=(0.0, 0.0, 0.0)),
                    "m2": Pose(
                        translation=(-1.2, 0.0, 0.0),
                        rotation_matrix=((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
                    ),
                    "m3": Pose(translation=(1.2, 0.0, 0.0)),
                },
                stacking_state="disabled",
            ),
            events=(
                ReactionEvent(
                    id="rxn1",
                    template_id="azine_bridge",
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id="hydrazine", motif_id="hyd_left"),
                        MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
                    ),
                ),
                ReactionEvent(
                    id="rxn2",
                    template_id="azine_bridge",
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id="hydrazine", motif_id="hyd_right"),
                        MotifRef(monomer_instance_id="m3", monomer_id="aldehyde", motif_id="ald1"),
                    ),
                ),
            ),
            metadata={"instance_to_monomer": {"m1": "hydrazine", "m2": "aldehyde", "m3": "aldehyde"}},
        )

        realizer = ReactionRealizer()
        result = realizer.realize(candidate, {"hydrazine": hydrazine, "aldehyde": aldehyde}, {"m1": "hydrazine", "m2": "aldehyde", "m3": "aldehyde"})

        self.assertIsNotNone(result)
        assert result is not None
        self.assertIn("coordinated bridge fit", " ".join(result.metadata["notes"]))

        hydrazine_atoms = {atom.atom_id: atom for atom in result.atoms_by_instance["m1"]}
        left_aldehyde_atoms = {atom.atom_id: atom for atom in result.atoms_by_instance["m2"]}
        right_aldehyde_atoms = {atom.atom_id: atom for atom in result.atoms_by_instance["m3"]}

        left_n_world = self._world_position(candidate.state.monomer_poses["m1"], hydrazine_atoms[0].local_position)
        right_n_world = self._world_position(candidate.state.monomer_poses["m1"], hydrazine_atoms[1].local_position)
        left_c_world = self._world_position(candidate.state.monomer_poses["m2"], left_aldehyde_atoms[0].local_position)
        right_c_world = self._world_position(candidate.state.monomer_poses["m3"], right_aldehyde_atoms[0].local_position)

        self.assertAlmostEqual(realizer._distance(left_n_world, right_n_world), 1.268, delta=0.06)
        self.assertAlmostEqual(realizer._distance(left_c_world, left_n_world), 1.3, delta=0.03)
        self.assertAlmostEqual(realizer._distance(right_c_world, right_n_world), 1.3, delta=0.03)
        left_h_world = self._world_position(candidate.state.monomer_poses["m2"], left_aldehyde_atoms[3].local_position)
        right_h_world = self._world_position(candidate.state.monomer_poses["m3"], right_aldehyde_atoms[3].local_position)
        self.assertAlmostEqual(realizer._distance(left_c_world, left_h_world), 0.8, delta=0.08)
        self.assertAlmostEqual(realizer._distance(right_c_world, right_h_world), 0.8, delta=0.08)
        self.assertLess(realizer._angle(left_c_world, left_n_world, right_n_world), 150.0)
        self.assertLess(realizer._angle(right_c_world, right_n_world, left_n_world), 150.0)

    def test_keto_enamine_realization_removes_oxygen_and_one_amine_hydrogen(self):
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
        candidate = _single_event_candidate(
            "keto_enamine_bridge",
            MotifRef(monomer_instance_id="m1", monomer_id="amine", motif_id="ami1"),
            MotifRef(monomer_instance_id="m2", monomer_id="keto_aldehyde", motif_id="kal1"),
            distance=1.36,
        )

        result = ReactionRealizer().realize(
            candidate,
            {"amine": amine, "keto_aldehyde": keto_aldehyde},
            {"m1": "amine", "m2": "keto_aldehyde"},
        )

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["applied_templates"], {"keto_enamine_bridge": 1})
        self.assertEqual(result.metadata["removed_atom_symbols"], {"H": 2, "O": 1})
        self.assertEqual(len(result.bonds), 2)
        inter_monomer_bond = next(bond for bond in result.bonds if {bond.label_1, bond.label_2} == {"m2_C1", "m1_N1"})
        tautomerized_carbonyl_bond = next(
            bond for bond in result.bonds if {bond.label_1, bond.label_2} == {"m2_C3", "m2_O4"}
        )
        self.assertAlmostEqual(inter_monomer_bond.distance, 1.36, places=6)
        self.assertEqual(inter_monomer_bond.bond_order, 2.0)
        self.assertAlmostEqual(tautomerized_carbonyl_bond.distance, 1.24, places=2)
        self.assertEqual(tautomerized_carbonyl_bond.bond_order, 2.0)
        keto_atoms = {atom.atom_id: atom for atom in result.atoms_by_instance["m2"]}
        self.assertNotIn(1, keto_atoms)
        self.assertNotIn(4, keto_atoms)
        self.assertIn(3, keto_atoms)
        self.assertAlmostEqual(keto_atoms[3].local_position[0], 1.2, places=6)
        self.assertAlmostEqual(keto_atoms[3].local_position[1], 1.24, places=6)
        self.assertEqual(result.metadata["hydrogen_cleanup"]["atom_labels"], ("m2_H6",))
        self._assert_retained_hydrogen_reoriented(
            result=result,
            monomer=keto_aldehyde,
            instance_id="m2",
            parent_atom_id=0,
            hydrogen_atom_id=5,
        )

    def test_boronate_ester_realization_removes_boronic_oxygens_and_four_hydrogens(self):
        boronic_acid = MonomerSpec(
            id="boronic_acid",
            name="minimal boronic acid",
            motifs=(
                ReactiveMotif(
                    id="bor1",
                    kind="boronic_acid",
                    atom_ids=(0, 1, 2, 3, 4, 5),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("boronate_ester_bridge",),
                    metadata={
                        "reactive_atom_id": 1,
                        "anchor_atom_id": 0,
                        "oxygen_atom_ids": (2, 3),
                        "hydrogen_atom_ids": (4, 5),
                    },
                ),
            ),
            atom_symbols=("C", "B", "O", "O", "H", "H"),
            atom_positions=(
                (-1.2, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                (0.9, 0.9, 0.0),
                (0.9, -0.9, 0.0),
                (1.6, 1.2, 0.0),
                (1.6, -1.2, 0.0),
            ),
            bonds=((0, 1, 1.0), (1, 2, 1.0), (1, 3, 1.0), (2, 4, 1.0), (3, 5, 1.0)),
        )
        catechol = MonomerSpec(
            id="catechol",
            name="minimal catechol",
            motifs=(
                ReactiveMotif(
                    id="cat1",
                    kind="catechol",
                    atom_ids=(0, 1, 2, 3, 4, 5),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("boronate_ester_bridge",),
                    metadata={
                        "reactive_atom_id": 0,
                        "anchor_atom_id": 4,
                        "reactive_atom_ids": (0, 1),
                        "anchor_atom_ids": (4, 5),
                        "hydrogen_atom_ids": (2, 3),
                    },
                ),
            ),
            atom_symbols=("O", "O", "H", "H", "C", "C"),
            atom_positions=(
                (0.0, 0.7, 0.0),
                (0.0, -0.7, 0.0),
                (-0.7, 1.2, 0.0),
                (-0.7, -1.2, 0.0),
                (-1.0, 0.7, 0.0),
                (-1.0, -0.7, 0.0),
            ),
            bonds=((0, 2, 1.0), (1, 3, 1.0), (0, 4, 1.0), (1, 5, 1.0)),
        )
        candidate = _single_event_candidate(
            "boronate_ester_bridge",
            MotifRef(monomer_instance_id="m1", monomer_id="boronic_acid", motif_id="bor1"),
            MotifRef(monomer_instance_id="m2", monomer_id="catechol", motif_id="cat1"),
            distance=1.4,
        )

        result = ReactionRealizer().realize(
            candidate,
            {"boronic_acid": boronic_acid, "catechol": catechol},
            {"m1": "boronic_acid", "m2": "catechol"},
        )

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["applied_templates"], {"boronate_ester_bridge": 1})
        self.assertEqual(result.metadata["removed_atom_symbols"], {"O": 2, "H": 4})
        self.assertEqual(len(result.bonds), 2)
        self.assertTrue(all(bond.label_1 == "m1_B2" for bond in result.bonds))
        self.assertEqual({bond.label_2 for bond in result.bonds}, {"m2_O1", "m2_O2"})
        self.assertNotIn("hydrogen_cleanup", result.metadata)

    def test_vinylene_realization_removes_aldehyde_oxygen_and_two_activated_hydrogens(self):
        activated_methylene = MonomerSpec(
            id="activated_methylene",
            name="minimal activated methylene",
            motifs=(
                ReactiveMotif(
                    id="act1",
                    kind="activated_methylene",
                    atom_ids=(0, 1, 2, 3, 4, 5, 6),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("vinylene_bridge",),
                    metadata={
                        "reactive_atom_id": 0,
                        "anchor_atom_id": 1,
                        "hydrogen_atom_ids": (5, 6),
                        "activator_atom_ids": (1, 3),
                    },
                ),
            ),
            atom_symbols=("C", "C", "N", "C", "N", "H", "H"),
            atom_positions=(
                (0.0, 0.0, 0.0),
                (1.2, 0.8, 0.0),
                (2.3, 1.2, 0.0),
                (1.2, -0.8, 0.0),
                (2.3, -1.2, 0.0),
                (-0.7, 0.9, 0.0),
                (-0.7, -0.9, 0.0),
            ),
            bonds=((0, 1, 1.0), (1, 2, 3.0), (0, 3, 1.0), (3, 4, 3.0), (0, 5, 1.0), (0, 6, 1.0)),
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
                    allowed_reaction_templates=("vinylene_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 2},
                ),
            ),
            atom_symbols=("C", "O", "C", "H"),
            atom_positions=((0.0, 0.0, 0.0), (0.0, 1.2, 0.0), (1.2, 0.0, 0.0), (-0.8, 0.0, 0.0)),
            bonds=((0, 1, 2.0), (0, 2, 1.0), (0, 3, 1.0)),
        )
        candidate = _single_event_candidate(
            "vinylene_bridge",
            MotifRef(monomer_instance_id="m1", monomer_id="activated_methylene", motif_id="act1"),
            MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
            distance=1.34,
        )

        result = ReactionRealizer().realize(
            candidate,
            {"activated_methylene": activated_methylene, "aldehyde": aldehyde},
            {"m1": "activated_methylene", "m2": "aldehyde"},
        )

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["applied_templates"], {"vinylene_bridge": 1})
        self.assertEqual(result.metadata["removed_atom_symbols"], {"H": 2, "O": 1})
        self.assertEqual(len(result.bonds), 1)
        self.assertEqual(result.bonds[0].label_1, "m2_C1")
        self.assertEqual(result.bonds[0].label_2, "m1_C1")
        self.assertAlmostEqual(result.bonds[0].distance, 1.34, places=6)
        self.assertEqual(result.metadata["hydrogen_cleanup"]["atom_labels"], ("m2_H4",))
        self._assert_retained_hydrogen_reoriented(
            result=result,
            monomer=aldehyde,
            instance_id="m2",
            parent_atom_id=0,
            hydrogen_atom_id=3,
        )

    def test_custom_realization_registry_can_override_builtin_handler(self):
        amine = MonomerSpec(
            id="amine",
            name="minimal amine",
            motifs=(
                ReactiveMotif(
                    id="ami1",
                    kind="amine",
                    atom_ids=(0, 1, 2, 3),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("imine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 1},
                ),
            ),
            atom_symbols=("N", "C", "H", "H"),
            atom_positions=((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.1, 1.0, 0.0), (0.1, -1.0, 0.0)),
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
                    allowed_reaction_templates=("imine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 2},
                ),
            ),
            atom_symbols=("C", "O", "C", "H"),
            atom_positions=((0.0, 0.0, 0.0), (0.0, 1.2, 0.0), (1.2, 0.0, 0.0), (-0.8, 0.0, 0.0)),
        )
        candidate = Candidate(
            id="imine-demo",
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
                    template_id="imine_bridge",
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id="amine", motif_id="ami1"),
                        MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
                    ),
                ),
            ),
            metadata={"instance_to_monomer": {"m1": "amine", "m2": "aldehyde"}},
        )

        registry = ReactionEventRealizationRegistry()

        def _custom_handler(realizer, event, candidate, monomer_specs):
            del realizer, event, candidate, monomer_specs
            return EventRealization(notes=("custom handler used",))

        registry.register("imine_bridge", _custom_handler)

        result = ReactionRealizer(registry=registry).realize(
            candidate,
            {"amine": amine, "aldehyde": aldehyde},
            {"m1": "amine", "m2": "aldehyde"},
        )

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["applied_event_count"], 1)
        self.assertEqual(result.metadata["notes"], ("custom handler used",))
        self.assertEqual(result.metadata["removed_atom_count"], 0)

    def test_imine_realization_reorients_retained_aldehydic_hydrogen_and_removes_n_hydrogens(self):
        amine = MonomerSpec(
            id="amine",
            name="bonded amine",
            motifs=(
                ReactiveMotif(
                    id="ami1",
                    kind="amine",
                    atom_ids=(0, 1, 2, 3),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("imine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 1},
                ),
            ),
            atom_symbols=("N", "C", "H", "H"),
            atom_positions=((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.1, 1.0, 0.0), (0.1, -1.0, 0.0)),
            bonds=((0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0)),
        )
        aldehyde = MonomerSpec(
            id="aldehyde",
            name="bonded aldehyde",
            motifs=(
                ReactiveMotif(
                    id="ald1",
                    kind="aldehyde",
                    atom_ids=(0, 1, 2, 3),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("imine_bridge",),
                    metadata={"reactive_atom_id": 0, "anchor_atom_id": 2},
                ),
            ),
            atom_symbols=("C", "O", "C", "H"),
            atom_positions=((0.0, 0.0, 0.0), (0.0, 1.2, 0.0), (1.2, 0.0, 0.0), (-0.8, 0.0, 0.0)),
            bonds=((0, 1, 2.0), (0, 2, 1.0), (0, 3, 1.0)),
        )
        candidate = Candidate(
            id="imine-demo",
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
                    template_id="imine_bridge",
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id="amine", motif_id="ami1"),
                        MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
                    ),
                ),
            ),
            metadata={"instance_to_monomer": {"m1": "amine", "m2": "aldehyde"}},
        )

        result = ReactionRealizer().realize(
            candidate,
            {"amine": amine, "aldehyde": aldehyde},
            {"m1": "amine", "m2": "aldehyde"},
        )

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["hydrogen_cleanup"]["n_overrides"], 1)

        amine_atoms = {atom.atom_id: atom for atom in result.atoms_by_instance["m1"]}
        aldehyde_atoms = {atom.atom_id: atom for atom in result.atoms_by_instance["m2"]}

        self.assertNotIn(2, amine_atoms)
        self.assertNotIn(3, amine_atoms)

        self.assertIn(3, aldehyde_atoms)
        self.assertAlmostEqual(aldehyde_atoms[3].local_position[0], 0.0, places=6)
        self.assertGreater(aldehyde_atoms[3].local_position[1], 0.7)


if __name__ == "__main__":
    unittest.main()
