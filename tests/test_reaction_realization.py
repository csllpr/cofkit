import contextlib
import io
import sys
import tempfile
import unittest
from math import acos, pi
from pathlib import Path


from cofkit import AssemblyState, Candidate, CIFWriter, Frame, MonomerSpec, Pose, ReactiveMotif, ReactionEvent, MotifRef, build_rdkit_monomer
from cofkit.bond_types import cif_type_to_bond_order
from cofkit.decompose import _build_bonded_mol, _eligible_beta_ketoenamine_single_bonds, read_periodic_cif_atoms
from cofkit.linkage_geometry import BORONATE_ESTER_BOND_TARGET_DISTANCE, BORONATE_ESTER_OBO_TARGET_ANGLE_DEG, required_bridge_span
from cofkit import reaction_realization
from cofkit.reaction_realization import EventRealization, ReactionEventRealizationRegistry, ReactionRealizer


def _boronate_ester_monomers() -> tuple[MonomerSpec, MonomerSpec]:
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
    return boronic_acid, catechol


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

    def _distance(self, first, second) -> float:
        return sum((a - b) ** 2 for a, b in zip(first, second)) ** 0.5

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

        self.assertAlmostEqual(realizer._distance(left_n_world, right_n_world), 1.408, delta=0.06)
        self.assertAlmostEqual(realizer._distance(left_c_world, left_n_world), 1.3, delta=0.03)
        self.assertAlmostEqual(realizer._distance(right_c_world, right_n_world), 1.3, delta=0.03)
        left_h_world = self._world_position(candidate.state.monomer_poses["m2"], left_aldehyde_atoms[3].local_position)
        right_h_world = self._world_position(candidate.state.monomer_poses["m3"], right_aldehyde_atoms[3].local_position)
        self.assertAlmostEqual(realizer._distance(left_c_world, left_h_world), 0.8, delta=0.08)
        self.assertAlmostEqual(realizer._distance(right_c_world, right_h_world), 0.8, delta=0.08)
        self.assertLess(realizer._angle(left_c_world, left_n_world, right_n_world), 150.0)
        self.assertLess(realizer._angle(right_c_world, right_n_world, left_n_world), 150.0)

    def test_imine_chain_closure_fit_failure_replaces_success_note(self):
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
        # Place the aldehyde so its anchor atom coincides with the amine anchor:
        # the degenerate bridge axis makes the chain-closure fit return None.
        candidate = _single_event_candidate(
            "imine_bridge",
            MotifRef(monomer_instance_id="m1", monomer_id="amine", motif_id="ami1"),
            MotifRef(monomer_instance_id="m2", monomer_id="aldehyde", motif_id="ald1"),
            distance=-2.2,
        )

        result = ReactionRealizer().realize(candidate, {"amine": amine, "aldehyde": aldehyde}, {"m1": "amine", "m2": "aldehyde"})

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["applied_event_count"], 1)
        notes = " ".join(result.metadata["notes"])
        self.assertNotIn("applies a local imine chain-closure fit", notes)
        self.assertIn("imine chain-closure fit could not run", notes)
        self.assertIn("near-collinear", notes)
        # No fitted overrides: the exported bond keeps the unfitted geometry.
        self.assertAlmostEqual(result.bonds[0].distance, 2.2, places=6)

    def test_azine_fit_skip_reason_recorded_for_single_event_hydrazine(self):
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
        notes = " ".join(result.metadata["notes"])
        self.assertIn("Azine bridge fit skipped for hydrazine instance m1", notes)
        self.assertIn("1 reaction events (expected 2)", notes)
        self.assertNotIn("coordinated bridge fit", notes)

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
        self.assertEqual(len(result.bonds), 3)
        inter_monomer_bond = next(bond for bond in result.bonds if {bond.label_1, bond.label_2} == {"m2_C1", "m1_N1"})
        enamine_double_bond = next(
            bond for bond in result.bonds if {bond.label_1, bond.label_2} == {"m2_C1", "m2_C3"}
        )
        tautomerized_carbonyl_bond = next(
            bond for bond in result.bonds if {bond.label_1, bond.label_2} == {"m2_C3", "m2_O4"}
        )
        self.assertAlmostEqual(inter_monomer_bond.distance, 1.36, places=6)
        self.assertEqual(inter_monomer_bond.bond_order, 1.0)
        self.assertEqual(enamine_double_bond.bond_order, 2.0)
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

    def test_beta_ketoenol_route_retains_preexisting_carbonyl(self):
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
            id="beta_ketoenol",
            name="minimal beta-ketoenol route precursor",
            motifs=(
                ReactiveMotif(
                    id="kal1",
                    kind="keto_aldehyde",
                    atom_ids=(0, 1, 2, 3, 4, 6, 7, 8),
                    frame=Frame(origin=(0.0, 0.0, 0.0), primary=(-1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)),
                    allowed_reaction_templates=("keto_enamine_bridge",),
                    metadata={
                        "precursor_route": "beta_ketoenol_michael_addition",
                        "reactive_atom_id": 0,
                        "anchor_atom_id": 2,
                        "aldehyde_oxygen_atom_id": 1,
                        "beta_ketoenol_alpha_carbon_atom_id": 2,
                        "beta_keto_carbonyl_carbon_atom_id": 3,
                        "beta_keto_carbonyl_oxygen_atom_id": 4,
                    },
                ),
            ),
            atom_symbols=("C", "O", "C", "C", "O", "C", "H", "H", "H"),
            atom_positions=(
                (0.0, 0.0, 0.0),
                (0.0, 1.2, 0.0),
                (1.3, 0.0, 0.0),
                (2.6, 0.0, 0.0),
                (2.6, 1.2, 0.0),
                (3.9, 0.0, 0.0),
                (-0.8, 0.0, 0.0),
                (1.3, 0.9, 0.0),
                (1.3, -0.9, 0.0),
            ),
            bonds=(
                (0, 1, 2.0),
                (0, 2, 1.0),
                (0, 6, 1.0),
                (2, 3, 1.0),
                (2, 7, 1.0),
                (2, 8, 1.0),
                (3, 4, 2.0),
                (3, 5, 1.0),
            ),
        )
        candidate = _single_event_candidate(
            "keto_enamine_bridge",
            MotifRef(monomer_instance_id="m1", monomer_id="amine", motif_id="ami1"),
            MotifRef(monomer_instance_id="m2", monomer_id="beta_ketoenol", motif_id="kal1"),
            distance=1.36,
        )

        result = ReactionRealizer().realize(
            candidate,
            {"amine": amine, "beta_ketoenol": keto_aldehyde},
            {"m1": "amine", "m2": "beta_ketoenol"},
        )

        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result.metadata["removed_atom_symbols"], {"H": 2, "O": 1})
        self.assertEqual(len(result.bonds), 2)
        self.assertTrue(
            any(
                {bond.label_1, bond.label_2} == {"m2_C1", "m2_C3"}
                and bond.bond_order == 2.0
                for bond in result.bonds
            )
        )
        retained_atoms = {atom.atom_id for atom in result.atoms_by_instance["m2"]}
        self.assertNotIn(1, retained_atoms)
        self.assertIn(4, retained_atoms)

    def test_boronate_ester_realization_removes_boronic_oxygens_and_four_hydrogens(self):
        boronic_acid, catechol = _boronate_ester_monomers()
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

    def test_boronate_ester_realization_closes_both_bonds_at_target_distance(self):
        boronic_acid, catechol = _boronate_ester_monomers()
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
        for bond in result.bonds:
            self.assertAlmostEqual(bond.distance, BORONATE_ESTER_BOND_TARGET_DISTANCE, places=6)
        # The closure fit keeps every anchor-to-reactive bond length fixed: the boron
        # moves on its C-B sphere and each catechol oxygen on its C-O sphere.
        boron_atom = next(atom for atom in result.atoms_by_instance["m1"] if atom.atom_id == 1)
        self.assertAlmostEqual(self._distance(boron_atom.local_position, (-1.2, 0.0, 0.0)), 1.2, places=6)
        oxygen_positions = {
            atom.atom_id: atom.local_position
            for atom in result.atoms_by_instance["m2"]
            if atom.atom_id in (0, 1)
        }
        self.assertAlmostEqual(self._distance(oxygen_positions[0], (-1.0, 0.7, 0.0)), 1.0, places=6)
        self.assertAlmostEqual(self._distance(oxygen_positions[1], (-1.0, -0.7, 0.0)), 1.0, places=6)
        # The ring contracts toward the measured baseline O-B-O angle instead of
        # keeping the free-catechol opening (~139 degrees in this fixture, whose
        # artificial bond lengths make the exact target infeasible).
        boron_world = boron_atom.local_position
        vector_1 = tuple(a - b for a, b in zip(oxygen_positions[0], boron_world))
        vector_2 = tuple(a - b for a, b in zip(oxygen_positions[1], boron_world))
        norm_1 = self._distance(vector_1, (0.0, 0.0, 0.0))
        norm_2 = self._distance(vector_2, (0.0, 0.0, 0.0))
        cosine = sum(a * b for a, b in zip(vector_1, vector_2)) / (norm_1 * norm_2)
        obo_angle_deg = acos(max(-1.0, min(1.0, cosine))) * 180.0 / pi
        self.assertLess(obo_angle_deg, 135.0)
        self.assertGreater(obo_angle_deg, BORONATE_ESTER_OBO_TARGET_ANGLE_DEG - 15.0)

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


def _real_keto_enamine_candidate(keto_smiles: str, *, site_count: int = 1):
    """Build a keto_enamine_bridge candidate from real RDKit-derived monomers.

    Each site pairs a fresh aniline instance with one keto-aldehyde motif of
    the keto monomer, placed so the reacting N...C atoms sit 1.36 angstrom
    apart along +x.
    """
    keto = build_rdkit_monomer("keto", "keto-aldehyde precursor", keto_smiles, "keto_aldehyde")
    amine = build_rdkit_monomer("aniline", "aniline", "Nc1ccccc1", "amine")
    amine_motif = amine.motifs[0]
    nitrogen_atom_id = amine_motif.metadata["reactive_atom_id"]
    nitrogen_position = amine.atom_positions[nitrogen_atom_id]

    poses = {"m_keto": Pose(translation=(15.0, 15.0, 7.5))}
    instance_to_monomer = {"m_keto": "keto"}
    events = []
    for index in range(site_count):
        keto_motif = keto.motifs[index]
        carbon_atom_id = keto_motif.metadata["reactive_atom_id"]
        carbon_position = keto.atom_positions[carbon_atom_id]
        amine_instance_id = f"m_amine{index}"
        instance_to_monomer[amine_instance_id] = "aniline"
        poses[amine_instance_id] = Pose(
            translation=(
                15.0 + carbon_position[0] - nitrogen_position[0] + 1.36,
                15.0 + carbon_position[1] - nitrogen_position[1],
                7.5 + carbon_position[2] - nitrogen_position[2],
            )
        )
        events.append(
            ReactionEvent(
                id=f"rxn{index}",
                template_id="keto_enamine_bridge",
                participants=(
                    MotifRef(monomer_instance_id=amine_instance_id, monomer_id="aniline", motif_id=amine_motif.id),
                    MotifRef(monomer_instance_id="m_keto", monomer_id="keto", motif_id=keto_motif.id),
                ),
            )
        )
    candidate = Candidate(
        id="keto-enamine-real-monomer-demo",
        score=0.0,
        state=AssemblyState(
            cell=((40.0, 0.0, 0.0), (0.0, 40.0, 0.0), (0.0, 0.0, 20.0)),
            monomer_poses=poses,
            stacking_state="disabled",
        ),
        events=tuple(events),
        metadata={"instance_to_monomer": instance_to_monomer},
    )
    specs = {"keto": keto, "aniline": amine}
    return candidate, specs, instance_to_monomer


def _exported_bond_order_map(cif_text: str) -> dict[frozenset[str], float]:
    """Parse the exported CIF bond loop into a label-pair -> bond-order map."""

    orders: dict[frozenset[str], float] = {}
    for line in cif_text.splitlines():
        tokens = line.split()
        if len(tokens) != 6 or "_" not in tokens[1]:
            continue
        order = cif_type_to_bond_order(tokens[5])
        assert order is not None, f"unparseable bond type in line: {line}"
        orders[frozenset((tokens[0], tokens[1]))] = order
    return orders


class KetoEnamineQuinoidRedistributionTests(unittest.TestCase):
    """Quinoid/cyclohexadienone ring re-bond-ordering for real aromatic monomers."""

    def test_quinoid_redistribution_fixes_carbonyl_valence_and_ring_pattern(self):
        candidate, specs, instance_to_monomer = _real_keto_enamine_candidate("O=Cc1ccccc1O")

        result = ReactionRealizer().realize(candidate, specs, instance_to_monomer)

        self.assertIsNotNone(result)
        assert result is not None
        keto = specs["keto"]
        motif = keto.motifs[0]
        linkage_anchor_atom_id = motif.metadata["anchor_atom_id"]
        carbonyl_anchor_atom_id = motif.metadata["ortho_hydroxyl_anchor_atom_id"]
        report = result.metadata["quinoid_redistribution"]
        self.assertEqual(report["site_count"], 1)
        self.assertEqual(report["redistributed_site_count"], 1)
        self.assertEqual(report["unresolved_sites"], ())
        self.assertEqual(report["changed_ring_bond_count"], 6)

        # The realization bond list now carries the rewritten ring bonds.
        label = lambda atom_id: ReactionRealizer.atom_label("m_keto", "C", atom_id)
        realized_orders = {
            frozenset((bond.label_1, bond.label_2)): bond.bond_order for bond in result.bonds
        }
        ring_edges = [
            (first, second)
            for first, second, order in keto.bonds
            if abs(order - 1.5) <= 1.0e-6
        ]
        self.assertEqual(len(ring_edges), 6)
        for first, second in ring_edges:
            order = realized_orders[frozenset((label(first), label(second)))]
            if carbonyl_anchor_atom_id in (first, second) or linkage_anchor_atom_id in (first, second):
                self.assertEqual(order, 1.0, f"anchor ring bond {first}-{second} must become single")
        residual_atom_ids = {
            atom_id for edge in ring_edges for atom_id in edge
        } - {carbonyl_anchor_atom_id, linkage_anchor_atom_id}
        self.assertEqual(len(residual_atom_ids), 4)
        for atom_id in residual_atom_ids:
            double_count = sum(
                1
                for first, second in ring_edges
                if atom_id in (first, second)
                and realized_orders[frozenset((label(first), label(second)))] == 2.0
            )
            self.assertEqual(double_count, 1, f"ring atom {atom_id} must take exactly one quinoid double bond")

        # The exported CIF (the real output contract) carries no aromatic "A"
        # type on the quinoid ring and a valence-4 carbonyl carbon.
        exported = CIFWriter().export_candidate(candidate, specs)
        exported_orders = _exported_bond_order_map(exported.text)
        keto_internal = {
            pair: order
            for pair, order in exported_orders.items()
            if all(label.startswith("m_keto_") for label in pair)
        }
        for first, second in ring_edges:
            self.assertIn(frozenset((label(first), label(second))), keto_internal)
        valences: dict[str, float] = {}
        for pair, order in exported_orders.items():
            for atom_label in pair:
                valences[atom_label] = valences.get(atom_label, 0.0) + order
        carbonyl_label = label(carbonyl_anchor_atom_id)
        self.assertAlmostEqual(valences[carbonyl_label], 4.0, places=6)
        for atom_id in residual_atom_ids | {linkage_anchor_atom_id}:
            self.assertAlmostEqual(valences[label(atom_id)], 4.0, places=6)
        for first, second in ring_edges:
            pair = frozenset((label(first), label(second)))
            self.assertIn(keto_internal[pair], (1.0, 2.0))

    def test_quinoid_redistribution_saturates_triformylphloroglucinol_ring(self):
        candidate, specs, instance_to_monomer = _real_keto_enamine_candidate(
            "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O",
            site_count=3,
        )

        result = ReactionRealizer().realize(candidate, specs, instance_to_monomer)

        self.assertIsNotNone(result)
        assert result is not None
        report = result.metadata["quinoid_redistribution"]
        self.assertEqual(report["site_count"], 3)
        self.assertEqual(report["redistributed_site_count"], 3)
        self.assertEqual(report["unresolved_sites"], ())
        keto = specs["keto"]
        ring_edges = [
            frozenset((first, second))
            for first, second, order in keto.bonds
            if abs(order - 1.5) <= 1.0e-6
        ]
        self.assertEqual(len(ring_edges), 6)
        realized_orders = {
            frozenset((bond.label_1, bond.label_2)): bond.bond_order for bond in result.bonds
        }
        label = lambda atom_id: ReactionRealizer.atom_label("m_keto", "C", atom_id)
        # Every ring atom carries an exocyclic double bond (C=O or C=C), so
        # the cyclohexanetrione-style ring is all single bonds.
        for edge in ring_edges:
            pair = frozenset((label(atom_id) for atom_id in sorted(edge)))
            self.assertEqual(realized_orders[pair], 1.0)

    def test_quinoid_redistribution_warns_and_reports_when_no_pattern_exists(self):
        # 5-hydroxyfurfural: the five-membered aromatic ring admits no
        # alternating double-bond pattern once both anchors are saturated.
        candidate, specs, instance_to_monomer = _real_keto_enamine_candidate("O=Cc1occc1O")

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            result = ReactionRealizer().realize(candidate, specs, instance_to_monomer)

        self.assertIsNotNone(result)
        assert result is not None
        self.assertIn("warning: keto-enamine quinoid redistribution unresolved", stderr.getvalue())
        report = result.metadata["quinoid_redistribution"]
        self.assertEqual(report["site_count"], 1)
        self.assertEqual(report["redistributed_site_count"], 0)
        self.assertEqual(report["changed_ring_bond_count"], 0)
        self.assertEqual(len(report["unresolved_sites"]), 1)
        self.assertIn("no alternating double-bond pattern", report["unresolved_sites"][0]["reason"])
        self.assertTrue(
            any("quinoid redistribution could not determine" in note for note in result.metadata["notes"])
        )
        # Documented degradation: the ring keeps precursor aromatic orders.
        exported = CIFWriter().export_candidate(candidate, specs)
        keto_ring_orders = [
            tokens[5]
            for line in exported.text.splitlines()
            if len((tokens := line.split())) == 6
            and "_" in tokens[1]
            and tokens[0].startswith("m_keto_")
            and tokens[1].startswith("m_keto_")
        ]
        self.assertIn("A", keto_ring_orders)

    def test_quinoid_product_is_recognized_by_decompose_side(self):
        candidate, specs, instance_to_monomer = _real_keto_enamine_candidate("O=Cc1ccccc1O")
        exported = CIFWriter().export_candidate(candidate, specs)

        with tempfile.TemporaryDirectory() as temp_dir:
            cif_path = Path(temp_dir) / "product.cif"
            cif_path.write_text(exported.text)
            built = _build_bonded_mol(read_periodic_cif_atoms(cif_path))

        # No valence overflow survives the export: every carbon is at most 4.
        for atom in built.mol.GetAtoms():
            valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
            if atom.GetAtomicNum() == 6:
                self.assertLessEqual(valence, 4.0 + 1.0e-6)
        # The decompose side's keto-enamine environment check recognizes the
        # written C-N single bond with its C=C anchor and ring carbonyl.
        eligible = _eligible_beta_ketoenamine_single_bonds(built.mol)
        self.assertEqual(len(eligible), 1)


def _minimal_imine_monomers() -> tuple[MonomerSpec, MonomerSpec]:
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
    return amine, aldehyde


class ImineBridgeConstructorTests(unittest.TestCase):
    """Single-owner closed-form imine bridge construction (audit #15 rework).

    Fixture geometry: amine N at the origin with its anchor at (-1, 0, 0);
    aldehyde translated by t along +x so its anchor sits at (t + 1.2, 0, 0).
    The anchor span is therefore d = t + 2.2 with side lengths l_C = 1.2,
    l_N = 1.0 and the C=N target 1.3.
    """

    def _candidate(self, translation: float) -> Candidate:
        return Candidate(
            id="imine-constructor-demo",
            score=0.0,
            state=AssemblyState(
                cell=((60.0, 0.0, 0.0), (0.0, 60.0, 0.0), (0.0, 0.0, 10.0)),
                monomer_poses={
                    "m1": Pose(translation=(0.0, 0.0, 0.0)),
                    "m2": Pose(translation=(translation, 0.0, 0.0)),
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

    def _realized_geometry(self, translation: float):
        amine, aldehyde = _minimal_imine_monomers()
        candidate = self._candidate(translation)
        realizer = ReactionRealizer()
        result = realizer.realize(candidate, {"amine": amine, "aldehyde": aldehyde}, {"m1": "amine", "m2": "aldehyde"})
        self.assertIsNotNone(result)
        assert result is not None
        realized = {
            instance_id: {atom.atom_id: atom.local_position for atom in atoms}
            for instance_id, atoms in result.atoms_by_instance.items()
        }
        pose_m2 = candidate.state.monomer_poses["m2"]
        carbon_world = realizer._world_position(pose_m2, realized["m2"][0])
        carbon_anchor_world = realizer._world_position(pose_m2, realized["m2"][2])
        nitrogen_world = realizer._world_position(candidate.state.monomer_poses["m1"], realized["m1"][0])
        nitrogen_anchor_world = realizer._world_position(candidate.state.monomer_poses["m1"], realized["m1"][1])
        return realizer, result, carbon_world, carbon_anchor_world, nitrogen_world, nitrogen_anchor_world

    def test_consistent_span_closes_exactly_at_priors(self):
        span = required_bridge_span(1.2, 1.3, 1.0, 120.0, 120.0)
        self.assertIsNotNone(span)
        assert span is not None
        realizer, result, carbon, anchor_c, nitrogen, anchor_n = self._realized_geometry(span - 2.2)

        carbon_angle = realizer._angle(anchor_c, carbon, nitrogen)
        nitrogen_angle = realizer._angle(carbon, nitrogen, anchor_n)
        self.assertAlmostEqual(carbon_angle, 120.0, delta=1e-9)
        self.assertAlmostEqual(nitrogen_angle, 120.0, delta=1e-9)
        # All three lengths exact: aryl-C-C, C=N (exported bond), N-aryl-C.
        self.assertAlmostEqual(realizer._distance(anchor_c, carbon), 1.2, delta=1e-9)
        self.assertAlmostEqual(result.bonds[0].distance, 1.3, delta=1e-9)
        self.assertAlmostEqual(realizer._distance(nitrogen, anchor_n), 1.0, delta=1e-9)
        # E-configuration preserved: the carbon stays on the (removed)
        # oxygen's side of the anchor axis, the nitrogen on the opposite side.
        self.assertGreater(carbon[1], 0.0)
        self.assertLess(nitrogen[1], 0.0)
        notes = " ".join(result.metadata["notes"])
        self.assertIn("closed exactly at the template angle priors", notes)
        self.assertNotIn("best-effort", notes)

    def test_inconsistent_span_is_best_effort_with_honest_diagnostic(self):
        span = required_bridge_span(1.2, 1.3, 1.0, 120.0, 120.0)
        assert span is not None
        realizer, result, carbon, anchor_c, nitrogen, anchor_n = self._realized_geometry(span - 2.2 + 0.1)

        # Lengths stay exact; only the angles absorb the span mismatch, and
        # they absorb it together (least-squares on both residuals).
        self.assertAlmostEqual(realizer._distance(anchor_c, carbon), 1.2, delta=1e-9)
        self.assertAlmostEqual(result.bonds[0].distance, 1.3, delta=1e-9)
        self.assertAlmostEqual(realizer._distance(nitrogen, anchor_n), 1.0, delta=1e-9)
        carbon_angle = realizer._angle(anchor_c, carbon, nitrogen)
        nitrogen_angle = realizer._angle(carbon, nitrogen, anchor_n)
        # +0.1 angstrom past the consistent span flattens the zig-zag:
        # independently measured optimum is 127.97 / 126.98 degrees.
        self.assertGreater(carbon_angle, 124.0)
        self.assertLess(carbon_angle, 132.0)
        self.assertGreater(nitrogen_angle, 124.0)
        self.assertLess(nitrogen_angle, 132.0)
        notes = " ".join(result.metadata["notes"])
        self.assertIn("best-effort imine bridge construction", notes)
        self.assertIn("deviates from its prior-consistent span", notes)
        self.assertNotIn("closed exactly", notes)

    def test_infeasible_span_falls_back_to_existing_failure_note(self):
        span = required_bridge_span(1.2, 1.3, 1.0, 120.0, 120.0)
        assert span is not None
        realizer, result, carbon, anchor_c, nitrogen, anchor_n = self._realized_geometry(20.0)

        notes = " ".join(result.metadata["notes"])
        self.assertIn("imine chain-closure fit could not run", notes)
        self.assertIn("near-collinear", notes)
        self.assertAlmostEqual(result.bonds[0].distance, 20.0, places=6)

    def test_priors_are_read_from_the_template_profile(self):
        from cofkit.reactions import BridgeGeometryPriors

        span_125 = required_bridge_span(1.2, 1.3, 1.0, 125.0, 125.0)
        self.assertIsNotNone(span_125)
        assert span_125 is not None
        translation = span_125 - 2.2

        # With the default 120-degree priors the same placement is only
        # best-effort...
        _, result_default, _, _, _, _ = self._realized_geometry(translation)
        self.assertIn("best-effort", " ".join(result_default.metadata["notes"]))

        # ...but with monkeypatched profile priors the identical geometry
        # closes exactly, proving the constructor reads the profile.
        original = reaction_realization.bridge_geometry_priors
        reaction_realization.bridge_geometry_priors = lambda template_id: BridgeGeometryPriors(
            carbon_angle_deg=125.0,
            nitrogen_angle_deg=125.0,
        )
        try:
            realizer, result, carbon, anchor_c, nitrogen, anchor_n = self._realized_geometry(translation)
        finally:
            reaction_realization.bridge_geometry_priors = original
        notes = " ".join(result.metadata["notes"])
        self.assertIn("closed exactly at the template angle priors (125.0/125.0", notes)
        self.assertAlmostEqual(realizer._angle(anchor_c, carbon, nitrogen), 125.0, delta=1e-9)
        self.assertAlmostEqual(realizer._angle(carbon, nitrogen, anchor_n), 125.0, delta=1e-9)


class AzineBridgeConstructorTests(unittest.TestCase):
    def _azine_group(self, aldehyde_translation: float):
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
            id="azine-constructor-demo",
            score=0.0,
            state=AssemblyState(
                cell=((60.0, 0.0, 0.0), (0.0, 60.0, 0.0), (0.0, 0.0, 10.0)),
                monomer_poses={
                    "m1": Pose(translation=(0.0, 0.0, 0.0)),
                    "m2": Pose(
                        translation=(-aldehyde_translation, 0.0, 0.0),
                        rotation_matrix=((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0)),
                    ),
                    "m3": Pose(translation=(aldehyde_translation, 0.0, 0.0)),
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
        specs = {"hydrazine": hydrazine, "aldehyde": aldehyde}
        return candidate, specs

    def test_infeasible_endpoint_triangle_is_best_effort_with_diagnostic(self):
        # Anchor-N distance 3.496 angstrom > l_C + C=N target (2.5): the
        # endpoint triangle cannot close, so the carbon goes on the
        # anchor<->N line at l_C and the exported C=N lands at the closest
        # feasible distance 3.496 - 1.2 = 2.296 with an honest note.
        candidate, specs = self._azine_group(3.0)
        realizer = ReactionRealizer()
        result = realizer.realize(candidate, specs, {"m1": "hydrazine", "m2": "aldehyde", "m3": "aldehyde"})

        self.assertIsNotNone(result)
        assert result is not None
        notes = " ".join(result.metadata["notes"])
        self.assertIn("best-effort", notes)
        self.assertIn("cannot close the endpoint triangle", notes)
        self.assertIn("2.296", notes)
        cn_bond_distances = sorted(
            bond.distance for bond in result.bonds if bond.bond_order == 2.0
        )
        self.assertEqual(len(cn_bond_distances), 2)
        for distance in cn_bond_distances:
            self.assertAlmostEqual(distance, 2.296, delta=0.01)
        # The coordinated N-N closure is unaffected by the endpoint fallback.
        hydrazine_atoms = {atom.atom_id: atom for atom in result.atoms_by_instance["m1"]}
        pose_m1 = candidate.state.monomer_poses["m1"]
        left_n = realizer._world_position(pose_m1, hydrazine_atoms[0].local_position)
        right_n = realizer._world_position(pose_m1, hydrazine_atoms[1].local_position)
        self.assertAlmostEqual(realizer._distance(left_n, right_n), 1.408, delta=0.01)


if __name__ == "__main__":
    unittest.main()
