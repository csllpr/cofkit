import json
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator
from cofkit.build_workflows.ring_forming import RingFormationConfig, RingFormingStructureGenerator
from cofkit.chem.rdkit import build_rdkit_monomer
from cofkit.cif import CIFWriter
from cofkit.cofid import generate_cofid
from cofkit.decompose import (
    _BOROXINE_SPEC,
    BondedMolBuildResult,
    _build_bonded_mol,
    decompose_cif_to_cofid,
)
from cofkit.decompose_cif import PeriodicCifAtoms, read_periodic_cif_atoms
from cofkit.decompose_events import (
    CHARGE_SCOPE_GUEST_LOCALIZED,
    CHARGE_SCOPE_LIMITED,
    CHARGE_SCOPE_NOT_PRESENT,
    CHARGE_SCOPE_PRESERVED,
    CHARGE_SCOPE_RESTORED_BY_HEURISTIC,
    EVENT_STATUS_UNEXPLAINED,
    STEREO_SCOPE_NOT_DETECTED,
    STEREO_SCOPE_UNSUPPORTED,
    EventDetectionResult,
    LinkageEvent,
    ReconstructionHypothesis,
    _cut_and_reconstruct,
    _select_event_result,
)
from cofkit.decompose_ledger import (
    build_atom_ledger,
    build_recovered_fragment_record,
    perceive_stereochemistry_evidence,
)


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
BDBA = "OB(O)c1ccc(B(O)O)cc1"


def _record(record_id: str, smiles: str, motif_kind: str, connectivity: int) -> BatchMonomerRecord:
    return BatchMonomerRecord(
        id=record_id,
        name=record_id,
        smiles=smiles,
        motif_kind=motif_kind,
        expected_connectivity=connectivity,
    )


def _write_imine_cif(root: Path) -> Path:
    generator = BatchStructureGenerator(
        BatchGenerationConfig(
            allowed_reactions=("imine_bridge",),
            single_node_topology_ids=("hcb",),
            topology_ids=("hcb",),
            enumerate_all_topologies=False,
            write_cif=True,
            rdkit_num_conformers=2,
        )
    )
    summary, _candidate = generator.generate_pair_candidate(
        _record("tapb", TAPB, "amine", 3),
        _record("tfb", TFB, "aldehyde", 3),
        out_dir=root,
        write_cif=True,
    )
    return Path(summary.cif_path)


def _write_boroxine_cif(root: Path) -> Path:
    monomer = build_rdkit_monomer(
        "cof1_precursor",
        "benzene-1,4-diboronic acid",
        BDBA,
        "boronic_acid",
        num_conformers=1,
    )
    generator = RingFormingStructureGenerator(RingFormationConfig(stacking_ids=()))
    candidate = generator.generate(monomer, "boroxine_trimerization")
    cofid = generate_cofid(candidate, {monomer.id: monomer})
    path = root / "boroxine.cif"
    CIFWriter().write_candidate(path, candidate, {monomer.id: monomer}, cofid=cofid)
    return path


def _cif_lines_without_cofid(source_path: str | Path) -> list[str]:
    lines = Path(source_path).read_text(encoding="utf-8").splitlines()
    if lines and lines[0].startswith("# COFid: "):
        lines = lines[1:]
    return lines


def _append_guest_atoms(
    source_path: str | Path,
    target_path: Path,
    atom_rows: tuple[str, ...],
    *,
    bond_rows: tuple[str, ...] = (),
    charge_by_label: dict[str, int] | None = None,
) -> Path:
    """Append guest atom rows (and optional bonds/charges) to a written CIF."""
    lines = _cif_lines_without_cofid(source_path)
    label_header_index = next(
        index for index, line in enumerate(lines) if line.strip() == "_atom_site_label"
    )
    loop_start = label_header_index
    while lines[loop_start].strip() != "loop_":
        loop_start -= 1
    header_end = label_header_index
    while lines[header_end + 1].lstrip().startswith("_"):
        header_end += 1
    row_end = header_end + 1
    while (
        row_end < len(lines)
        and lines[row_end].strip()
        and not lines[row_end].lstrip().startswith(("_", "loop_", "data_", "#"))
    ):
        row_end += 1
    headers = [line.strip() for line in lines[loop_start + 1 : header_end + 1]]
    rows = list(lines[header_end + 1 : row_end])
    guest_rows = list(atom_rows)
    if charge_by_label is not None:
        headers.append("_atom_site_pdbx_formal_charge")
        label_column = headers.index("_atom_site_label")
        rows = [
            f"{row} {charge_by_label.get(row.split()[label_column], 0)}"
            for row in rows
        ]
        guest_rows = [
            f"{row} {charge_by_label.get(row.split()[0], 0)}" for row in guest_rows
        ]
    updated = (
        lines[: loop_start + 1]
        + headers
        + rows
        + guest_rows
        + lines[row_end:]
        + list(bond_rows)
    )
    target_path.write_text("\n".join(updated) + "\n", encoding="utf-8")
    return target_path


class EventAtomLedgerBalanceTests(unittest.TestCase):
    def test_imine_ledger_balances_with_reaction_stoichiometry(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cif_path = _write_imine_cif(temp_path)
            stripped = temp_path / "imine_stripped.cif"
            stripped.write_text(
                "\n".join(_cif_lines_without_cofid(cif_path)) + "\n", encoding="utf-8"
            )
            input_atom_count = len(read_periodic_cif_atoms(stripped))

            result = decompose_cif_to_cofid(stripped, topology="hcb", linkage="imine")

            self.assertTrue(result.ok, result.reason)
            ledger = result.metadata.get("atom_ledger")
            self.assertIsNotNone(ledger)
            n_imine_bonds = int(result.metadata["n_events"])
            self.assertEqual(n_imine_bonds, 3)
            # One aldehyde oxygen restored per imine condensation (H2O per bond).
            self.assertEqual(ledger["reaction_additions"], {"O": n_imine_bonds})
            self.assertEqual(ledger["guest_molecule_count"], 0)
            self.assertEqual(ledger["guest_atom_count"], 0)
            self.assertEqual(ledger["residue_atom_count"], 0)
            self.assertEqual(ledger["unaccounted_atom_count"], 0)
            self.assertEqual(ledger["input_atom_count"], input_atom_count)
            balance = ledger["balance"]
            self.assertTrue(balance["balanced"], balance)
            self.assertTrue(
                all(residual == 0 for residual in balance["per_element_residuals"].values())
            )
            self.assertEqual(
                balance["input_atom_count"] + sum(balance["reaction_additions"].values()),
                balance["recovered_precursor_explicit_atom_count"]
                + sum(balance["reaction_deletions"].values())
                + balance["residue_atom_count"]
                + balance["guest_atom_count"]
                + balance["unaccounted_atom_count"],
            )
            self.assertGreater(ledger["recovered_precursor_implicit_hydrogen_count"], 0)
            self.assertEqual(ledger["charge_scope"]["status"], CHARGE_SCOPE_NOT_PRESENT)
            self.assertEqual(ledger["stereo_scope"]["status"], STEREO_SCOPE_NOT_DETECTED)
            # The full result must remain JSON-serializable with the ledger attached.
            payload = result.to_dict()
            json.dumps(payload)
            self.assertTrue(payload["metadata"]["atom_ledger"]["balance"]["balanced"])
            hypothesis_summary = payload["metadata"]["hypotheses"][0]["atom_ledger_summary"]
            self.assertTrue(hypothesis_summary["balanced"])
            self.assertEqual(hypothesis_summary["guest_molecule_count"], 0)

    def test_boroxine_ledger_balances_with_ring_residue_stoichiometry(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cif_path = _write_boroxine_cif(temp_path)
            stripped = temp_path / "boroxine_stripped.cif"
            stripped.write_text(
                "\n".join(_cif_lines_without_cofid(cif_path)) + "\n", encoding="utf-8"
            )

            result = decompose_cif_to_cofid(stripped, topology="hcb", linkage="boroxine")

            self.assertTrue(result.ok, result.reason)
            ledger = result.metadata.get("atom_ledger")
            self.assertIsNotNone(ledger)
            n_rings = int(result.metadata["n_events"])
            self.assertEqual(n_rings, 2)
            # Boroxine condensation: 3 ring oxygens per ring become residue and
            # two hydroxyl oxygens per boron are restored during repair
            # (net 3 H2O eliminated per ring).
            self.assertEqual(ledger["residue_atom_count"], 3 * n_rings)
            self.assertEqual(ledger["residue_elements"], {"O": 3 * n_rings})
            self.assertIn("boroxine", ledger["residue_policy"])
            self.assertEqual(ledger["reaction_additions"], {"O": 6 * n_rings})
            self.assertEqual(ledger["guest_molecule_count"], 0)
            self.assertEqual(ledger["unaccounted_atom_count"], 0)
            balance = ledger["balance"]
            self.assertTrue(balance["balanced"], balance)
            self.assertTrue(
                all(residual == 0 for residual in balance["per_element_residuals"].values())
            )


class EventAtomLedgerGuestTests(unittest.TestCase):
    def test_guest_molecule_and_atom_counts_are_distinct_and_balanced(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cif_path = _write_imine_cif(temp_path)
            guest_cif = _append_guest_atoms(
                cif_path,
                temp_path / "imine_guests.cif",
                (
                    "gst_O1 O 0.100000 0.100000 0.900000 1.00",
                    "gst_H1 H 0.151000 0.100000 0.900000 1.00",
                    "gst_H2 H 0.100000 0.151000 0.900000 1.00",
                    "gst_Na1 Na 0.600000 0.600000 0.900000 1.00",
                ),
                bond_rows=(
                    "gst_O1 gst_H1 . . 0.958 S",
                    "gst_O1 gst_H2 . . 0.958 S",
                ),
            )

            result = decompose_cif_to_cofid(guest_cif, topology="hcb", linkage="imine")

            self.assertTrue(result.ok, result.reason)
            ledger = result.metadata["atom_ledger"]
            # Two guest molecules (H2O, Na) with distinct molecule vs atom counts.
            self.assertEqual(ledger["guest_molecule_count"], 2)
            self.assertEqual(ledger["guest_atom_count"], 4)
            formulas = sorted(guest["molecular_formula"] for guest in ledger["guests"])
            self.assertEqual(formulas, ["H2O", "Na"])
            self.assertEqual(
                sorted(guest["atom_count"] for guest in ledger["guests"]), [1, 3]
            )
            # Guests enter the atom balance as their own bucket.
            balance = ledger["balance"]
            self.assertEqual(balance["guest_atom_count"], 4)
            self.assertTrue(balance["balanced"], balance)
            # Reaction stoichiometry is unaffected by guests.
            self.assertEqual(ledger["reaction_additions"], {"O": 3})

    def test_charged_guest_is_accounted_as_guest_localized(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cif_path = _write_imine_cif(temp_path)
            guest_cif = _append_guest_atoms(
                cif_path,
                temp_path / "imine_charged_guest.cif",
                ("gst_Na1 Na 0.600000 0.600000 0.900000 1.00",),
                charge_by_label={"gst_Na1": 1},
            )

            result = decompose_cif_to_cofid(guest_cif, topology="hcb", linkage="imine")

            self.assertTrue(result.ok, result.reason)
            ledger = result.metadata["atom_ledger"]
            charge_scope = ledger["charge_scope"]
            self.assertTrue(charge_scope["formal_charge_column_present"])
            self.assertEqual(charge_scope["input_net_formal_charge"], 1)
            self.assertEqual(charge_scope["input_charged_atom_count"], 1)
            self.assertEqual(charge_scope["guest_net_formal_charge"], 1)
            self.assertEqual(charge_scope["recovered_precursor_net_formal_charge"], 0)
            self.assertEqual(charge_scope["status"], CHARGE_SCOPE_GUEST_LOCALIZED)
            self.assertTrue(charge_scope["charge_conserved"])
            guest = ledger["guests"][0]
            self.assertEqual(guest["net_formal_charge"], 1)
            self.assertTrue(ledger["balance"]["balanced"])


class ChargeStereoScopeTests(unittest.TestCase):
    def _neutral_evidence(self) -> dict[str, object]:
        return {
            "formal_charge_column_present": True,
            "stereochemistry": {
                "perception_status": "perceived",
                "chiral_atom_count": 0,
                "stereo_bond_count": 0,
            },
        }

    def _ledger_for(self, mol, recovered_smiles: str):
        fragment = build_recovered_fragment_record(
            fragment_id=0,
            role="amine",
            mol=mol,
            atom_indices=tuple(range(mol.GetNumAtoms())),
            recovered_canonical_smiles=recovered_smiles,
        )
        return build_atom_ledger(
            family="imine",
            mol=mol,
            recovered_fragments=(fragment,),
            guests=(),
            residue_atom_indices=(),
            unaccounted=(),
            identity_evidence=self._neutral_evidence(),
        )

    def test_input_charge_preserved_in_recovered_precursor(self):
        mol = Chem.MolFromSmiles("C[NH3+]")
        self.assertIsNotNone(mol)

        ledger = self._ledger_for(mol, "C[NH3+]")

        self.assertEqual(ledger.charge_scope.status, CHARGE_SCOPE_PRESERVED)
        self.assertEqual(ledger.charge_scope.input_net_formal_charge, 1)
        self.assertEqual(ledger.charge_scope.recovered_precursor_net_formal_charge, 1)
        self.assertTrue(ledger.charge_scope.charge_conserved)
        self.assertTrue(ledger.balance.balanced)

    def test_charge_lost_in_recovery_is_explicitly_limited(self):
        mol = Chem.MolFromSmiles("C[NH3+]")
        self.assertIsNotNone(mol)

        ledger = self._ledger_for(mol, "CN")

        self.assertEqual(ledger.charge_scope.status, CHARGE_SCOPE_LIMITED)
        self.assertFalse(ledger.charge_scope.charge_conserved)
        self.assertNotEqual(ledger.charge_scope.net_charge_delta, 0)
        self.assertTrue(
            any("unsupported" in note for note in ledger.charge_scope.notes),
            ledger.charge_scope.notes,
        )

    def test_charge_created_by_repair_is_marked_heuristic(self):
        mol = Chem.MolFromSmiles("CN")
        self.assertIsNotNone(mol)

        ledger = self._ledger_for(mol, "C[NH3+]")

        self.assertEqual(ledger.charge_scope.status, CHARGE_SCOPE_RESTORED_BY_HEURISTIC)
        self.assertEqual(ledger.charge_scope.heuristic_restored_charge, 1)
        self.assertTrue(ledger.charge_scope.charge_conserved)

    def test_chiral_input_is_detected_and_marked_unsupported(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](F)Cl"))
        self.assertIsNotNone(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xC0F)
        conformer = mol.GetConformer()
        positions = tuple(
            tuple(float(value) for value in conformer.GetAtomPosition(index))
            for index in range(mol.GetNumAtoms())
        )

        evidence = perceive_stereochemistry_evidence(mol, positions)
        self.assertEqual(evidence["perception_status"], "perceived")
        self.assertEqual(evidence["chiral_atom_count"], 1)

        fragment = build_recovered_fragment_record(
            fragment_id=0,
            role="amine",
            mol=mol,
            atom_indices=tuple(range(mol.GetNumAtoms())),
            recovered_canonical_smiles="CC(F)Cl",
        )
        ledger = build_atom_ledger(
            family="imine",
            mol=mol,
            recovered_fragments=(fragment,),
            guests=(),
            residue_atom_indices=(),
            unaccounted=(),
            identity_evidence={"stereochemistry": evidence},
        )
        self.assertEqual(ledger.stereo_scope.status, STEREO_SCOPE_UNSUPPORTED)
        self.assertFalse(ledger.stereo_scope.preserved)
        self.assertIn("isomericSmiles=False", ledger.stereo_scope.policy)

    def test_build_bonded_mol_records_chiral_input_evidence(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](F)Cl"))
        self.assertIsNotNone(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xC0F)
        conformer = mol.GetConformer()
        cartesian = tuple(
            tuple(float(value) for value in conformer.GetAtomPosition(index))
            for index in range(mol.GetNumAtoms())
        )
        basis = ((20.0, 0.0, 0.0), (0.0, 20.0, 0.0), (0.0, 0.0, 20.0))
        atoms = PeriodicCifAtoms(
            symbols=tuple(atom.GetSymbol() for atom in mol.GetAtoms()),
            fractional_positions=tuple(
                tuple(coord / 20.0 for coord in position) for position in cartesian
            ),
            cartesian_positions=cartesian,
            cell_basis=basis,
            info={
                "_atom_site_label": tuple(
                    f"{atom.GetSymbol()}{index}"
                    for index, atom in enumerate(mol.GetAtoms(), start=1)
                )
            },
        )

        build_result = _build_bonded_mol(atoms, bond_mode="distance")

        evidence = build_result.metadata["input_identity_evidence"]["stereochemistry"]
        self.assertEqual(evidence["perception_status"], "perceived")
        self.assertEqual(evidence["chiral_atom_count"], 1)
        self.assertFalse(
            build_result.metadata["input_identity_evidence"]["formal_charge_column_present"]
        )


class UnaccountedAtomLedgerTests(unittest.TestCase):
    def _boroxine_ring_with_stray_heteroatom(self):
        mol = Chem.MolFromSmiles("c1ccc(B2SB(c3ccccc3)OB(c3ccccc3)O2)cc1")
        self.assertIsNotNone(mol)
        Chem.GetSymmSSSR(mol)
        ring = next(
            tuple(int(atom_idx) for atom_idx in ring_atoms)
            for ring_atoms in mol.GetRingInfo().AtomRings()
            if len(ring_atoms) == 6
            and sum(1 for atom_idx in ring_atoms if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 5)
            == 3
        )
        ring_bonds = []
        for index, atom_idx in enumerate(ring):
            bond = mol.GetBondBetweenAtoms(atom_idx, ring[(index + 1) % 6])
            ring_bonds.append(int(bond.GetIdx()))
        anchors = tuple(
            atom_idx
            for atom_idx in ring
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 5
        )
        event = LinkageEvent(
            event_id="boroxine:ring:test",
            family="boroxine",
            atoms=ring,
            bonds=tuple(ring_bonds),
            cut_bonds=tuple(ring_bonds),
            instance_ids=tuple(None for _ in ring),
            confidence="high",
            endpoint_roles=tuple((atom_idx, "boronic_acid") for atom_idx in anchors),
            site_id="boroxine:ring:test",
            metadata={
                "ring_atoms": list(ring),
                "anchor_atom_indices": list(anchors),
                "anchor_images": [[0, 0, 0] for _ in anchors],
                "ring_size": 6,
            },
        )
        return mol, event

    def test_stray_framework_atom_is_explicitly_unaccounted(self):
        mol, event = self._boroxine_ring_with_stray_heteroatom()
        build_result = BondedMolBuildResult(mol=mol, metadata={}, candidates=())

        cut_result = _cut_and_reconstruct(build_result, (event,), _BOROXINE_SPEC)

        # The sulfur ring member is not allowed boroxine residue (oxygen only),
        # so it must surface as an explicit unaccounted atom, never vanish.
        self.assertEqual(cut_result.status, EVENT_STATUS_UNEXPLAINED)
        self.assertIsNotNone(cut_result.atom_ledger)
        ledger = cut_result.atom_ledger
        assert ledger is not None
        self.assertEqual(ledger.unaccounted_atom_count, 1)
        self.assertEqual(ledger.unaccounted[0].elements, {"S": 1})
        self.assertIn("not assigned", ledger.unaccounted[0].reason)
        self.assertEqual(ledger.residue_atom_count, 2)
        self.assertEqual(ledger.residue_elements, {"O": 2})
        self.assertTrue(ledger.balance.balanced, ledger.balance)

        # The ledger reaches the selection-level result metadata on failure.
        hypothesis = ReconstructionHypothesis(
            hypothesis_id="boroxine:test",
            events=(event,),
            monomers=cut_result.monomers,
            status=cut_result.status,
            validation_errors=cut_result.errors,
            metadata=cut_result.metadata,
            atom_ledger=ledger,
        )
        result = _select_event_result(
            Path("synthetic.cif"),
            requested_family=None,
            topology=None,
            detection=EventDetectionResult(events=(event,)),
            hypotheses=(hypothesis,),
            generation_metadata={},
        )
        self.assertEqual(result.metadata["event_status"], EVENT_STATUS_UNEXPLAINED)
        result_ledger = result.metadata["atom_ledger"]
        self.assertEqual(result_ledger["unaccounted_atom_count"], 1)
        self.assertEqual(result_ledger["unaccounted"][0]["elements"], {"S": 1})
        json.dumps(result.to_dict())


if __name__ == "__main__":
    unittest.main()
