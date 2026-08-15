import contextlib
import io
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator, CoarseValidationThresholds
from cofkit.build_workflows.ring_forming import RingFormationConfig, RingFormingStructureGenerator
from cofkit.bond_types import bond_order_to_cif_type
from cofkit.chem.rdkit import build_rdkit_monomer
from cofkit.cif import CIFWriter
from cofkit.cli import main as cli_main
from cofkit.cofid import generate_cofid
from cofkit.decompose import (
    BondedMolBuildResult,
    DecomposedMonomer,
    _IMINE_SPEC,
    _classify_nitrogen_linkage_bonds,
    _classify_vinylene_linkage_bonds,
    _explicit_bond_candidates,
    _finalize_repaired_fragment,
    _mark_imine_linkage_bonds,
    _minimum_image_bond_geometry,
    _periodic_dimension_hint,
    _periodic_edge_gains_match,
    _validate_recovered_precursors,
    _validate_supported_cif_periodicity,
    decompose_cif_to_cofid,
)
from cofkit.decompose_events import EVENT_STATUS_TRIAZINE_MOTIF, detect_linkage_events
from cofkit.decompose_cif import PeriodicCifAtoms
from cofkit.reactions import ReactionLibrary
from cofkit.validate import validate_cif_against_cofid


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
PDA = "O=Cc1ccc(C=O)cc1"
COF42_HYDRAZIDE = "CCOc1cc(C(=O)NN)cc(C(=O)NN)c1OCC"
HYDRAZINE = "NN"
TP = "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O"
PPD = "Nc1ccc(N)cc1"
BDBA = "OB(O)c1ccc(B(O)O)cc1"
HHTP = "OC1=C(O)C=C2C(=C1)C1=CC(O)=C(O)C=C1C1=CC(O)=C(O)C=C21"
TMT = "Cc1nc(C)nc(C)n1"
MELAMINE = "Nc1nc(N)nc(N)n1"
BEX_D2H_ALDEHYDE = (
    "C1C=C(N(C2C=CC(C3C4C(=NON=4)C(C4C=CC(N(C5C=CC(C=O)=CC=5)C5C=CC(C=O)=CC=5)=CC=4)=CC=3)=CC=2)"
    "C2C=CC(C=O)=CC=2)C=CC=1C=O"
)
BEX_D2H_AMINE = (
    "Nc1ccc(-c2ccc3c(c2)C(=C2c4cc(-c5ccc(N)cc5)ccc4-c4ccc(-c5ccc(N)cc5)cc42)c2cc(-c4ccc(N)cc4)ccc2-3)cc1"
)
TEREPHTHALALDEHYDE = "O=Cc1ccc(C=O)cc1"
TETRA_AMINE = "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1"
TETRA_ALDEHYDE = "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1"
HEXA_AMINE = "Nc1cc2c3cc(N)c(N)cc3c3cc(N)c(N)cc3c2cc1N"


def _generator(template_id: str = "imine_bridge") -> BatchStructureGenerator:
    return BatchStructureGenerator(
        BatchGenerationConfig(
            allowed_reactions=(template_id,),
            single_node_topology_ids=("hcb",),
            topology_ids=("hcb",),
            enumerate_all_topologies=False,
            write_cif=True,
            rdkit_num_conformers=2,
        )
    )


def _record(
    record_id: str,
    smiles: str,
    motif_kind: str,
    connectivity: int,
) -> BatchMonomerRecord:
    return BatchMonomerRecord(
        id=record_id,
        name=record_id,
        smiles=smiles,
        motif_kind=motif_kind,
        expected_connectivity=connectivity,
    )


ALL_BINARY_LINKAGE_ROUND_TRIP_CASES = (
    (
        "hydrazone_bridge",
        "hydrazone",
        _record("hdz", COF42_HYDRAZIDE, "hydrazide", 2),
        _record("tfb", TFB, "aldehyde", 3),
        6,
    ),
    (
        "azine_bridge",
        "azine",
        _record("hyd", HYDRAZINE, "hydrazine", 2),
        _record("tfb", TFB, "aldehyde", 3),
        6,
    ),
    (
        "keto_enamine_bridge",
        "bken",
        _record("tp", TP, "keto_aldehyde", 3),
        _record("ppd", PPD, "amine", 2),
        6,
    ),
    (
        "boronate_ester_bridge",
        "boest",
        _record("bdba", BDBA, "boronic_acid", 2),
        _record("hhtp", HHTP, "catechol", 3),
        12,
    ),
    (
        "vinylene_bridge",
        "vinylene",
        _record("tmt", TMT, "activated_methylene", 3),
        _record("pda", PDA, "aldehyde", 2),
        6,
    ),
)
ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES = (
    (
        "imine_bridge",
        "imine",
        _record("tapb", TAPB, "amine", 3),
        _record("tfb", TFB, "aldehyde", 3),
        3,
    ),
) + ALL_BINARY_LINKAGE_ROUND_TRIP_CASES

NITROGEN_LINKAGE_ROUND_TRIP_CASES = (
    ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES[0],
    ALL_BINARY_LINKAGE_ROUND_TRIP_CASES[0],
    ALL_BINARY_LINKAGE_ROUND_TRIP_CASES[1],
    ALL_BINARY_LINKAGE_ROUND_TRIP_CASES[2],
)

TOPOLOGY_ENABLED_ROUND_TRIP_POOL = (
    (
        "hcb_3_3_imine",
        "imine_bridge",
        "hcb",
        _record("tapb", TAPB, "amine", 3),
        _record("tfb", TFB, "aldehyde", 3),
        False,
    ),
    (
        "hcb_3_2_imine",
        "imine_bridge",
        "hcb",
        _record("tapb", TAPB, "amine", 3),
        _record("pda", PDA, "aldehyde", 2),
        False,
    ),
    (
        "hcb_2_3_hydrazone",
        "hydrazone_bridge",
        "hcb",
        _record("hdz", COF42_HYDRAZIDE, "hydrazide", 2),
        _record("tfb", TFB, "aldehyde", 3),
        False,
    ),
    (
        "hcb_3_2_bken",
        "keto_enamine_bridge",
        "hcb",
        _record("tp", TP, "keto_aldehyde", 3),
        _record("ppd", PPD, "amine", 2),
        False,
    ),
    (
        "sql_4_2_imine",
        "imine_bridge",
        "sql",
        _record("tetra_amine", TETRA_AMINE, "amine", 4),
        _record("tpal", TEREPHTHALALDEHYDE, "aldehyde", 2),
        True,
    ),
    (
        "sql_4_4_imine",
        "imine_bridge",
        "sql",
        _record("tetra_amine", TETRA_AMINE, "amine", 4),
        _record("tetra_aldehyde", TETRA_ALDEHYDE, "aldehyde", 4),
        True,
    ),
    (
        "hxl_6_2_imine",
        "imine_bridge",
        "hxl",
        _record("hexa_amine", HEXA_AMINE, "amine", 6),
        _record("tpal", TEREPHTHALALDEHYDE, "aldehyde", 2),
        False,
    ),
    (
        "bex_4_4_imine",
        "imine_bridge",
        "bex",
        _record("bex_d2h_amine", BEX_D2H_AMINE, "amine", 4),
        _record("bex_d2h_aldehyde", BEX_D2H_ALDEHYDE, "aldehyde", 4),
        True,
    ),
)


def _cif_lines_without_cofid(source_path: str | Path) -> list[str]:
    lines = Path(source_path).read_text(encoding="utf-8").splitlines()
    if lines and lines[0].startswith("# COFid: "):
        lines = lines[1:]
    return lines


def _write_cif_lines(lines: list[str], target_path: Path) -> Path:
    target_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target_path


def _without_cofid_comment(source_path: str | Path, target_path: Path) -> Path:
    return _write_cif_lines(_cif_lines_without_cofid(source_path), target_path)


def _without_cif_bond_type_column(source_path: str | Path, target_path: Path) -> Path:
    lines = _cif_lines_without_cofid(source_path)
    output: list[str] = []
    index = 0
    while index < len(lines):
        if lines[index].strip() == "loop_":
            next_index = index + 1
            headers: list[str] = []
            while next_index < len(lines) and lines[next_index].lstrip().startswith("_"):
                headers.append(lines[next_index].strip())
                next_index += 1
            if "_geom_bond_atom_site_label_1" in headers and "_geom_bond_atom_site_label_2" in headers:
                remove_indices = {
                    header_index
                    for header_index, header in enumerate(headers)
                    if header in {"_ccdc_geom_bond_type", "_geom_bond_type"}
                }
                output.append(lines[index])
                output.extend(
                    header for header_index, header in enumerate(headers) if header_index not in remove_indices
                )
                while (
                    next_index < len(lines)
                    and lines[next_index].strip()
                    and lines[next_index].strip() != "loop_"
                    and not lines[next_index].lstrip().startswith("_")
                ):
                    row = lines[next_index].split()
                    output.append(
                        " ".join(value for value_index, value in enumerate(row) if value_index not in remove_indices)
                    )
                    next_index += 1
                index = next_index
                continue
        output.append(lines[index])
        index += 1
    return _write_cif_lines(output, target_path)


def _without_cif_bond_loop(source_path: str | Path, target_path: Path) -> Path:
    lines = _cif_lines_without_cofid(source_path)
    output: list[str] = []
    index = 0
    while index < len(lines):
        if lines[index].strip() == "loop_":
            next_index = index + 1
            headers: list[str] = []
            while next_index < len(lines) and lines[next_index].lstrip().startswith("_"):
                headers.append(lines[next_index].strip())
                next_index += 1
            if "_geom_bond_atom_site_label_1" in headers and "_geom_bond_atom_site_label_2" in headers:
                while (
                    next_index < len(lines)
                    and lines[next_index].strip()
                    and lines[next_index].strip() != "loop_"
                    and not lines[next_index].lstrip().startswith("_")
                ):
                    next_index += 1
                index = next_index
                continue
        output.append(lines[index])
        index += 1
    return _write_cif_lines(output, target_path)


def _with_generic_atom_labels(source_path: str | Path, target_path: Path) -> Path:
    lines = Path(source_path).read_text(encoding="utf-8").splitlines()
    output: list[str] = []
    label_map: dict[str, str] = {}
    index = 0
    while index < len(lines):
        if lines[index].strip() == "loop_":
            next_index = index + 1
            headers: list[str] = []
            while next_index < len(lines) and lines[next_index].lstrip().startswith("_"):
                headers.append(lines[next_index].strip())
                next_index += 1
            if "_atom_site_label" in headers:
                label_index = headers.index("_atom_site_label")
                type_index = headers.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in headers else None
                output.append(lines[index])
                output.extend(headers)
                atom_index = 1
                while (
                    next_index < len(lines)
                    and lines[next_index].strip()
                    and lines[next_index].strip() != "loop_"
                    and not lines[next_index].lstrip().startswith("_")
                ):
                    row = lines[next_index].split()
                    symbol = row[type_index] if type_index is not None else "X"
                    old_label = row[label_index]
                    row[label_index] = f"{symbol}{atom_index}"
                    label_map[old_label] = row[label_index]
                    output.append(" ".join(row))
                    atom_index += 1
                    next_index += 1
                index = next_index
                continue
            if "_geom_bond_atom_site_label_1" in headers and "_geom_bond_atom_site_label_2" in headers:
                first_label_index = headers.index("_geom_bond_atom_site_label_1")
                second_label_index = headers.index("_geom_bond_atom_site_label_2")
                output.append(lines[index])
                output.extend(headers)
                while (
                    next_index < len(lines)
                    and lines[next_index].strip()
                    and lines[next_index].strip() != "loop_"
                    and not lines[next_index].lstrip().startswith("_")
                ):
                    row = lines[next_index].split()
                    row[first_label_index] = label_map.get(row[first_label_index], row[first_label_index])
                    row[second_label_index] = label_map.get(row[second_label_index], row[second_label_index])
                    output.append(" ".join(row))
                    next_index += 1
                index = next_index
                continue
        output.append(lines[index])
        index += 1
    return _write_cif_lines(output, target_path)


def _write_molecular_probe_cif(smiles: str, target_path: Path) -> Path:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid molecular probe SMILES: {smiles!r}")
    rdDepictor.Compute2DCoords(mol)
    conformer = mol.GetConformer()
    labels = tuple(
        f"{atom.GetSymbol()}{atom.GetIdx() + 1}"
        for atom in mol.GetAtoms()
    )
    lines = [
        "data_molecular_probe",
        "_space_group_name_H-M_alt 'P 1'",
        "_space_group_IT_number 1",
        "_cell_length_a 30.0",
        "_cell_length_b 30.0",
        "_cell_length_c 30.0",
        "_cell_angle_alpha 90.0",
        "_cell_angle_beta 90.0",
        "_cell_angle_gamma 90.0",
        "",
        "loop_",
        "_space_group_symop_operation_xyz",
        "'x,y,z'",
        "",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
    ]
    for atom, label in zip(mol.GetAtoms(), labels):
        position = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{label} {atom.GetSymbol()} "
            f"{(position.x + 10.0) / 30.0:.8f} "
            f"{(position.y + 10.0) / 30.0:.8f} 0.33333333 1.00"
        )
    lines.extend([
        "",
        "loop_",
        "_geom_bond_atom_site_label_1",
        "_geom_bond_atom_site_label_2",
        "_geom_bond_distance",
        "_geom_bond_site_symmetry_1",
        "_geom_bond_site_symmetry_2",
        "_ccdc_geom_bond_type",
    ])
    for bond in mol.GetBonds():
        first = conformer.GetAtomPosition(bond.GetBeginAtomIdx())
        second = conformer.GetAtomPosition(bond.GetEndAtomIdx())
        distance = ((first.x - second.x) ** 2 + (first.y - second.y) ** 2) ** 0.5
        lines.append(
            f"{labels[bond.GetBeginAtomIdx()]} {labels[bond.GetEndAtomIdx()]} "
            f"{distance:.6f} . . {bond_order_to_cif_type(bond.GetBondTypeAsDouble())}"
        )
    return _write_cif_lines(lines, target_path)


class DecomposeRoundTripTests(unittest.TestCase):
    def test_precursor_validation_rejects_incomplete_or_same_role_binary_output(self):
        amine = DecomposedMonomer(1, "amine", "CN")
        aldehyde = DecomposedMonomer(1, "aldehyde", "CC=O")

        for monomers in ((amine,), (aldehyde, aldehyde), (amine, aldehyde, aldehyde)):
            with self.subTest(roles=tuple(monomer.reactive_group for monomer in monomers)):
                reason, metadata = _validate_recovered_precursors(monomers, _IMINE_SPEC)
                self.assertIsNotNone(reason)
                self.assertEqual(metadata["status"], "incomplete")

        reason, metadata = _validate_recovered_precursors((amine, aldehyde), _IMINE_SPEC)
        self.assertIsNone(reason)
        self.assertEqual(metadata["status"], "valid")

        invalid_amine = DecomposedMonomer(1, "amine", "N")
        reason, metadata = _validate_recovered_precursors((invalid_amine, aldehyde), _IMINE_SPEC)
        self.assertIn("connectivity mismatch", reason)
        self.assertEqual(metadata["status"], "chemically_invalid")

    def test_nitrogen_cross_branch_conflict_only_withholds_conflicting_bond(self):
        mol = Chem.MolFromSmiles("CC=NN=CC1C(=O)C=CC=C1.CC=NC")
        self.assertIsNotNone(mol)
        for atom in mol.GetAtoms():
            atom.SetProp("instance_id", "")

        classification = _classify_nitrogen_linkage_bonds(mol)
        marked = _mark_imine_linkage_bonds(mol)

        self.assertEqual(len(classification.cross_branch_conflict_bond_indices), 1)
        self.assertEqual(marked, classification.imine_bond_indices)
        self.assertEqual(len(marked), 1)

    def test_minimum_image_search_handles_skew_cell(self):
        basis = (
            (3.0, 0.0, 0.0),
            (-1.5, 2.598076211, 0.0),
            (0.0, 0.0, 10.0),
        )
        fractional = ((0.0, 0.49, 0.0), (0.49, 0.0, 0.0))

        def cartesian(position):
            return tuple(
                sum(position[axis] * basis[axis][component] for axis in range(3))
                for component in range(3)
            )

        atoms = PeriodicCifAtoms(
            symbols=("C", "C"),
            fractional_positions=fractional,
            cartesian_positions=tuple(cartesian(position) for position in fractional),
            cell_basis=basis,
            info={"_atom_site_label": ("C1", "C2")},
        )

        image, distance = _minimum_image_bond_geometry(atoms, 0, 1)

        self.assertEqual(image, (0, 1, 0))
        self.assertAlmostEqual(distance, 1.50089973, places=6)

    def test_explicit_bond_candidates_preserve_distinct_periodic_images(self):
        basis = ((10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0))
        atoms = PeriodicCifAtoms(
            symbols=("C", "N"),
            fractional_positions=((0.0, 0.0, 0.0), (0.13, 0.0, 0.0)),
            cartesian_positions=((0.0, 0.0, 0.0), (1.3, 0.0, 0.0)),
            cell_basis=basis,
            info={
                "_geom_bond_distance": ("1.3", "1.3", "1.3"),
                "_ccdc_geom_bond_type": ("D", "D", "D"),
                "_geom_bond_site_symmetry_1": (".", ".", "."),
                "_geom_bond_site_symmetry_2": (".", "1_655", "1_455"),
            },
        )

        candidates, _missing = _explicit_bond_candidates(
            atoms,
            {"C1": 0, "N1": 1},
            ("C1", "C1", "N1"),
            ("N1", "N1", "C1"),
        )

        self.assertEqual(len(candidates), 2)
        self.assertEqual({candidate.periodic_image for candidate in candidates}, {(0, 0, 0), (1, 0, 0)})

    def test_non_p1_periodicity_is_rejected_explicitly(self):
        atoms = PeriodicCifAtoms(
            symbols=("C",),
            fractional_positions=((0.1, 0.2, 0.3),),
            cartesian_positions=((1.0, 2.0, 3.0),),
            cell_basis=((10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)),
            info={"_space_group_symop_operation_xyz": ("x,y,z", "-x,-y,-z")},
        )

        with self.assertRaisesRegex(ValueError, "P1-expanded"):
            _validate_supported_cif_periodicity(atoms)

    def test_fragment_finalization_preserves_formal_charge(self):
        fragment = Chem.MolFromSmiles("[NH4+]")
        self.assertIsNotNone(fragment)
        fragment.GetAtomWithIdx(0).SetProp("cofkit_decompose_role", "amine")

        monomer = _finalize_repaired_fragment(fragment, "amine")

        self.assertIsNotNone(monomer)
        self.assertEqual(monomer.canonical_smiles, "[NH4+]")

    def test_periodic_dimension_uses_cycle_gain_rank_not_z_axis(self):
        rotated_2d_edges = (
            (0, 0, (0, 1, 0)),
            (0, 0, (0, 0, 1)),
        )
        three_dimensional_edges = rotated_2d_edges + ((0, 0, (1, 0, 0)),)

        self.assertEqual(_periodic_dimension_hint(1, rotated_2d_edges), "2D")
        self.assertEqual(_periodic_dimension_hint(1, three_dimensional_edges), "3D")
        self.assertIsNone(_periodic_dimension_hint(1, ()))
        self.assertIsNone(_periodic_dimension_hint(1, ((0, 0, (1, 0, 0)),)))

    def test_molecular_motifs_cannot_masquerade_as_periodic_frameworks(self):
        cases = (
            ("CC=NC", "imine"),
            ("CC(=O)C=Cc1ccccc1", "vinylene"),
            ("B1(OCCO1)c1ccccc1", "boest"),
            (TMT, "triazine"),
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            for index, (smiles, linkage) in enumerate(cases):
                with self.subTest(linkage=linkage):
                    cif_path = _write_molecular_probe_cif(
                        smiles,
                        temp_path / f"probe_{index}.cif",
                    )

                    result = decompose_cif_to_cofid(
                        cif_path,
                        topology="hcb",
                        linkage=linkage,
                    )

                    self.assertFalse(result.ok)
                    self.assertEqual(result.status, "skipped")
                    if linkage in {"imine", "triazine"}:
                        self.assertEqual(
                            result.metadata["topology_validation"]["observed_periodic_rank"],
                            0,
                        )
                    else:
                        self.assertEqual(
                            result.metadata["precursor_recovery_validation"]["status"],
                            "chemically_invalid",
                        )

    def test_conventional_beta_ketoenamine_bond_orders_are_classified_as_bken(self):
        mol = Chem.MolFromSmiles("CNC=C1C(=O)C=CC=C1")
        self.assertIsNotNone(mol)
        for atom in mol.GetAtoms():
            atom.SetProp("instance_id", "")

        metadata = _classify_nitrogen_linkage_bonds(mol).to_metadata()

        self.assertEqual(metadata["beta_ketoenamine_single_bond_count"], 1)
        self.assertEqual(metadata["assigned_bond_counts"]["bken"], 1)
        self.assertEqual(metadata["assigned_bond_counts"]["imine"], 0)

    def test_periodic_gain_matching_allows_vertex_switching_and_axis_permutation(self):
        expected = (
            (0, 1, (0, 0, 0)),
            (0, 1, (1, 0, 0)),
            (0, 1, (0, 1, 0)),
        )
        observed = (
            (0, 1, (3, -2, 0)),
            (0, 1, (3, -1, 0)),
            (0, 1, (3, -2, 1)),
        )
        sheared = (
            (0, 1, (0, 0, 0)),
            (0, 1, (1, 0, 0)),
            (0, 1, (-1, 1, 0)),
        )
        expected_simple = (
            (0, 1, (0, 0, 0)),
            (0, 2, (0, 0, 0)),
            (0, 3, (0, 0, 0)),
            (1, 2, (1, 0, 0)),
            (1, 3, (0, 1, 0)),
            (2, 3, (1, 1, 0)),
        )
        sheared_simple = (
            (0, 1, (0, 0, 0)),
            (0, 2, (0, 0, 0)),
            (0, 3, (0, 0, 0)),
            (1, 2, (1, 0, 0)),
            (1, 3, (-1, 1, 0)),
            (2, 3, (0, 1, 0)),
        )
        non_unimodular = (
            (0, 1, (0, 0, 0)),
            (0, 1, (2, 0, 0)),
            (0, 1, (0, 1, 0)),
        )

        self.assertTrue(_periodic_edge_gains_match(observed, expected))
        self.assertTrue(_periodic_edge_gains_match(sheared, expected))
        self.assertTrue(_periodic_edge_gains_match(sheared_simple, expected_simple))
        self.assertFalse(_periodic_edge_gains_match(non_unimodular, expected))

    def test_decompose_normalizes_known_topology_and_rejects_invalid_token(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            normalized = decompose_cif_to_cofid(input_cif, topology=" HCB ")
            invalid = decompose_cif_to_cofid(input_cif, topology="bad&&token")

        self.assertTrue(normalized.ok, normalized.reason)
        self.assertEqual(normalized.topology, "hcb")
        self.assertEqual(invalid.status, "error")
        self.assertFalse(invalid.ok)

    def test_supplied_topology_must_match_recovered_graph_and_connectivities(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(
                input_cif,
                topology="sql",
                linkage="imine",
            )

        self.assertFalse(result.ok)
        self.assertEqual(result.metadata["topology_validation"]["status"], "incompatible")

    def test_incompatible_embedded_topology_falls_back_to_periodic_graph(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            lines = Path(summary.cif_path).read_text(encoding="utf-8").splitlines()
            self.assertTrue(lines[0].startswith("# COFid: "))
            lines[0] = lines[0].replace("&&hcb&&", "&&sql&&")
            input_cif = _write_cif_lines(lines, temp_path / "wrong_comment.cif")

            result = decompose_cif_to_cofid(input_cif, linkage="imine")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.topology, "hcb")
        detection = result.metadata["topology_detection"]
        self.assertEqual(detection["metadata"]["source"], "periodic_linkage_graph")
        self.assertEqual(
            detection["metadata"]["embedded_cofid_comment"]["status"],
            "rejected",
        )

    def test_decompose_generated_hcb_three_three_imine_returns_original_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 3)
        self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_decompose_auto_detects_topology_from_cofid_comment(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )

            result = decompose_cif_to_cofid(summary.cif_path)

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.topology, "hcb")
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        detection = result.metadata["topology_detection"]
        self.assertEqual(detection["confidence"], "exact")
        self.assertEqual(detection["metadata"]["source"], "cofid_comment")

    def test_decompose_generated_hcb_three_two_imine_returns_original_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("pda", PDA, "aldehyde", 2),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 6)
        self.assertEqual(result.metadata["n_unique_monomers"], 2)
        self.assertEqual(result.metadata["topology_graph"]["periodic_rank"], 2)
        self.assertEqual(result.metadata["linkage_detection"]["successful_linkages"], ["imine"])

    def test_decompose_generated_decorated_bex_imine_recovers_monomers(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                topology_ids=("bex",),
                stacking_ids=("AA",),
                enumerate_all_topologies=False,
                write_cif=True,
                rdkit_num_conformers=1,
                hard_hard_max_bridge_distance=10.0,
                validation_thresholds=CoarseValidationThresholds(hard_hard_max_bridge_distance=10.0),
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = generator.generate_pair_candidate(
                _record("bex_d2h_amine", BEX_D2H_AMINE, "amine", 4),
                _record("bex_d2h_aldehyde", BEX_D2H_ALDEHYDE, "aldehyde", 4),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped_bex.cif")

            result = decompose_cif_to_cofid(input_cif, topology="bex")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 8)
        self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_analyze_decompose_cli_prints_recovered_cofid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_main(["analyze", "decompose", str(input_cif), "--topology", "hcb"])

        self.assertEqual(buffer.getvalue().strip(), summary.metadata["cofid"])

    def test_analyze_decompose_cli_auto_detects_topology_when_omitted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_main(["analyze", "decompose", str(input_cif)])

        self.assertEqual(buffer.getvalue().strip(), summary.metadata["cofid"])

    def test_analyze_decompose_cli_exposes_event_mode_json(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_main(
                    [
                        "analyze",
                        "decompose",
                        str(input_cif),
                        "--topology",
                        "hcb",
                        "--decomposition-mode",
                        "event",
                        "--json",
                    ]
                )

        payload = json.loads(buffer.getvalue())
        self.assertEqual(payload["cofid"], summary.metadata["cofid"])
        self.assertEqual(payload["metadata"]["decomposition_mode"], "event")
        self.assertEqual(payload["metadata"]["event_status"], "SUCCESS_COMPLETE")

    def test_decompose_generated_hcb_cofs_for_all_build_binary_linkages(self):
        for template_id, linkage, first, second, expected_linkage_bonds in ALL_BINARY_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

                result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage=linkage)

                self.assertTrue(result.ok, result.reason)
                self.assertEqual(result.cofid, summary.metadata["cofid"])
                self.assertEqual(result.metadata[f"n_{linkage}_linkage_bonds"], expected_linkage_bonds)
                self.assertEqual(result.metadata["n_unique_monomers"], 2)
                if linkage == "vinylene":
                    detection = result.metadata["vinylene_linkage_detection"]
                    self.assertEqual(detection["accepted_bond_count"], expected_linkage_bonds)
                    self.assertEqual(detection["small_ring_rejected_bond_count"], 0)
                    self.assertEqual(detection["endpoint_rejected_bond_count"], 0)

    def test_event_mode_round_trips_all_buildable_binary_linkages(self):
        for template_id, linkage, first, second, expected_linkage_bonds in ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

                legacy_result = decompose_cif_to_cofid(
                    input_cif,
                    topology="hcb",
                    linkage=linkage,
                )
                event_result = decompose_cif_to_cofid(
                    input_cif,
                    topology="hcb",
                    linkage=linkage,
                    decomposition_mode="event",
                )
                automatic_event_result = decompose_cif_to_cofid(
                    input_cif,
                    topology="hcb",
                    decomposition_mode="event",
                )

                self.assertTrue(legacy_result.ok, legacy_result.reason)
                self.assertNotIn("decomposition_mode", legacy_result.metadata)
                self.assertTrue(event_result.ok, event_result.reason)
                self.assertEqual(event_result.cofid, summary.metadata["cofid"])
                self.assertEqual(event_result.metadata["decomposition_mode"], "event")
                self.assertEqual(event_result.metadata["event_status"], "SUCCESS_COMPLETE")
                self.assertTrue(automatic_event_result.ok, automatic_event_result.reason)
                self.assertEqual(automatic_event_result.linkage, linkage)
                self.assertEqual(automatic_event_result.cofid, summary.metadata["cofid"])
                family_events = [
                    event
                    for event in event_result.metadata["event_detection"]["events"]
                    if event["family"] == linkage
                ]
                expected_event_count = (
                    expected_linkage_bonds // 2
                    if linkage in {"azine", "boest"}
                    else expected_linkage_bonds
                )
                self.assertEqual(len(family_events), expected_event_count)
                if linkage in {"azine", "boest"}:
                    self.assertTrue(all(len(event["cut_bonds"]) == 2 for event in family_events))

    def test_event_detector_branches_tied_vinylene_orientations(self):
        mol = Chem.MolFromSmiles("CC=CC")
        self.assertIsNotNone(mol)
        for atom in mol.GetAtoms():
            atom.SetProp("instance_id", "")

        detection = detect_linkage_events(BondedMolBuildResult(mol=mol))
        vinylene_events = [event for event in detection.events if event.family == "vinylene"]

        self.assertEqual(len(vinylene_events), 2)
        self.assertEqual(len({event.site_id for event in vinylene_events}), 1)
        self.assertEqual({event.confidence for event in vinylene_events}, {"low"})
        self.assertEqual(
            detection.diagnostics["detectors"]["vinylene"]["ambiguous_orientation_site_count"],
            1,
        )

    def test_event_mode_globally_classifies_triazine_as_an_intact_imine_monomer_motif(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator("imine_bridge").generate_pair_candidate(
                _record("melamine", MELAMINE, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            automatic = decompose_cif_to_cofid(
                input_cif,
                topology="hcb",
                decomposition_mode="event",
            )
            triazine = decompose_cif_to_cofid(
                input_cif,
                topology="hcb",
                linkage="triazine",
                decomposition_mode="event",
            )

        self.assertTrue(automatic.ok, automatic.reason)
        self.assertEqual(automatic.linkage, "imine")
        self.assertFalse(triazine.ok)
        self.assertEqual(triazine.metadata["event_status"], EVENT_STATUS_TRIAZINE_MOTIF)
        self.assertIn("intact monomer motif", triazine.reason)

    def test_event_mode_is_explicit_and_invalid_modes_are_rejected(self):
        result = decompose_cif_to_cofid(
            "unused.cif",
            decomposition_mode="not-a-mode",
        )

        self.assertEqual(result.status, "error")
        self.assertIn("legacy", result.reason)
        self.assertIn("event", result.reason)

    def test_nitrogen_linkage_priority_branches_make_generated_matches_exclusive(self):
        requested_linkages = ("azine", "hydrazone", "bken", "imine")
        for template_id, expected_linkage, first, second, expected_linkage_bonds in NITROGEN_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

                for requested_linkage in requested_linkages:
                    result = decompose_cif_to_cofid(
                        input_cif,
                        topology="hcb",
                        linkage=requested_linkage,
                    )
                    detection = result.metadata["nitrogen_linkage_detection"]
                    assigned_counts = detection["assigned_bond_counts"]
                    self.assertEqual(detection["resolution_status"], "resolved")
                    self.assertEqual(assigned_counts[expected_linkage], expected_linkage_bonds)
                    self.assertEqual(sum(assigned_counts.values()), expected_linkage_bonds)
                    if requested_linkage == expected_linkage:
                        self.assertTrue(result.ok, result.reason)
                        self.assertEqual(result.cofid, summary.metadata["cofid"])
                    else:
                        self.assertFalse(result.ok)
                        self.assertEqual(result.metadata[f"n_{requested_linkage}_linkage_bonds"], 0)

    def test_nitrogen_linkage_classifier_reports_cross_branch_conflict(self):
        mol = Chem.MolFromSmiles("CC=NN=CC1C(=O)C=CC=C1")
        self.assertIsNotNone(mol)
        for atom in mol.GetAtoms():
            atom.SetProp("instance_id", "")

        classification = _classify_nitrogen_linkage_bonds(mol)
        metadata = classification.to_metadata()

        self.assertEqual(metadata["resolution_status"], "cross_branch_ambiguous")
        self.assertEqual(metadata["cross_branch_conflict_bond_count"], 1)

    def test_specific_nitrogen_motif_does_not_promote_separate_imine_bond(self):
        mol = Chem.MolFromSmiles("CC=NN=CC.CC=NC")
        self.assertIsNotNone(mol)
        for atom in mol.GetAtoms():
            atom.SetProp("instance_id", "")

        metadata = _classify_nitrogen_linkage_bonds(mol).to_metadata()

        self.assertEqual(metadata["resolution_status"], "resolved")
        self.assertEqual(
            metadata["assigned_bond_counts"],
            {"azine": 2, "hydrazone": 0, "bken": 0, "imine": 1},
        )

    def test_vinylene_classifier_rejects_small_ring_and_unactivated_double_bonds(self):
        benzene = Chem.MolFromSmiles("c1ccccc1")
        Chem.Kekulize(benzene, clearAromaticFlags=True)
        cyclohexene = Chem.MolFromSmiles("C1=CCCCC1")
        unactivated = Chem.MolFromSmiles("CC=CC")
        valid = Chem.MolFromSmiles("CC(=O)C=Cc1ccccc1")
        for mol in (benzene, cyclohexene, unactivated, valid):
            self.assertIsNotNone(mol)
            for atom in mol.GetAtoms():
                atom.SetProp("instance_id", "")

        benzene_detection = _classify_vinylene_linkage_bonds(benzene).to_metadata()
        cyclohexene_detection = _classify_vinylene_linkage_bonds(cyclohexene).to_metadata()
        unactivated_detection = _classify_vinylene_linkage_bonds(unactivated).to_metadata()
        valid_detection = _classify_vinylene_linkage_bonds(valid).to_metadata()

        self.assertEqual(benzene_detection["small_ring_rejected_bond_count"], 3)
        self.assertEqual(benzene_detection["accepted_bond_count"], 0)
        self.assertEqual(cyclohexene_detection["small_ring_rejected_bond_count"], 1)
        self.assertEqual(cyclohexene_detection["accepted_bond_count"], 0)
        self.assertEqual(unactivated_detection["endpoint_rejected_bond_count"], 1)
        self.assertEqual(unactivated_detection["accepted_bond_count"], 0)
        self.assertEqual(valid_detection["accepted_bond_count"], 1)

    def test_recognized_boron_linkage_does_not_suppress_separate_vinylene_match(self):
        cases = (
            (
                "boest",
                "CC(=O)C=Cc1ccccc1.B(Oc1ccccc1)(Oc1ccccc1)c1ccccc1",
            ),
            (
                "boroxine",
                "CC(=O)C=Cc1ccccc1.c1ccc(B2OB(c3ccccc3)OB(c3ccccc3)O2)cc1",
            ),
        )
        for expected_boron_linkage, smiles in cases:
            with self.subTest(boron_linkage=expected_boron_linkage):
                mol = Chem.MolFromSmiles(smiles)
                self.assertIsNotNone(mol)
                for atom in mol.GetAtoms():
                    atom.SetProp("instance_id", "")

                detection = _classify_vinylene_linkage_bonds(mol).to_metadata()

                self.assertEqual(detection["recognized_boron_linkages"], [expected_boron_linkage])
                self.assertFalse(detection["boron_linkage_override_applied"])
                self.assertEqual(detection["boron_linkage_rejected_bond_count"], 0)
                self.assertEqual(detection["accepted_bond_count"], 1)

    def test_imine_and_vinylene_decompositions_make_triazine_assignment_ambiguous(self):
        cases = (
            (
                "imine_bridge",
                "imine",
                _record("melamine", MELAMINE, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
            ),
            (
                "vinylene_bridge",
                "vinylene",
                _record("tmt", TMT, "activated_methylene", 3),
                _record("pda", PDA, "aldehyde", 2),
            ),
        )
        for template_id, expected_linkage, first, second in cases:
            with self.subTest(linkage=expected_linkage), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

                preferred = decompose_cif_to_cofid(
                    input_cif,
                    topology="hcb",
                    linkage=expected_linkage,
                )
                triazine = decompose_cif_to_cofid(
                    input_cif,
                    topology="hcb",
                    linkage="triazine",
                )

                self.assertTrue(preferred.ok, preferred.reason)
                self.assertFalse(triazine.ok)
                self.assertIn("assignment is ambiguous", triazine.reason)
                resolution = triazine.metadata["triazine_linkage_resolution"]
                self.assertTrue(resolution["ambiguous"])
                self.assertFalse(resolution["override_applied"])
                self.assertEqual(
                    resolution["coexisting_linkages"],
                    [expected_linkage],
                )
                self.assertEqual(
                    resolution["evaluations"][expected_linkage]["status"],
                    "recognized",
                )
                self.assertGreater(triazine.metadata["n_triazine_rings"], 0)

    def test_vinylene_decomposition_survives_generic_labels_with_explicit_bonds(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator("vinylene_bridge").generate_pair_candidate(
                _record("tmt", TMT, "activated_methylene", 3),
                _record("pda", PDA, "aldehyde", 2),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")
            generic_labels = _with_generic_atom_labels(input_cif, temp_path / "generic_labels.cif")

            result = decompose_cif_to_cofid(
                generic_labels,
                topology="hcb",
                linkage="vinylene",
            )

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["bond_source"], "explicit_cif")
        self.assertEqual(result.metadata["vinylene_linkage_detection"]["accepted_bond_count"], 6)

    def test_decompose_infers_bond_orders_when_cif_bond_type_column_is_absent(self):
        for template_id, linkage, first, second, expected_linkage_bonds in ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cif_bond_type_column(summary.cif_path, temp_path / "without_bond_type.cif")

                result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage=linkage)

                self.assertTrue(result.ok, result.reason)
                self.assertEqual(result.cofid, summary.metadata["cofid"])
                self.assertEqual(result.metadata["bond_source"], "explicit_cif")
                self.assertGreater(result.metadata["n_missing_cif_bond_orders"], 0)
                self.assertGreater(result.metadata["n_bond_orders_inferred"], 0)
                self.assertEqual(result.metadata[f"n_{linkage}_linkage_bonds"], expected_linkage_bonds)

    def test_decompose_infers_bonds_when_cif_bond_loop_is_absent(self):
        for template_id, linkage, first, second, expected_linkage_bonds in ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cif_bond_loop(summary.cif_path, temp_path / "without_bond_loop.cif")

                result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage=linkage)

                self.assertTrue(result.ok, result.reason)
                self.assertEqual(result.cofid, summary.metadata["cofid"])
                self.assertEqual(result.metadata["bond_source"], "distance_inferred")
                self.assertGreater(result.metadata["n_distance_inferred_bonds"], 0)
                self.assertGreater(result.metadata["n_bond_orders_inferred"], 0)
                self.assertEqual(result.metadata[f"n_{linkage}_linkage_bonds"], expected_linkage_bonds)

    def test_decompose_auto_detects_topology_without_cofid_bonds_or_instance_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            no_bonds = _without_cif_bond_loop(summary.cif_path, temp_path / "without_bond_loop.cif")
            generic_labels = _with_generic_atom_labels(no_bonds, temp_path / "generic_labels.cif")

            result = decompose_cif_to_cofid(generic_labels)

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.topology, "hcb")
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["bond_source"], "distance_inferred")
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 3)
        detection = result.metadata["topology_detection"]
        self.assertEqual(detection["selected_topology"], "hcb")
        self.assertIn(detection["confidence"], {"exact", "high"})

    def test_generated_decompose_cases_cover_buildable_binary_linkages(self):
        library = ReactionLibrary.builtin()
        buildable = {
            template_id
            for template_id, profile in library.linkage_profiles.items()
            if profile.supports_binary_bridge_pair_generation and profile.supports_atomistic_realization
        }
        covered = {
            template_id
            for template_id, _linkage, _first, _second, _expected in ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES
        }

        self.assertEqual(covered, buildable)

    def test_decompose_accepts_template_id_linkage_alias(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator("hydrazone_bridge").generate_pair_candidate(
                _record("hdz", COF42_HYDRAZIDE, "hydrazide", 2),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage="hydrazone_bridge")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.linkage, "hydrazone")
        self.assertEqual(result.cofid, summary.metadata["cofid"])

    def test_random_topology_enabled_round_trips_cover_connectivity_combinations(self):
        selected_cases = random.Random(20260625).sample(
            TOPOLOGY_ENABLED_ROUND_TRIP_POOL,
            k=min(8, len(TOPOLOGY_ENABLED_ROUND_TRIP_POOL)),
        )
        observed_connectivity_pairs: set[tuple[int, int]] = set()
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            for case_name, template_id, topology, first, second, allow_repairable_geometry in selected_cases:
                with self.subTest(case=case_name):
                    config_kwargs = {
                        "allowed_reactions": (template_id,),
                        "topology_ids": (topology,),
                        "single_node_topology_ids": (topology,),
                        "enumerate_all_topologies": False,
                        "write_cif": True,
                        "rdkit_num_conformers": 1,
                        "retain_top_results": 1,
                    }
                    if allow_repairable_geometry:
                        config_kwargs.update(
                            {
                                "hard_hard_max_bridge_distance": 10.0,
                                "validation_thresholds": CoarseValidationThresholds(
                                    hard_hard_max_bridge_distance=10.0
                                ),
                            }
                        )
                    generator = BatchStructureGenerator(BatchGenerationConfig(**config_kwargs))
                    summary, _candidate = generator.generate_pair_candidate(
                        first,
                        second,
                        out_dir=temp_path / case_name,
                        write_cif=True,
                    )
                    self.assertEqual(summary.status, "ok")
                    self.assertEqual(summary.topology_id, topology)
                    self.assertIsNotNone(summary.cif_path)
                    input_cif = _without_cofid_comment(
                        summary.cif_path,
                        temp_path / case_name / "stripped.cif",
                    )

                    result = decompose_cif_to_cofid(input_cif, topology=topology, linkage=template_id)

                    if case_name == "sql_4_4_imine":
                        self.assertFalse(result.ok)
                        self.assertIn("periodic rank 3", result.reason)
                        self.assertEqual(
                            result.metadata["topology_validation"]["expected_periodic_rank"],
                            2,
                        )
                        self.assertEqual(
                            result.metadata["topology_validation"]["observed_periodic_rank"],
                            3,
                        )
                        continue
                    self.assertTrue(result.ok, result.reason)
                    self.assertEqual(result.cofid, summary.metadata["cofid"])
                    observed_connectivity_pairs.add(
                        tuple(sorted((first.expected_connectivity, second.expected_connectivity), reverse=True))
                    )

        self.assertGreaterEqual(
            observed_connectivity_pairs,
            {(3, 3), (3, 2), (4, 2), (4, 4), (6, 2)},
        )


class RingDecomposeRoundTripTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ring_cases = (
            (
                "boroxine_trimerization",
                "boroxine",
                build_rdkit_monomer(
                    "cof1_precursor",
                    "benzene-1,4-diboronic acid",
                    "OB(O)c1ccc(B(O)O)cc1",
                    "boronic_acid",
                    num_conformers=1,
                ),
            ),
            (
                "triazine_trimerization",
                "triazine",
                build_rdkit_monomer(
                    "ctf1_precursor",
                    "terephthalonitrile",
                    "N#Cc1ccc(C#N)cc1",
                    "nitrile",
                    num_conformers=1,
                ),
            ),
        )

    def _write_ring_candidate(
        self,
        root: Path,
        template_id: str,
        monomer,
        *,
        stacking_id: str | None = None,
    ) -> tuple[Path, str]:
        generator = RingFormingStructureGenerator(
            RingFormationConfig(stacking_ids=() if stacking_id is None else (stacking_id,))
        )
        candidate = generator.generate(monomer, template_id)
        cofid = generate_cofid(candidate, {monomer.id: monomer})
        path = root / f"{template_id}{'' if stacking_id is None else f'__{stacking_id}'}.cif"
        CIFWriter().write_candidate(
            path,
            candidate,
            {monomer.id: monomer},
            cofid=cofid,
            cofid_comment_suffix=(
                None if stacking_id is None else candidate.metadata["stacking"]["comment_suffix"]
            ),
        )
        return path, cofid

    def test_ring_decomposition_recovers_precursors_and_topology_without_comment(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            for template_id, linkage, monomer in self.ring_cases:
                with self.subTest(linkage=linkage):
                    cif_path, expected_cofid = self._write_ring_candidate(temp_path, template_id, monomer)
                    input_cif = _without_cofid_comment(cif_path, temp_path / f"{linkage}_stripped.cif")

                    result = decompose_cif_to_cofid(input_cif, linkage=linkage)

                    self.assertTrue(result.ok, result.reason)
                    self.assertEqual(result.cofid, expected_cofid)
                    self.assertEqual(result.topology, "hcb")
                    self.assertEqual(result.metadata["decomposition_strategy"], "ring")
                    self.assertEqual(result.metadata[f"n_{linkage}_rings"], 2)
                    self.assertEqual(result.metadata[f"n_{linkage}_ring_bonds"], 12)
                    self.assertEqual(result.metadata["n_topology_components"], 1)
                    self.assertEqual(result.metadata["topology_graph"]["node_connectivities"], [3, 3])
                    self.assertEqual(result.monomers[0].amount, 3)
                    if linkage == "triazine":
                        resolution = result.metadata["triazine_linkage_resolution"]
                        self.assertFalse(resolution["override_applied"])
                        self.assertEqual(resolution["recognized_higher_priority_linkages"], [])

    def test_event_mode_round_trips_ring_linkages_with_atomic_events(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            for template_id, linkage, monomer in self.ring_cases:
                with self.subTest(linkage=linkage):
                    cif_path, expected_cofid = self._write_ring_candidate(temp_path, template_id, monomer)
                    input_cif = _without_cofid_comment(
                        cif_path,
                        temp_path / f"{linkage}_event_stripped.cif",
                    )

                    result = decompose_cif_to_cofid(
                        input_cif,
                        linkage=linkage,
                        decomposition_mode="event",
                    )
                    automatic = decompose_cif_to_cofid(
                        input_cif,
                        decomposition_mode="event",
                    )

                    self.assertTrue(result.ok, result.reason)
                    self.assertEqual(result.cofid, expected_cofid)
                    self.assertTrue(automatic.ok, automatic.reason)
                    self.assertEqual(automatic.linkage, linkage)
                    self.assertEqual(automatic.cofid, expected_cofid)
                    self.assertEqual(result.metadata["event_status"], "SUCCESS_COMPLETE")
                    family_events = [
                        event
                        for event in result.metadata["event_detection"]["events"]
                        if event["family"] == linkage
                    ]
                    if linkage == "boroxine":
                        self.assertEqual(len(family_events), 2)
                        self.assertTrue(all(len(event["cut_bonds"]) == 6 for event in family_events))
                    else:
                        self.assertEqual(len(family_events), 4)
                        self.assertEqual(len({event["site_id"] for event in family_events}), 2)
                        self.assertTrue(all(len(event["cut_bonds"]) == 3 for event in family_events))
                        self.assertIn("N#C", result.monomers[0].canonical_smiles)

    def test_ring_decomposition_works_without_bond_loop_or_instance_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            for template_id, linkage, monomer in self.ring_cases:
                with self.subTest(linkage=linkage):
                    cif_path, expected_cofid = self._write_ring_candidate(temp_path, template_id, monomer)
                    no_bonds = _without_cif_bond_loop(cif_path, temp_path / f"{linkage}_no_bonds.cif")
                    generic = _with_generic_atom_labels(no_bonds, temp_path / f"{linkage}_generic.cif")

                    result = decompose_cif_to_cofid(generic, topology="hcb", linkage=template_id)

                    self.assertTrue(result.ok, result.reason)
                    self.assertEqual(result.cofid, expected_cofid)
                    self.assertEqual(result.linkage, linkage)
                    self.assertEqual(result.metadata["bond_source"], "distance_inferred")
                    self.assertGreater(result.metadata["n_bond_orders_inferred"], 0)

    def test_ring_decomposition_retains_higher_connected_precursor_nodes(self):
        node_cases = (
            (
                "boroxine_trimerization",
                "boroxine",
                "OB(O)c1cc(B(O)O)cc(B(O)O)c1",
                "boronic_acid",
            ),
            (
                "triazine_trimerization",
                "triazine",
                "N#Cc1cc(C#N)cc(C#N)c1",
                "nitrile",
            ),
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            for template_id, linkage, smiles, motif_kind in node_cases:
                with self.subTest(linkage=linkage):
                    monomer = build_rdkit_monomer(
                        f"{linkage}_node",
                        f"three-connected {linkage} precursor",
                        smiles,
                        motif_kind,
                        num_conformers=1,
                    )
                    cif_path, expected_cofid = self._write_ring_candidate(temp_path, template_id, monomer)
                    input_cif = _without_cofid_comment(cif_path, temp_path / f"{linkage}_node_stripped.cif")

                    result = decompose_cif_to_cofid(input_cif, linkage=linkage)

                    self.assertTrue(result.ok, result.reason)
                    self.assertEqual(result.cofid, expected_cofid)
                    self.assertEqual(result.monomers[0].connectivity, 3)
                    self.assertEqual(result.metadata["topology_graph"]["node_connectivities"], [3, 3])
                    self.assertEqual(result.metadata["topology_graph"]["n_edges"], 3)

    def test_ring_decomposition_normalizes_all_stacked_layer_registries(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            for template_id, linkage, monomer in self.ring_cases:
                for stacking_id in ("AA", "AB", "slipped"):
                    with self.subTest(linkage=linkage, stacking=stacking_id):
                        cif_path, expected_cofid = self._write_ring_candidate(
                            temp_path,
                            template_id,
                            monomer,
                            stacking_id=stacking_id,
                        )
                        input_cif = _without_cofid_comment(
                            cif_path,
                            temp_path / f"{linkage}_{stacking_id}_stripped.cif",
                        )

                        result = decompose_cif_to_cofid(input_cif, linkage=linkage)
                        validation = validate_cif_against_cofid(expected_cofid, cif_path)

                        self.assertTrue(result.ok, result.reason)
                        self.assertEqual(result.cofid, expected_cofid)
                        self.assertEqual(result.metadata[f"n_{linkage}_rings"], 4)
                        self.assertEqual(result.metadata["n_topology_components"], 2)
                        self.assertEqual(result.metadata["topology_graph"]["n_nodes"], 2)
                        self.assertEqual(result.metadata["topology_graph"]["n_edges"], 3)
                        self.assertEqual(result.monomers[0].amount, 6)
                        self.assertTrue(validation.ok, validation.reason)


if __name__ == "__main__":
    unittest.main()
