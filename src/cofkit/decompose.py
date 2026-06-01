"""CIF-to-COFid decomposition helpers.

The decomposition approach here is adapted from the
deCOFpose project: https://github.com/r-fedorov/deCOFpose
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from math import sqrt
from pathlib import Path
from typing import Callable, Mapping

from .bond_types import cif_type_to_bond_order, is_aromatic_bond_order
from .cofid import COFidMonomer, canonicalize_smiles, serialize_cofid
from .decompose_cif import PeriodicCifAtoms, read_periodic_cif_atoms

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - project dependency guard
    Chem = None


@dataclass(frozen=True)
class DecomposedMonomer:
    connectivity: int
    reactive_group: str
    canonical_smiles: str
    amount: int = 1

    def to_cofid_monomer(self) -> COFidMonomer:
        return COFidMonomer(
            connectivity=self.connectivity,
            reactive_group=self.reactive_group,
            canonical_smiles=self.canonical_smiles,
        )


@dataclass(frozen=True)
class CifDecompositionResult:
    status: str
    input_cif: str
    topology: str
    linkage: str
    cofid: str | None = None
    monomers: tuple[DecomposedMonomer, ...] = ()
    reason: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.status == "success" and self.cofid is not None

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "input_cif": self.input_cif,
            "topology": self.topology,
            "linkage": self.linkage,
            "cofid": self.cofid,
            "monomers": [
                {
                    "connectivity": monomer.connectivity,
                    "reactive_group": monomer.reactive_group,
                    "canonical_smiles": monomer.canonical_smiles,
                    "amount": monomer.amount,
                }
                for monomer in self.monomers
            ],
            "reason": self.reason,
            "metadata": dict(self.metadata),
        }


FragmentRepairer = Callable[[object], object]
LinkageMarker = Callable[[object], tuple[int, ...]]
ConnectivityCounter = Callable[[object, str], int]
Vec3 = tuple[float, float, float]
_CIF_NUMBER_RE = re.compile(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)(?:\(\d+\))?$")


@dataclass(frozen=True)
class LinkageDecompositionSpec:
    linkage_code: str
    template_id: str
    roles: tuple[str, ...]
    marker: LinkageMarker
    repairers: Mapping[str, FragmentRepairer] = field(default_factory=dict)
    connectivity_counter: ConnectivityCounter | None = None

    @property
    def metadata_key(self) -> str:
        return f"n_{self.linkage_code}_linkage_bonds"


@dataclass(frozen=True)
class BondCandidate:
    atom_idx_1: int
    atom_idx_2: int
    distance: float
    explicit_order: float | None = None


@dataclass(frozen=True)
class BondedMolBuildResult:
    mol: object
    metadata: Mapping[str, object] = field(default_factory=dict)


def decompose_cif_to_cofid(
    cif_path: str | Path,
    *,
    topology: str,
    linkage: str = "imine",
) -> CifDecompositionResult:
    input_path = Path(cif_path)
    spec = _resolve_linkage_spec(linkage)
    if spec is None:
        return CifDecompositionResult(
            status="skipped",
            input_cif=str(input_path),
            topology=topology,
            linkage=linkage,
            reason=(
                f"unsupported linkage {linkage!r}; current CIF decomposition supports "
                f"{tuple(sorted(_DECOMPOSITION_LINKAGE_ALIASES))!r}"
            ),
        )

    try:
        atoms = read_periodic_cif_atoms(input_path)
        monomers, metadata = _decompose_linkage_atoms(atoms, spec)
        if not monomers:
            return CifDecompositionResult(
                status="skipped",
                input_cif=str(input_path),
                topology=topology,
                linkage=spec.linkage_code,
                reason=f"no {spec.linkage_code} monomers were recovered",
                metadata=metadata,
            )
        cofid_monomers = tuple(
            sorted(
                (monomer.to_cofid_monomer() for monomer in monomers),
                key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
            )
        )
        cofid = serialize_cofid(monomers=cofid_monomers, topology=topology, linkage=spec.linkage_code)
        return CifDecompositionResult(
            status="success",
            input_cif=str(input_path),
            topology=topology,
            linkage=spec.linkage_code,
            cofid=cofid,
            monomers=monomers,
            metadata=metadata,
        )
    except Exception as exc:
        return CifDecompositionResult(
            status="skipped",
            input_cif=str(input_path),
            topology=topology,
            linkage=spec.linkage_code,
            reason=f"{type(exc).__name__}: {exc}",
        )


def _decompose_linkage_atoms(
    atoms: PeriodicCifAtoms,
    spec: LinkageDecompositionSpec,
) -> tuple[tuple[DecomposedMonomer, ...], dict[str, object]]:
    if Chem is None:
        raise RuntimeError("RDKit is required for CIF decomposition.")
    build_result = _build_bonded_mol(atoms)
    mol = build_result.mol
    linkage_bond_indices = spec.marker(mol)
    if not linkage_bond_indices:
        return (), {
            **dict(build_result.metadata),
            "n_atoms": mol.GetNumAtoms(),
            "n_bonds": mol.GetNumBonds(),
            spec.metadata_key: 0,
        }

    editable = Chem.RWMol(mol)
    for bond_idx in sorted(linkage_bond_indices, reverse=True):
        bond = editable.GetBondWithIdx(bond_idx)
        editable.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    fragments = Chem.GetMolFrags(editable.GetMol(), asMols=True, sanitizeFrags=False)
    monomers_by_key: dict[tuple[str, str, int], tuple[DecomposedMonomer, int]] = {}
    skipped_fragments = 0
    for fragment in fragments:
        monomer = _repair_linkage_fragment_to_monomer(fragment, spec)
        if monomer is None:
            skipped_fragments += 1
            continue
        key = (monomer.reactive_group, monomer.canonical_smiles, monomer.connectivity)
        previous = monomers_by_key.get(key)
        if previous is None:
            monomers_by_key[key] = (monomer, 1)
        else:
            monomers_by_key[key] = (previous[0], previous[1] + 1)

    monomers = tuple(
        sorted(
            (
                DecomposedMonomer(
                    connectivity=monomer.connectivity,
                    reactive_group=monomer.reactive_group,
                    canonical_smiles=monomer.canonical_smiles,
                    amount=amount,
                )
                for monomer, amount in monomers_by_key.values()
            ),
            key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
        )
    )
    return monomers, {
        **dict(build_result.metadata),
        "n_atoms": mol.GetNumAtoms(),
        "n_bonds": mol.GetNumBonds(),
        spec.metadata_key: len(linkage_bond_indices),
        "n_fragments_after_cut": len(fragments),
        "n_unique_monomers": len(monomers),
        "n_skipped_fragments": skipped_fragments,
    }


def _build_bonded_mol(atoms: PeriodicCifAtoms) -> BondedMolBuildResult:
    labels = tuple(atoms.info.get("_atom_site_label", ()))
    if len(labels) != len(atoms):
        raise ValueError("decomposition requires CIF atom labels aligned with atom sites")

    rw_mol = Chem.RWMol()
    for label, symbol in zip(labels, atoms.symbols):
        atom = Chem.Atom(symbol)
        atom.SetProp("cif_label", str(label))
        atom.SetProp("instance_id", _instance_id(str(label)))
        rw_mol.AddAtom(atom)

    candidates, metadata = _bond_candidates_from_cif_or_geometry(atoms, labels)
    if not candidates:
        raise ValueError("decomposition could not infer any covalent bonds from CIF atom positions")

    added_pairs: set[frozenset[int]] = set()
    for candidate in candidates:
        key = frozenset((candidate.atom_idx_1, candidate.atom_idx_2))
        if key in added_pairs:
            continue
        rw_mol.AddBond(candidate.atom_idx_1, candidate.atom_idx_2, Chem.BondType.SINGLE)
        added_pairs.add(key)

    mol = rw_mol.GetMol()
    ring_bonds = _ring_bond_keys(mol)
    inferred_order_count = 0
    for candidate in candidates:
        bond = mol.GetBondBetweenAtoms(candidate.atom_idx_1, candidate.atom_idx_2)
        if bond is None:
            continue
        order = candidate.explicit_order
        if order is None:
            order = _infer_bond_order(mol, candidate, ring_bonds)
            inferred_order_count += 1
        _apply_bond_order(mol, bond, order)

    mol.UpdatePropertyCache(strict=False)
    return BondedMolBuildResult(
        mol=mol,
        metadata={
            **metadata,
            "n_bond_orders_inferred": inferred_order_count,
        },
    )


def _bond_candidates_from_cif_or_geometry(
    atoms: PeriodicCifAtoms,
    labels: tuple[str, ...],
) -> tuple[tuple[BondCandidate, ...], dict[str, object]]:
    label_to_idx = {label: index for index, label in enumerate(labels)}
    bond_labels_1 = tuple(atoms.info.get("_geom_bond_atom_site_label_1", ()))
    bond_labels_2 = tuple(atoms.info.get("_geom_bond_atom_site_label_2", ()))
    if bond_labels_1 and bond_labels_2:
        candidates, n_missing_orders = _explicit_bond_candidates(atoms, label_to_idx, bond_labels_1, bond_labels_2)
        return candidates, {
            "bond_source": "explicit_cif",
            "n_explicit_cif_bonds": len(candidates),
            "n_missing_cif_bond_orders": n_missing_orders,
            "n_distance_inferred_bonds": 0,
        }

    candidates = _distance_inferred_bond_candidates(atoms)
    return candidates, {
        "bond_source": "distance_inferred",
        "n_explicit_cif_bonds": 0,
        "n_missing_cif_bond_orders": 0,
        "n_distance_inferred_bonds": len(candidates),
    }


def _explicit_bond_candidates(
    atoms: PeriodicCifAtoms,
    label_to_idx: Mapping[str, int],
    bond_labels_1: tuple[str, ...],
    bond_labels_2: tuple[str, ...],
) -> tuple[tuple[BondCandidate, ...], int]:
    bond_types = tuple(atoms.info.get("_ccdc_geom_bond_type") or atoms.info.get("_geom_bond_type") or ())
    bond_distances = tuple(atoms.info.get("_geom_bond_distance") or ())
    candidates: list[BondCandidate] = []
    added_pairs: set[frozenset[int]] = set()
    for row_idx, (label_1, label_2) in enumerate(zip(bond_labels_1, bond_labels_2)):
        if label_1 not in label_to_idx or label_2 not in label_to_idx:
            continue
        idx_1 = label_to_idx[label_1]
        idx_2 = label_to_idx[label_2]
        if idx_1 == idx_2:
            continue
        key = frozenset((idx_1, idx_2))
        if key in added_pairs:
            continue
        order = cif_type_to_bond_order(bond_types[row_idx] if row_idx < len(bond_types) else None)
        distance = _explicit_bond_distance(atoms, idx_1, idx_2, bond_distances, row_idx)
        candidates.append(BondCandidate(idx_1, idx_2, distance, explicit_order=order))
        added_pairs.add(key)
    return tuple(candidates), sum(1 for candidate in candidates if candidate.explicit_order is None)


def _explicit_bond_distance(
    atoms: PeriodicCifAtoms,
    idx_1: int,
    idx_2: int,
    bond_distances: tuple[str, ...],
    row_idx: int,
) -> float:
    if row_idx < len(bond_distances):
        try:
            return _parse_cif_float(bond_distances[row_idx])
        except ValueError:
            pass
    return _minimum_image_distance(atoms, idx_1, idx_2)


def _distance_inferred_bond_candidates(atoms: PeriodicCifAtoms) -> tuple[BondCandidate, ...]:
    candidates: list[BondCandidate] = []
    for idx_1, symbol_1 in enumerate(atoms.symbols):
        for idx_2 in range(idx_1 + 1, len(atoms.symbols)):
            symbol_2 = atoms.symbols[idx_2]
            distance = _minimum_image_distance(atoms, idx_1, idx_2)
            if _is_plausible_bond_distance(symbol_1, symbol_2, distance):
                candidates.append(BondCandidate(idx_1, idx_2, distance, explicit_order=None))
    return tuple(candidates)


def _ring_bond_keys(mol) -> set[frozenset[int]]:
    ring_bonds: set[frozenset[int]] = set()
    Chem.GetSymmSSSR(mol)
    for bond_ring in mol.GetRingInfo().BondRings():
        if len(bond_ring) not in {5, 6}:
            continue
        for bond_idx in bond_ring:
            bond = mol.GetBondWithIdx(int(bond_idx))
            ring_bonds.add(frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
    return ring_bonds


def _infer_bond_order(mol, candidate: BondCandidate, ring_bonds: set[frozenset[int]]) -> float:
    atom_1 = mol.GetAtomWithIdx(candidate.atom_idx_1)
    atom_2 = mol.GetAtomWithIdx(candidate.atom_idx_2)
    z1 = atom_1.GetAtomicNum()
    z2 = atom_2.GetAtomicNum()
    pair = frozenset((z1, z2))
    distance = candidate.distance
    ring_key = frozenset((candidate.atom_idx_1, candidate.atom_idx_2))

    same_instance = atom_1.GetProp("instance_id") == atom_2.GetProp("instance_id")
    if same_instance and ring_key in ring_bonds and pair in {frozenset((6, 6)), frozenset((6, 7)), frozenset((6, 16))}:
        if 1.20 <= distance <= 1.50:
            return 1.5

    if 1 in pair:
        return 1.0
    if 5 in pair:
        return 1.0
    if pair == frozenset((6, 8)):
        return 2.0 if distance <= 1.30 else 1.0
    if pair == frozenset((6, 7)):
        if distance <= 1.22:
            return 3.0
        if not same_instance and distance <= 1.48:
            return 2.0
        return 2.0 if distance <= 1.28 else 1.0
    if pair == frozenset((6, 6)):
        if distance <= 1.24:
            return 3.0
        if not same_instance and distance <= 1.46:
            return 2.0
        return 2.0 if distance <= 1.37 else 1.0
    return 1.0


def _apply_bond_order(mol, bond, order: float) -> None:
    bond.SetBondType(_rdkit_bond_type(order))
    if is_aromatic_bond_order(order):
        bond.SetIsAromatic(True)
        bond.GetBeginAtom().SetIsAromatic(True)
        bond.GetEndAtom().SetIsAromatic(True)


def _mark_carbon_hetero_double_linkage_bonds(
    mol,
    *,
    carbon_role: str,
    hetero_role: str,
    hetero_atomic_num: int,
) -> tuple[int, ...]:
    bond_indices: list[int] = []
    for bond in mol.GetBonds():
        if abs(float(bond.GetBondTypeAsDouble()) - 2.0) > 1.0e-6:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        atomic_nums = {first.GetAtomicNum(), second.GetAtomicNum()}
        if atomic_nums != {6, hetero_atomic_num}:
            continue
        if first.GetProp("instance_id") == second.GetProp("instance_id"):
            continue
        carbon = first if first.GetAtomicNum() == 6 else second
        hetero = first if first.GetAtomicNum() == hetero_atomic_num else second
        carbon.SetProp("cofkit_decompose_role", carbon_role)
        hetero.SetProp("cofkit_decompose_role", hetero_role)
        bond_indices.append(bond.GetIdx())
    return tuple(bond_indices)


def _mark_imine_linkage_bonds(mol) -> tuple[int, ...]:
    return _mark_carbon_hetero_double_linkage_bonds(
        mol,
        carbon_role="aldehyde",
        hetero_role="amine",
        hetero_atomic_num=7,
    )


def _mark_hydrazone_linkage_bonds(mol) -> tuple[int, ...]:
    return _mark_carbon_hetero_double_linkage_bonds(
        mol,
        carbon_role="aldehyde",
        hetero_role="hydrazide",
        hetero_atomic_num=7,
    )


def _mark_azine_linkage_bonds(mol) -> tuple[int, ...]:
    return _mark_carbon_hetero_double_linkage_bonds(
        mol,
        carbon_role="aldehyde",
        hetero_role="hydrazine",
        hetero_atomic_num=7,
    )


def _mark_keto_enamine_linkage_bonds(mol) -> tuple[int, ...]:
    return _mark_carbon_hetero_double_linkage_bonds(
        mol,
        carbon_role="keto_aldehyde",
        hetero_role="amine",
        hetero_atomic_num=7,
    )


def _mark_vinylene_linkage_bonds(mol) -> tuple[int, ...]:
    bond_indices: list[int] = []
    for bond in mol.GetBonds():
        if abs(float(bond.GetBondTypeAsDouble()) - 2.0) > 1.0e-6:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        if first.GetAtomicNum() != 6 or second.GetAtomicNum() != 6:
            continue
        if first.GetProp("instance_id") == second.GetProp("instance_id"):
            continue
        aldehyde, activated = _orient_vinylene_atoms(first, second)
        aldehyde.SetProp("cofkit_decompose_role", "aldehyde")
        activated.SetProp("cofkit_decompose_role", "activated_methylene")
        bond_indices.append(bond.GetIdx())
    return tuple(bond_indices)


def _mark_boronate_ester_linkage_bonds(mol) -> tuple[int, ...]:
    bond_indices: list[int] = []
    for bond in mol.GetBonds():
        if abs(float(bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        atomic_nums = {first.GetAtomicNum(), second.GetAtomicNum()}
        if atomic_nums != {5, 8}:
            continue
        if first.GetProp("instance_id") == second.GetProp("instance_id"):
            continue
        boron = first if first.GetAtomicNum() == 5 else second
        oxygen = first if first.GetAtomicNum() == 8 else second
        boron.SetProp("cofkit_decompose_role", "boronic_acid")
        oxygen.SetProp("cofkit_decompose_role", "catechol")
        bond_indices.append(bond.GetIdx())
    return tuple(bond_indices)


def _repair_linkage_fragment_to_monomer(
    fragment,
    spec: LinkageDecompositionSpec,
) -> DecomposedMonomer | None:
    endpoint_roles = {
        atom.GetProp("cofkit_decompose_role")
        for atom in fragment.GetAtoms()
        if atom.HasProp("cofkit_decompose_role")
    }
    if len(endpoint_roles) != 1:
        return None
    reactive_group = next(iter(endpoint_roles))
    if reactive_group not in spec.roles:
        return None
    repaired = spec.repairers.get(reactive_group, _identity_fragment_repairer)(fragment)
    return _finalize_repaired_fragment(
        repaired,
        reactive_group,
        connectivity_counter=spec.connectivity_counter,
    )


def _identity_fragment_repairer(fragment):
    return fragment


def _restore_aldehyde_oxygens(fragment):
    return _restore_double_oxygens(fragment, "aldehyde")


def _restore_keto_aldehyde_fragment(fragment):
    editable = Chem.RWMol(fragment)
    for bond in tuple(editable.GetBonds()):
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        if {first.GetAtomicNum(), second.GetAtomicNum()} != {6, 8}:
            continue
        carbon = first if first.GetAtomicNum() == 6 else second
        if carbon.GetIsAromatic():
            bond.SetBondType(Chem.BondType.SINGLE)
    return _restore_double_oxygens(editable.GetMol(), "keto_aldehyde")


def _restore_boronic_acid_oxygens(fragment):
    editable = Chem.RWMol(fragment)
    endpoint_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == "boronic_acid"
    ]
    for boron_idx in endpoint_indices:
        boron = editable.GetAtomWithIdx(boron_idx)
        oxygen_neighbor_count = sum(1 for neighbor in boron.GetNeighbors() if neighbor.GetAtomicNum() == 8)
        for _ in range(max(0, 2 - oxygen_neighbor_count)):
            oxygen_idx = editable.AddAtom(Chem.Atom("O"))
            editable.AddBond(boron_idx, oxygen_idx, Chem.BondType.SINGLE)
    return editable.GetMol()


def _restore_double_oxygens(fragment, role: str):
    editable = Chem.RWMol(fragment)
    endpoint_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == role
    ]
    for carbon_idx in endpoint_indices:
        carbon = editable.GetAtomWithIdx(carbon_idx)
        if _has_double_bonded_oxygen(carbon):
            continue
        oxygen_idx = editable.AddAtom(Chem.Atom("O"))
        editable.AddBond(carbon_idx, oxygen_idx, Chem.BondType.DOUBLE)
    return editable.GetMol()


def _finalize_repaired_fragment(
    fragment,
    reactive_group: str,
    *,
    connectivity_counter: ConnectivityCounter | None = None,
) -> DecomposedMonomer | None:
    endpoint_count = (
        connectivity_counter(fragment, reactive_group)
        if connectivity_counter is not None
        else _default_connectivity_count(fragment, reactive_group)
    )
    if endpoint_count <= 0:
        return None

    mol = Chem.Mol(fragment)
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNumRadicalElectrons(0)
        if atom.HasProp("cofkit_decompose_role"):
            atom.ClearProp("cofkit_decompose_role")
        if atom.HasProp("cif_label"):
            atom.ClearProp("cif_label")
        if atom.HasProp("instance_id"):
            atom.ClearProp("instance_id")
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    Chem.SanitizeMol(mol)
    return DecomposedMonomer(
        connectivity=endpoint_count,
        reactive_group=reactive_group,
        canonical_smiles=canonicalize_smiles(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)),
    )


def _default_connectivity_count(fragment, reactive_group: str) -> int:
    return sum(
        1
        for atom in fragment.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == reactive_group
    )


def _connectivity_count(fragment, reactive_group: str) -> int:
    count = _default_connectivity_count(fragment, reactive_group)
    if reactive_group == "catechol":
        if count % 2 != 0:
            raise ValueError("catechol decomposition found an odd number of marked catechol oxygens")
        return count // 2
    return count


def _has_double_bonded_oxygen(atom) -> bool:
    return any(
        bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        for bond in atom.GetBonds()
    )


def _minimum_image_distance(atoms: PeriodicCifAtoms, idx_1: int, idx_2: int) -> float:
    frac_1 = atoms.fractional_positions[idx_1]
    frac_2 = atoms.fractional_positions[idx_2]
    delta = tuple(frac_2[axis] - frac_1[axis] for axis in range(3))
    minimum_delta = tuple(value - round(value) for value in delta)
    cartesian = _fractional_delta_to_cartesian(atoms.cell_basis, minimum_delta)  # type: ignore[arg-type]
    return _norm(cartesian)


def _fractional_delta_to_cartesian(cell_basis: tuple[Vec3, Vec3, Vec3], delta: Vec3) -> Vec3:
    return (
        cell_basis[0][0] * delta[0] + cell_basis[1][0] * delta[1] + cell_basis[2][0] * delta[2],
        cell_basis[0][1] * delta[0] + cell_basis[1][1] * delta[1] + cell_basis[2][1] * delta[2],
        cell_basis[0][2] * delta[0] + cell_basis[1][2] * delta[1] + cell_basis[2][2] * delta[2],
    )


def _norm(vector: Vec3) -> float:
    return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])


def _is_plausible_bond_distance(symbol_1: str, symbol_2: str, distance: float) -> bool:
    if distance < 0.35:
        return False
    z1 = _atomic_number(symbol_1)
    z2 = _atomic_number(symbol_2)
    if z1 == 1 and z2 == 1:
        return False
    if {z1, z2} == {5, 8}:
        return distance <= 2.05
    radius_sum = _covalent_radius(z1) + _covalent_radius(z2)
    tolerance = 0.30 if 1 in {z1, z2} else 0.45
    return distance <= radius_sum + tolerance


def _atomic_number(symbol: str) -> int:
    periodic_table = Chem.GetPeriodicTable()
    atomic_number = int(periodic_table.GetAtomicNumber(str(symbol)))
    if atomic_number <= 0:
        raise ValueError(f"cannot infer atomic number for element symbol {symbol!r}")
    return atomic_number


def _covalent_radius(atomic_number: int) -> float:
    radius = float(Chem.GetPeriodicTable().GetRcovalent(atomic_number))
    if radius > 0:
        return radius
    fallback = {
        1: 0.31,
        5: 0.84,
        6: 0.76,
        7: 0.71,
        8: 0.66,
        9: 0.57,
        15: 1.07,
        16: 1.05,
        17: 1.02,
    }
    return fallback.get(atomic_number, 0.75)


def _parse_cif_float(value: object) -> float:
    text = str(value).strip().strip("'\"")
    match = _CIF_NUMBER_RE.match(text)
    if match is None:
        raise ValueError(f"invalid CIF numeric value {value!r}")
    return float(match.group(1))


def _rdkit_bond_type(order: float):
    if is_aromatic_bond_order(order):
        return Chem.BondType.AROMATIC
    if abs(order - 2.0) <= 1.0e-6:
        return Chem.BondType.DOUBLE
    if abs(order - 3.0) <= 1.0e-6:
        return Chem.BondType.TRIPLE
    return Chem.BondType.SINGLE


def _instance_id(label: str) -> str:
    if "_" not in label:
        return ""
    return label.split("_", 1)[0]


def _orient_vinylene_atoms(first, second):
    first_score = _vinylene_aldehyde_score(first)
    second_score = _vinylene_aldehyde_score(second)
    if first_score > second_score:
        return first, second
    if second_score > first_score:
        return second, first
    first_instance = first.GetProp("instance_id")
    second_instance = second.GetProp("instance_id")
    return (first, second) if first_instance > second_instance else (second, first)


def _vinylene_aldehyde_score(atom) -> int:
    score = 0
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            score += 4
        elif neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
            score += 2
        elif neighbor.GetAtomicNum() in {7, 8, 15, 16}:
            score -= 1
    return score


def _resolve_linkage_spec(linkage: str) -> LinkageDecompositionSpec | None:
    return _DECOMPOSITION_LINKAGE_ALIASES.get(str(linkage).strip().lower())


_IMINE_SPEC = LinkageDecompositionSpec(
    linkage_code="imine",
    template_id="imine_bridge",
    roles=("amine", "aldehyde"),
    marker=_mark_imine_linkage_bonds,
    repairers={"aldehyde": _restore_aldehyde_oxygens},
    connectivity_counter=_connectivity_count,
)

_HYDRAZONE_SPEC = LinkageDecompositionSpec(
    linkage_code="hydrazone",
    template_id="hydrazone_bridge",
    roles=("hydrazide", "aldehyde"),
    marker=_mark_hydrazone_linkage_bonds,
    repairers={"aldehyde": _restore_aldehyde_oxygens},
    connectivity_counter=_connectivity_count,
)

_AZINE_SPEC = LinkageDecompositionSpec(
    linkage_code="azine",
    template_id="azine_bridge",
    roles=("hydrazine", "aldehyde"),
    marker=_mark_azine_linkage_bonds,
    repairers={"aldehyde": _restore_aldehyde_oxygens},
    connectivity_counter=_connectivity_count,
)

_BORONATE_ESTER_SPEC = LinkageDecompositionSpec(
    linkage_code="boest",
    template_id="boronate_ester_bridge",
    roles=("boronic_acid", "catechol"),
    marker=_mark_boronate_ester_linkage_bonds,
    repairers={"boronic_acid": _restore_boronic_acid_oxygens},
    connectivity_counter=_connectivity_count,
)

_KETO_ENAMINE_SPEC = LinkageDecompositionSpec(
    linkage_code="bken",
    template_id="keto_enamine_bridge",
    roles=("amine", "keto_aldehyde"),
    marker=_mark_keto_enamine_linkage_bonds,
    repairers={"keto_aldehyde": _restore_keto_aldehyde_fragment},
    connectivity_counter=_connectivity_count,
)

_VINYLENE_SPEC = LinkageDecompositionSpec(
    linkage_code="vinylene",
    template_id="vinylene_bridge",
    roles=("activated_methylene", "aldehyde"),
    marker=_mark_vinylene_linkage_bonds,
    repairers={"aldehyde": _restore_aldehyde_oxygens},
    connectivity_counter=_connectivity_count,
)

_DECOMPOSITION_LINKAGE_ALIASES: Mapping[str, LinkageDecompositionSpec] = {
    "imine": _IMINE_SPEC,
    "imine_bridge": _IMINE_SPEC,
    "hydrazone": _HYDRAZONE_SPEC,
    "hydrazone_bridge": _HYDRAZONE_SPEC,
    "azine": _AZINE_SPEC,
    "azine_bridge": _AZINE_SPEC,
    "boest": _BORONATE_ESTER_SPEC,
    "boronate_ester": _BORONATE_ESTER_SPEC,
    "boronate_ester_bridge": _BORONATE_ESTER_SPEC,
    "bken": _KETO_ENAMINE_SPEC,
    "beta-ketoenamine": _KETO_ENAMINE_SPEC,
    "beta_ketoenamine": _KETO_ENAMINE_SPEC,
    "keto_enamine_bridge": _KETO_ENAMINE_SPEC,
    "vinylene": _VINYLENE_SPEC,
    "vinylene_bridge": _VINYLENE_SPEC,
}


__all__ = [
    "CifDecompositionResult",
    "DecomposedMonomer",
    "decompose_cif_to_cofid",
]
