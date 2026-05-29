from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

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


def decompose_cif_to_cofid(
    cif_path: str | Path,
    *,
    topology: str,
    linkage: str = "imine",
) -> CifDecompositionResult:
    input_path = Path(cif_path)
    if linkage != "imine":
        return CifDecompositionResult(
            status="skipped",
            input_cif=str(input_path),
            topology=topology,
            linkage=linkage,
            reason=f"unsupported linkage {linkage!r}; current CIF decomposition supports 'imine'",
        )

    try:
        atoms = read_periodic_cif_atoms(input_path)
        monomers, metadata = _decompose_explicit_imine_atoms(atoms)
        if not monomers:
            return CifDecompositionResult(
                status="skipped",
                input_cif=str(input_path),
                topology=topology,
                linkage=linkage,
                reason="no imine monomers were recovered",
                metadata=metadata,
            )
        cofid_monomers = tuple(
            sorted(
                (monomer.to_cofid_monomer() for monomer in monomers),
                key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
            )
        )
        cofid = serialize_cofid(monomers=cofid_monomers, topology=topology, linkage=linkage)
        return CifDecompositionResult(
            status="success",
            input_cif=str(input_path),
            topology=topology,
            linkage=linkage,
            cofid=cofid,
            monomers=monomers,
            metadata=metadata,
        )
    except Exception as exc:
        return CifDecompositionResult(
            status="skipped",
            input_cif=str(input_path),
            topology=topology,
            linkage=linkage,
            reason=f"{type(exc).__name__}: {exc}",
        )


def _decompose_explicit_imine_atoms(
    atoms: PeriodicCifAtoms,
) -> tuple[tuple[DecomposedMonomer, ...], dict[str, object]]:
    if Chem is None:
        raise RuntimeError("RDKit is required for CIF decomposition.")
    mol = _build_explicit_bond_mol(atoms)
    imine_bond_indices = _mark_imine_linkage_bonds(mol)
    if not imine_bond_indices:
        return (), {
            "n_atoms": mol.GetNumAtoms(),
            "n_bonds": mol.GetNumBonds(),
            "n_imine_linkage_bonds": 0,
        }

    editable = Chem.RWMol(mol)
    for bond_idx in sorted(imine_bond_indices, reverse=True):
        bond = editable.GetBondWithIdx(bond_idx)
        editable.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    fragments = Chem.GetMolFrags(editable.GetMol(), asMols=True, sanitizeFrags=False)
    monomers_by_key: dict[tuple[str, str, int], tuple[DecomposedMonomer, int]] = {}
    skipped_fragments = 0
    for fragment in fragments:
        monomer = _repair_imine_fragment_to_monomer(fragment)
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
        "n_atoms": mol.GetNumAtoms(),
        "n_bonds": mol.GetNumBonds(),
        "n_imine_linkage_bonds": len(imine_bond_indices),
        "n_fragments_after_cut": len(fragments),
        "n_unique_monomers": len(monomers),
        "n_skipped_fragments": skipped_fragments,
    }


def _build_explicit_bond_mol(atoms: PeriodicCifAtoms):
    labels = tuple(atoms.info.get("_atom_site_label", ()))
    if len(labels) != len(atoms):
        raise ValueError("explicit decomposition requires CIF atom labels aligned with atom sites")

    rw_mol = Chem.RWMol()
    for label, symbol in zip(labels, atoms.symbols):
        atom = Chem.Atom(symbol)
        atom.SetProp("cif_label", str(label))
        atom.SetProp("instance_id", _instance_id(str(label)))
        rw_mol.AddAtom(atom)

    label_to_idx = {label: index for index, label in enumerate(labels)}
    bond_labels_1 = tuple(atoms.info.get("_geom_bond_atom_site_label_1", ()))
    bond_labels_2 = tuple(atoms.info.get("_geom_bond_atom_site_label_2", ()))
    if not bond_labels_1 or not bond_labels_2:
        raise ValueError("explicit decomposition requires _geom_bond_atom_site_label_1/2 loops")
    bond_types = tuple(atoms.info.get("_ccdc_geom_bond_type") or atoms.info.get("_geom_bond_type") or ())
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
        order = cif_type_to_bond_order(bond_types[row_idx] if row_idx < len(bond_types) else None) or 1.0
        bond_type = _rdkit_bond_type(order)
        rw_mol.AddBond(idx_1, idx_2, bond_type)
        bond = rw_mol.GetBondBetweenAtoms(idx_1, idx_2)
        if is_aromatic_bond_order(order):
            bond.SetIsAromatic(True)
            rw_mol.GetAtomWithIdx(idx_1).SetIsAromatic(True)
            rw_mol.GetAtomWithIdx(idx_2).SetIsAromatic(True)
        added_pairs.add(key)

    mol = rw_mol.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol


def _mark_imine_linkage_bonds(mol) -> tuple[int, ...]:
    imine_bond_indices: list[int] = []
    for bond in mol.GetBonds():
        if abs(float(bond.GetBondTypeAsDouble()) - 2.0) > 1.0e-6:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        atomic_nums = {first.GetAtomicNum(), second.GetAtomicNum()}
        if atomic_nums != {6, 7}:
            continue
        if first.GetProp("instance_id") == second.GetProp("instance_id"):
            continue
        carbon = first if first.GetAtomicNum() == 6 else second
        nitrogen = first if first.GetAtomicNum() == 7 else second
        carbon.SetProp("cofkit_decompose_role", "aldehyde")
        nitrogen.SetProp("cofkit_decompose_role", "amine")
        imine_bond_indices.append(bond.GetIdx())
    return tuple(imine_bond_indices)


def _repair_imine_fragment_to_monomer(fragment) -> DecomposedMonomer | None:
    endpoint_roles = {
        atom.GetProp("cofkit_decompose_role")
        for atom in fragment.GetAtoms()
        if atom.HasProp("cofkit_decompose_role")
    }
    if endpoint_roles == {"amine"}:
        return _finalize_repaired_fragment(fragment, "amine")
    if endpoint_roles == {"aldehyde"}:
        return _finalize_repaired_fragment(_restore_aldehyde_oxygens(fragment), "aldehyde")
    return None


def _restore_aldehyde_oxygens(fragment):
    editable = Chem.RWMol(fragment)
    endpoint_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == "aldehyde"
    ]
    for carbon_idx in endpoint_indices:
        carbon = editable.GetAtomWithIdx(carbon_idx)
        if _has_double_bonded_oxygen(carbon):
            continue
        oxygen_idx = editable.AddAtom(Chem.Atom("O"))
        editable.AddBond(carbon_idx, oxygen_idx, Chem.BondType.DOUBLE)
    return editable.GetMol()


def _finalize_repaired_fragment(fragment, reactive_group: str) -> DecomposedMonomer | None:
    endpoint_count = sum(
        1
        for atom in fragment.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == reactive_group
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


def _has_double_bonded_oxygen(atom) -> bool:
    return any(
        bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        for bond in atom.GetBonds()
    )


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


__all__ = [
    "CifDecompositionResult",
    "DecomposedMonomer",
    "decompose_cif_to_cofid",
]
