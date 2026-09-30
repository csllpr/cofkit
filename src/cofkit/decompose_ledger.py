"""Decomposition atom ledger and charge/stereo scope reporting.

This module owns the atom-accounting schema produced by the event
decomposition engine (:mod:`cofkit.decompose_events`).  Every atom of the
input CIF is classified into exactly one ledger bucket:

- ``recovered_fragments`` — product atoms mapped to a fragment that was
  repaired into a recovered precursor;
- ``guests`` — disconnected components with no linkage-event atoms
  (one ledger entry per guest *molecule*, each with its own atom count);
- ``residue`` — product-side atoms excluded from precursor recovery by the
  linkage chemistry (currently boroxine ring oxygens);
- ``unaccounted`` — framework fragments that could not be assigned to a
  precursor role (an explicit failure signal, never silently dropped).

Precursor repair changes atom counts relative to the input (condensation
byproduct atoms are restored, product-state hydrogens are re-expressed as
implicit hydrogens), so the ledger records those changes as
``reaction_additions`` / ``reaction_deletions`` and checks the explicit-atom
balance

    input + reaction_additions
        = recovered_precursor_atoms(explicit) + reaction_deletions
          + residue + guests + unaccounted

per element.  Recovered-precursor implicit hydrogens are reported as an
informational total so the full precursor formula content remains visible.

Charge/stereo scope is deliberately conservative: formal charges that survive
into recovered precursor SMILES are verified per fragment and reported with a
typed status; stereochemistry is never preserved by decomposition because
COFid canonicalization uses ``isomericSmiles=False``, so detected stereo
evidence yields an explicit ``unsupported`` marker rather than silent loss.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Mapping, Sequence

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - project dependency guard
    Chem = None


# Formal-charge scope statuses, recorded under
# ``metadata["atom_ledger"]["charge_scope"]["status"]``:
#   NOT_PRESENT — no formal charges on input atoms or recovered precursors;
#   PRESERVED — every charged input atom maps to a recovered fragment whose
#       precursor SMILES carries the same fragment net charge;
#   GUEST_LOCALIZED — formal charges exist only on guest molecules; they are
#       reported per guest but guests are not part of the serialized COFid;
#   RESTORED_BY_HEURISTIC — the input carried no formal charge but precursor
#       repair restored one from local valence (the tetravalent iminium/
#       pyridinium nitrogen rule in decompose._finalize_repaired_fragment);
#   LIMITED — charge information was present but could not be fully accounted
#       (charge altered or lost on a recovered fragment, or charge sitting on
#       residue/unaccounted atoms); an explicit unsupported/limited marker,
#       never silent loss.
CHARGE_SCOPE_NOT_PRESENT = "not_present"
CHARGE_SCOPE_PRESERVED = "preserved"
CHARGE_SCOPE_GUEST_LOCALIZED = "guest_localized"
CHARGE_SCOPE_RESTORED_BY_HEURISTIC = "restored_by_heuristic"
CHARGE_SCOPE_LIMITED = "limited"

# Stereochemistry scope statuses, recorded under
# ``metadata["atom_ledger"]["stereo_scope"]["status"]``.  Decomposition never
# preserves stereochemistry (canonical SMILES are non-isomeric), so the only
# question is whether the input showed stereo evidence that is now explicitly
# reported as unrecovered.
STEREO_SCOPE_NOT_DETECTED = "not_detected"
STEREO_SCOPE_UNSUPPORTED = "unsupported"
STEREO_SCOPE_UNDETERMINED = "undetermined"

_STEREO_SERIALIZATION_POLICY = (
    "COFid and decomposition canonical SMILES are serialized with "
    "isomericSmiles=False; tetrahedral chirality and E/Z bond geometry are "
    "not preserved or certified by decomposition"
)

_RESIDUE_POLICY_BY_FAMILY = {
    "boroxine": (
        "boroxine ring oxygens are product-side residue excluded from "
        "precursor recovery; the boronic-acid hydroxyl oxygens are re-added "
        "during fragment repair (net 3 H2O eliminated per ring)"
    ),
}

_BALANCE_EQUATION = (
    "input_atoms + reaction_additions = recovered_precursor_explicit_atoms "
    "+ reaction_deletions + residue_atoms + guest_atoms + unaccounted_atoms"
)


@dataclass(frozen=True)
class LedgerGuest:
    """One disconnected non-framework component (a guest molecule)."""

    fragment_id: int
    atom_indices: tuple[int, ...]
    atom_count: int
    elements: Mapping[str, int]
    molecular_formula: str
    net_formal_charge: int

    def to_dict(self) -> dict[str, object]:
        return {
            "fragment_id": self.fragment_id,
            "atom_indices": list(self.atom_indices),
            "atom_count": self.atom_count,
            "elements": dict(self.elements),
            "molecular_formula": self.molecular_formula,
            "net_formal_charge": self.net_formal_charge,
        }


@dataclass(frozen=True)
class LedgerRecoveredFragment:
    """Product atoms mapped to one fragment repaired into a precursor."""

    fragment_id: int
    role: str
    atom_indices: tuple[int, ...]
    atom_count: int
    elements: Mapping[str, int]
    net_formal_charge: int
    recovered_canonical_smiles: str
    recovered_elements: Mapping[str, int] | None
    recovered_implicit_hydrogen_count: int | None
    recovered_net_formal_charge: int | None
    added_elements: Mapping[str, int]
    removed_elements: Mapping[str, int]

    @property
    def recovered_explicit_atom_count(self) -> int | None:
        if self.recovered_elements is None:
            return None
        return sum(self.recovered_elements.values())

    def to_dict(self) -> dict[str, object]:
        return {
            "fragment_id": self.fragment_id,
            "role": self.role,
            "atom_indices": list(self.atom_indices),
            "atom_count": self.atom_count,
            "elements": dict(self.elements),
            "net_formal_charge": self.net_formal_charge,
            "recovered_canonical_smiles": self.recovered_canonical_smiles,
            "recovered_elements": (
                None
                if self.recovered_elements is None
                else dict(self.recovered_elements)
            ),
            "recovered_explicit_atom_count": self.recovered_explicit_atom_count,
            "recovered_implicit_hydrogen_count": self.recovered_implicit_hydrogen_count,
            "recovered_net_formal_charge": self.recovered_net_formal_charge,
            "added_elements": dict(self.added_elements),
            "removed_elements": dict(self.removed_elements),
        }


@dataclass(frozen=True)
class LedgerUnaccountedFragment:
    """A framework fragment not assignable to any recovered precursor role."""

    fragment_id: int
    atom_indices: tuple[int, ...]
    atom_count: int
    elements: Mapping[str, int]
    net_formal_charge: int
    reason: str

    def to_dict(self) -> dict[str, object]:
        return {
            "fragment_id": self.fragment_id,
            "atom_indices": list(self.atom_indices),
            "atom_count": self.atom_count,
            "elements": dict(self.elements),
            "net_formal_charge": self.net_formal_charge,
            "reason": self.reason,
        }


@dataclass(frozen=True)
class LedgerBalance:
    """Per-element explicit-atom balance over the ledger buckets."""

    equation: str
    input_atom_count: int
    reaction_additions: Mapping[str, int]
    recovered_precursor_explicit_atom_count: int | None
    reaction_deletions: Mapping[str, int]
    residue_atom_count: int
    guest_atom_count: int
    unaccounted_atom_count: int
    per_element_residuals: Mapping[str, int] | None
    balanced: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "equation": self.equation,
            "input_atom_count": self.input_atom_count,
            "reaction_additions": dict(self.reaction_additions),
            "recovered_precursor_explicit_atom_count": (
                self.recovered_precursor_explicit_atom_count
            ),
            "reaction_deletions": dict(self.reaction_deletions),
            "residue_atom_count": self.residue_atom_count,
            "guest_atom_count": self.guest_atom_count,
            "unaccounted_atom_count": self.unaccounted_atom_count,
            "per_element_residuals": (
                None
                if self.per_element_residuals is None
                else dict(self.per_element_residuals)
            ),
            "balanced": self.balanced,
        }


@dataclass(frozen=True)
class ChargeScopeReport:
    """How formal charges in the input were accounted for."""

    status: str
    input_net_formal_charge: int
    input_charged_atom_count: int
    formal_charge_column_present: bool | None
    recovered_precursor_net_formal_charge: int
    guest_net_formal_charge: int
    residue_net_formal_charge: int
    unaccounted_net_formal_charge: int
    heuristic_restored_charge: int
    net_charge_delta: int
    charge_conserved: bool
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "input_net_formal_charge": self.input_net_formal_charge,
            "input_charged_atom_count": self.input_charged_atom_count,
            "formal_charge_column_present": self.formal_charge_column_present,
            "recovered_precursor_net_formal_charge": (
                self.recovered_precursor_net_formal_charge
            ),
            "guest_net_formal_charge": self.guest_net_formal_charge,
            "residue_net_formal_charge": self.residue_net_formal_charge,
            "unaccounted_net_formal_charge": self.unaccounted_net_formal_charge,
            "heuristic_restored_charge": self.heuristic_restored_charge,
            "net_charge_delta": self.net_charge_delta,
            "charge_conserved": self.charge_conserved,
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class StereoScopeReport:
    """Stereochemistry scope: decomposition never preserves stereo."""

    status: str
    preserved: bool
    perception_status: str
    chiral_atom_count: int
    stereo_bond_count: int
    policy: str = _STEREO_SERIALIZATION_POLICY

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "preserved": self.preserved,
            "perception_status": self.perception_status,
            "chiral_atom_count": self.chiral_atom_count,
            "stereo_bond_count": self.stereo_bond_count,
            "policy": self.policy,
        }


@dataclass(frozen=True)
class AtomLedger:
    """Full atom accounting of one event-mode cut/reconstruction."""

    family: str
    input_atom_count: int
    input_elements: Mapping[str, int]
    recovered_fragments: tuple[LedgerRecoveredFragment, ...]
    guests: tuple[LedgerGuest, ...]
    residue_atom_indices: tuple[int, ...]
    residue_elements: Mapping[str, int]
    residue_policy: str | None
    unaccounted: tuple[LedgerUnaccountedFragment, ...]
    reaction_additions: Mapping[str, int]
    reaction_deletions: Mapping[str, int]
    recovered_precursor_implicit_hydrogen_count: int | None
    recovered_precursor_total_atom_count: int | None
    balance: LedgerBalance
    charge_scope: ChargeScopeReport
    stereo_scope: StereoScopeReport

    @property
    def guest_molecule_count(self) -> int:
        return len(self.guests)

    @property
    def guest_atom_count(self) -> int:
        return sum(guest.atom_count for guest in self.guests)

    @property
    def residue_atom_count(self) -> int:
        return len(self.residue_atom_indices)

    @property
    def unaccounted_atom_count(self) -> int:
        return sum(record.atom_count for record in self.unaccounted)

    def to_dict(self) -> dict[str, object]:
        return {
            "family": self.family,
            "input_atom_count": self.input_atom_count,
            "input_elements": dict(self.input_elements),
            "recovered_fragments": [
                fragment.to_dict() for fragment in self.recovered_fragments
            ],
            "guests": [guest.to_dict() for guest in self.guests],
            "guest_molecule_count": self.guest_molecule_count,
            "guest_atom_count": self.guest_atom_count,
            "residue_atom_count": self.residue_atom_count,
            "residue_atom_indices": list(self.residue_atom_indices),
            "residue_elements": dict(self.residue_elements),
            "residue_policy": self.residue_policy,
            "unaccounted": [record.to_dict() for record in self.unaccounted],
            "unaccounted_atom_count": self.unaccounted_atom_count,
            "reaction_additions": dict(self.reaction_additions),
            "reaction_deletions": dict(self.reaction_deletions),
            "recovered_precursor_implicit_hydrogen_count": (
                self.recovered_precursor_implicit_hydrogen_count
            ),
            "recovered_precursor_total_atom_count": (
                self.recovered_precursor_total_atom_count
            ),
            "balance": self.balance.to_dict(),
            "charge_scope": self.charge_scope.to_dict(),
            "stereo_scope": self.stereo_scope.to_dict(),
        }

    def summary_dict(self) -> dict[str, object]:
        """Compact per-hypothesis record (the full ledger rides on the result)."""
        return {
            "family": self.family,
            "input_atom_count": self.input_atom_count,
            "balanced": self.balance.balanced,
            "guest_molecule_count": self.guest_molecule_count,
            "guest_atom_count": self.guest_atom_count,
            "residue_atom_count": self.residue_atom_count,
            "unaccounted_atom_count": self.unaccounted_atom_count,
            "reaction_additions": dict(self.reaction_additions),
            "reaction_deletions": dict(self.reaction_deletions),
            "charge_scope_status": self.charge_scope.status,
            "stereo_scope_status": self.stereo_scope.status,
        }


def element_counts_for_atoms(mol, atom_indices: Sequence[int]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for atom_idx in atom_indices:
        counts[mol.GetAtomWithIdx(int(atom_idx)).GetSymbol()] += 1
    return dict(sorted(counts.items()))


def net_formal_charge_for_atoms(mol, atom_indices: Sequence[int]) -> int:
    return sum(
        int(mol.GetAtomWithIdx(int(atom_idx)).GetFormalCharge())
        for atom_idx in atom_indices
    )


def hill_formula(elements: Mapping[str, int]) -> str:
    symbols = [symbol for symbol, count in elements.items() if count > 0]
    if "C" in symbols:
        ordered = (
            ["C"]
            + (["H"] if "H" in symbols else [])
            + sorted(symbol for symbol in symbols if symbol not in {"C", "H"})
        )
    else:
        ordered = sorted(symbols)
    return "".join(
        f"{symbol}{'' if elements[symbol] == 1 else elements[symbol]}"
        for symbol in ordered
    )


def smiles_atom_accounting(smiles: str) -> tuple[dict[str, int], int, int] | None:
    """Explicit element counts, implicit-H count, and net charge of a SMILES."""
    if Chem is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    explicit: Counter[str] = Counter(atom.GetSymbol() for atom in mol.GetAtoms())
    implicit_hydrogens = sum(int(atom.GetTotalNumHs()) for atom in mol.GetAtoms())
    net_charge = sum(int(atom.GetFormalCharge()) for atom in mol.GetAtoms())
    return dict(sorted(explicit.items())), implicit_hydrogens, net_charge


def build_recovered_fragment_record(
    *,
    fragment_id: int,
    role: str,
    mol,
    atom_indices: Sequence[int],
    recovered_canonical_smiles: str,
) -> LedgerRecoveredFragment:
    """Account one role-bearing fragment against its repaired precursor.

    ``added_elements`` / ``removed_elements`` are the per-element explicit-atom
    deltas of precursor repair: additions are atoms restored from condensation
    byproducts (e.g. one O per imine aldehyde endpoint), removals are
    product-state atoms dropped during normalization (explicit hydrogens
    re-expressed as implicit hydrogens).
    """
    indices = tuple(int(atom_idx) for atom_idx in atom_indices)
    pre_elements = element_counts_for_atoms(mol, indices)
    pre_charge = net_formal_charge_for_atoms(mol, indices)
    accounting = smiles_atom_accounting(recovered_canonical_smiles)
    if accounting is None:
        return LedgerRecoveredFragment(
            fragment_id=fragment_id,
            role=role,
            atom_indices=indices,
            atom_count=len(indices),
            elements=pre_elements,
            net_formal_charge=pre_charge,
            recovered_canonical_smiles=recovered_canonical_smiles,
            recovered_elements=None,
            recovered_implicit_hydrogen_count=None,
            recovered_net_formal_charge=None,
            added_elements={},
            removed_elements={},
        )
    post_elements, implicit_hydrogens, post_charge = accounting
    added = {
        element: post_elements.get(element, 0) - pre_elements.get(element, 0)
        for element in set(post_elements) | set(pre_elements)
        if post_elements.get(element, 0) > pre_elements.get(element, 0)
    }
    removed = {
        element: pre_elements.get(element, 0) - post_elements.get(element, 0)
        for element in set(post_elements) | set(pre_elements)
        if pre_elements.get(element, 0) > post_elements.get(element, 0)
    }
    return LedgerRecoveredFragment(
        fragment_id=fragment_id,
        role=role,
        atom_indices=indices,
        atom_count=len(indices),
        elements=pre_elements,
        net_formal_charge=pre_charge,
        recovered_canonical_smiles=recovered_canonical_smiles,
        recovered_elements=post_elements,
        recovered_implicit_hydrogen_count=implicit_hydrogens,
        recovered_net_formal_charge=post_charge,
        added_elements=added,
        removed_elements=removed,
    )


def perceive_stereochemistry_evidence(
    mol,
    cartesian_positions: Sequence[Sequence[float]],
) -> dict[str, object]:
    """Detect input stereochemical evidence from 3D coordinates.

    Only tetrahedral atom chirality is perceived (this RDKit line has no
    bond-stereo-from-3D perception); E/Z evidence is therefore not detectable
    here and is covered by the standing ``unsupported`` policy instead.
    """
    unavailable = {
        "perception_status": "unavailable",
        "chiral_atom_count": 0,
        "stereo_bond_count": 0,
        "bond_stereo_perception": "unsupported_by_rdkit",
    }
    if Chem is None:
        return {**unavailable, "reason": "RDKit is not available"}
    try:
        probe = Chem.Mol(mol)
        if probe.GetNumAtoms() != len(cartesian_positions):
            raise ValueError("atom count does not match coordinate count")
        conformer = Chem.Conformer(probe.GetNumAtoms())
        for idx, position in enumerate(cartesian_positions):
            conformer.SetAtomPosition(idx, tuple(float(value) for value in position))
        probe.AddConformer(conformer, assignId=True)
        Chem.AssignStereochemistryFrom3D(probe, replaceExistingTags=True)
        chiral_atom_count = sum(
            1
            for atom in probe.GetAtoms()
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        )
        return {
            "perception_status": "perceived",
            "chiral_atom_count": chiral_atom_count,
            "stereo_bond_count": 0,
            "bond_stereo_perception": "unsupported_by_rdkit",
        }
    except Exception as exc:
        return {**unavailable, "reason": f"{type(exc).__name__}: {exc}"}


def _classify_charge_scope(
    *,
    mol,
    recovered_fragments: tuple[LedgerRecoveredFragment, ...],
    guests: tuple[LedgerGuest, ...],
    residue_atom_indices: tuple[int, ...],
    unaccounted: tuple[LedgerUnaccountedFragment, ...],
    formal_charge_column_present: bool | None,
) -> ChargeScopeReport:
    all_indices = range(mol.GetNumAtoms())
    input_net = net_formal_charge_for_atoms(mol, all_indices)
    input_charged = sum(
        1
        for atom_idx in all_indices
        if int(mol.GetAtomWithIdx(atom_idx).GetFormalCharge()) != 0
    )
    recovered_net = sum(
        fragment.recovered_net_formal_charge or 0 for fragment in recovered_fragments
    )
    guest_net = sum(guest.net_formal_charge for guest in guests)
    residue_net = net_formal_charge_for_atoms(mol, residue_atom_indices)
    unaccounted_net = sum(record.net_formal_charge for record in unaccounted)

    altered_fragments = []
    heuristic_restored_charge = 0
    preserved_fragment_charge = 0
    unverifiable_fragments = 0
    for fragment in recovered_fragments:
        if fragment.recovered_net_formal_charge is None:
            unverifiable_fragments += 1
            continue
        delta = fragment.recovered_net_formal_charge - fragment.net_formal_charge
        if delta == 0:
            preserved_fragment_charge += abs(fragment.net_formal_charge)
        elif fragment.net_formal_charge == 0 and delta != 0:
            # Charge created during repair by the valence-driven tetravalent-N
            # restoration rule, not carried from the input.
            heuristic_restored_charge += delta
        else:
            altered_fragments.append(fragment.fragment_id)

    net_charge_delta = (
        input_net + heuristic_restored_charge
        - (recovered_net + guest_net + residue_net + unaccounted_net)
    )
    charge_conserved = net_charge_delta == 0 and not altered_fragments

    notes: list[str] = []
    if formal_charge_column_present is False:
        notes.append(
            "the input CIF carried no _atom_site_pdbx_formal_charge column; "
            "input charges default to zero"
        )
    if guest_net != 0:
        notes.append(
            "guest molecules carry formal charge; guests are reported in this "
            "ledger but are not part of the serialized COFid"
        )
    if residue_net != 0:
        notes.append("formal charge sits on linkage residue atoms excluded from recovery")
    if unaccounted_net != 0:
        notes.append("formal charge sits on unaccounted framework atoms")

    if unverifiable_fragments or altered_fragments or residue_net != 0 or unaccounted_net != 0:
        status = CHARGE_SCOPE_LIMITED
        if altered_fragments:
            notes.append(
                "formal charge was altered or lost on recovered fragments "
                f"{altered_fragments!r}; recovery of the charged identity is unsupported"
            )
    elif preserved_fragment_charge > 0:
        status = CHARGE_SCOPE_PRESERVED
        if heuristic_restored_charge:
            notes.append(
                "additional charge was restored by the valence-driven "
                "tetravalent-nitrogen heuristic during repair"
            )
    elif heuristic_restored_charge:
        status = CHARGE_SCOPE_RESTORED_BY_HEURISTIC
        notes.append(
            "the input carried no formal charge; recovered precursor charge "
            "was restored from local valence (tetravalent iminium/pyridinium "
            "nitrogen heuristic), not from input data"
        )
    elif guest_net != 0:
        status = CHARGE_SCOPE_GUEST_LOCALIZED
    else:
        status = CHARGE_SCOPE_NOT_PRESENT

    return ChargeScopeReport(
        status=status,
        input_net_formal_charge=input_net,
        input_charged_atom_count=input_charged,
        formal_charge_column_present=formal_charge_column_present,
        recovered_precursor_net_formal_charge=recovered_net,
        guest_net_formal_charge=guest_net,
        residue_net_formal_charge=residue_net,
        unaccounted_net_formal_charge=unaccounted_net,
        heuristic_restored_charge=heuristic_restored_charge,
        net_charge_delta=net_charge_delta,
        charge_conserved=charge_conserved,
        notes=tuple(notes),
    )


def _stereo_scope_from_evidence(
    identity_evidence: Mapping[str, object] | None,
) -> StereoScopeReport:
    evidence = identity_evidence.get("stereochemistry") if identity_evidence else None
    if not isinstance(evidence, Mapping):
        return StereoScopeReport(
            status=STEREO_SCOPE_UNDETERMINED,
            preserved=False,
            perception_status="unavailable",
            chiral_atom_count=0,
            stereo_bond_count=0,
        )
    perception_status = str(evidence.get("perception_status", "unavailable"))
    chiral_atom_count = int(evidence.get("chiral_atom_count", 0))
    stereo_bond_count = int(evidence.get("stereo_bond_count", 0))
    if perception_status != "perceived":
        status = STEREO_SCOPE_UNDETERMINED
    elif chiral_atom_count or stereo_bond_count:
        status = STEREO_SCOPE_UNSUPPORTED
    else:
        status = STEREO_SCOPE_NOT_DETECTED
    return StereoScopeReport(
        status=status,
        preserved=False,
        perception_status=perception_status,
        chiral_atom_count=chiral_atom_count,
        stereo_bond_count=stereo_bond_count,
    )


def build_atom_ledger(
    *,
    family: str,
    mol,
    recovered_fragments: tuple[LedgerRecoveredFragment, ...],
    guests: tuple[LedgerGuest, ...],
    residue_atom_indices: tuple[int, ...],
    unaccounted: tuple[LedgerUnaccountedFragment, ...],
    identity_evidence: Mapping[str, object] | None = None,
) -> AtomLedger:
    """Assemble the ledger and check the explicit-atom balance per element."""
    input_atom_count = mol.GetNumAtoms()
    input_elements = element_counts_for_atoms(mol, range(input_atom_count))
    residue_elements = element_counts_for_atoms(mol, residue_atom_indices)

    additions: Counter[str] = Counter()
    deletions: Counter[str] = Counter()
    for fragment in recovered_fragments:
        additions.update(fragment.added_elements)
        deletions.update(fragment.removed_elements)

    guest_elements: Counter[str] = Counter()
    for guest in guests:
        guest_elements.update(guest.elements)
    unaccounted_elements: Counter[str] = Counter()
    for record in unaccounted:
        unaccounted_elements.update(record.elements)

    accounting_complete = all(
        fragment.recovered_elements is not None for fragment in recovered_fragments
    )
    if accounting_complete:
        precursor_explicit: Counter[str] = Counter()
        implicit_hydrogens = 0
        for fragment in recovered_fragments:
            assert fragment.recovered_elements is not None
            precursor_explicit.update(fragment.recovered_elements)
            implicit_hydrogens += int(fragment.recovered_implicit_hydrogen_count or 0)
        precursor_explicit_count: int | None = sum(precursor_explicit.values())
        per_element_residuals: dict[str, int] = {}
        all_elements = (
            set(input_elements)
            | set(additions)
            | set(deletions)
            | set(precursor_explicit)
            | set(residue_elements)
            | set(guest_elements)
            | set(unaccounted_elements)
        )
        for element in sorted(all_elements):
            left = input_elements.get(element, 0) + additions.get(element, 0)
            right = (
                precursor_explicit.get(element, 0)
                + deletions.get(element, 0)
                + residue_elements.get(element, 0)
                + guest_elements.get(element, 0)
                + unaccounted_elements.get(element, 0)
            )
            per_element_residuals[element] = left - right
        balanced = all(residual == 0 for residual in per_element_residuals.values())
        total_atom_count: int | None = precursor_explicit_count + implicit_hydrogens
        implicit_count: int | None = implicit_hydrogens
        residuals: Mapping[str, int] | None = per_element_residuals
    else:
        precursor_explicit_count = None
        total_atom_count = None
        implicit_count = None
        residuals = None
        balanced = False

    balance = LedgerBalance(
        equation=_BALANCE_EQUATION,
        input_atom_count=input_atom_count,
        reaction_additions=dict(sorted(additions.items())),
        recovered_precursor_explicit_atom_count=precursor_explicit_count,
        reaction_deletions=dict(sorted(deletions.items())),
        residue_atom_count=len(residue_atom_indices),
        guest_atom_count=sum(guest.atom_count for guest in guests),
        unaccounted_atom_count=sum(record.atom_count for record in unaccounted),
        per_element_residuals=residuals,
        balanced=balanced,
    )

    formal_charge_column_present = None
    if identity_evidence is not None:
        raw = identity_evidence.get("formal_charge_column_present")
        formal_charge_column_present = None if raw is None else bool(raw)
    charge_scope = _classify_charge_scope(
        mol=mol,
        recovered_fragments=recovered_fragments,
        guests=guests,
        residue_atom_indices=tuple(residue_atom_indices),
        unaccounted=unaccounted,
        formal_charge_column_present=formal_charge_column_present,
    )
    stereo_scope = _stereo_scope_from_evidence(identity_evidence)

    return AtomLedger(
        family=family,
        input_atom_count=input_atom_count,
        input_elements=input_elements,
        recovered_fragments=recovered_fragments,
        guests=guests,
        residue_atom_indices=tuple(int(idx) for idx in residue_atom_indices),
        residue_elements=residue_elements,
        residue_policy=_RESIDUE_POLICY_BY_FAMILY.get(family),
        unaccounted=unaccounted,
        reaction_additions=dict(sorted(additions.items())),
        reaction_deletions=dict(sorted(deletions.items())),
        recovered_precursor_implicit_hydrogen_count=implicit_count,
        recovered_precursor_total_atom_count=total_atom_count,
        balance=balance,
        charge_scope=charge_scope,
        stereo_scope=stereo_scope,
    )


__all__ = [
    "CHARGE_SCOPE_GUEST_LOCALIZED",
    "CHARGE_SCOPE_LIMITED",
    "CHARGE_SCOPE_NOT_PRESENT",
    "CHARGE_SCOPE_PRESERVED",
    "CHARGE_SCOPE_RESTORED_BY_HEURISTIC",
    "STEREO_SCOPE_NOT_DETECTED",
    "STEREO_SCOPE_UNDETERMINED",
    "STEREO_SCOPE_UNSUPPORTED",
    "AtomLedger",
    "ChargeScopeReport",
    "LedgerBalance",
    "LedgerGuest",
    "LedgerRecoveredFragment",
    "LedgerUnaccountedFragment",
    "StereoScopeReport",
    "build_atom_ledger",
    "build_recovered_fragment_record",
    "element_counts_for_atoms",
    "hill_formula",
    "net_formal_charge_for_atoms",
    "perceive_stereochemistry_evidence",
    "smiles_atom_accounting",
]
