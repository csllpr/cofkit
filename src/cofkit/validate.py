from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping

from .cofid import COFidMonomer, ParsedCOFid, parse_cofid
from .decompose import CifDecompositionResult, decompose_cif_to_cofid
from .lammps import LammpsOptimizationResult, optimize_cif_with_lammps


@dataclass(frozen=True)
class COFidValidationResult:
    status: str
    input_cofid: str
    input_cif: str
    checked_cif: str
    expected_monomers: tuple[COFidMonomer, ...]
    expected_linkage: str
    expected_topology: str
    topology_compared: bool = False
    recovered_cofid: str | None = None
    recovered_monomers: tuple[COFidMonomer, ...] = ()
    recovered_linkage: str | None = None
    recovered_topology: str | None = None
    monomers_match: bool = False
    linkage_match: bool = False
    reason: str | None = None
    decomposition: CifDecompositionResult | None = None
    optimized_cif: str | None = None
    lammps_result: LammpsOptimizationResult | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.status == "match"

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "ok": self.ok,
            "input_cofid": self.input_cofid,
            "input_cif": self.input_cif,
            "checked_cif": self.checked_cif,
            "expected_monomers": [_monomer_to_dict(monomer) for monomer in self.expected_monomers],
            "expected_linkage": self.expected_linkage,
            "expected_topology": self.expected_topology,
            "topology_compared": self.topology_compared,
            "recovered_cofid": self.recovered_cofid,
            "recovered_monomers": [_monomer_to_dict(monomer) for monomer in self.recovered_monomers],
            "recovered_linkage": self.recovered_linkage,
            "recovered_topology": self.recovered_topology,
            "monomers_match": self.monomers_match,
            "linkage_match": self.linkage_match,
            "reason": self.reason,
            "decomposition": None if self.decomposition is None else self.decomposition.to_dict(),
            "optimized_cif": self.optimized_cif,
            "lammps_result": None if self.lammps_result is None else self.lammps_result.to_dict(),
            "metadata": dict(self.metadata),
        }


def validate_cif_against_cofid(
    cofid: str,
    cif_path: str | Path,
) -> COFidValidationResult:
    expected = parse_cofid(cofid)
    input_path = str(Path(cif_path))
    decomposition = decompose_cif_to_cofid(
        cif_path,
        topology=expected.topology,
        linkage=expected.linkage,
        bond_mode="distance",
    )
    return _validation_result_from_decomposition(
        cofid,
        expected,
        input_cif=input_path,
        checked_cif=input_path,
        decomposition=decomposition,
    )


def validate_lammps_optimized_cif_against_cofid(
    cofid: str,
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    lmp_path: str | Path | None = None,
    eqeq_path: str | Path | None = None,
    timeout_seconds: float = 300.0,
    eqeq_timeout_seconds: float | None = 300.0,
) -> COFidValidationResult:
    expected = parse_cofid(cofid)
    lammps_result = optimize_cif_with_lammps(
        cif_path,
        output_dir=output_dir,
        lmp_path=lmp_path,
        eqeq_path=eqeq_path,
        timeout_seconds=timeout_seconds,
        eqeq_timeout_seconds=eqeq_timeout_seconds,
    )
    optimized_cif = lammps_result.optimized_cif
    decomposition = decompose_cif_to_cofid(
        optimized_cif,
        topology=expected.topology,
        linkage=expected.linkage,
        bond_mode="distance",
    )
    return _validation_result_from_decomposition(
        cofid,
        expected,
        input_cif=str(Path(cif_path)),
        checked_cif=optimized_cif,
        decomposition=decomposition,
        optimized_cif=optimized_cif,
        lammps_result=lammps_result,
    )


def _validation_result_from_decomposition(
    cofid: str,
    expected: ParsedCOFid,
    *,
    input_cif: str,
    checked_cif: str,
    decomposition: CifDecompositionResult,
    optimized_cif: str | None = None,
    lammps_result: LammpsOptimizationResult | None = None,
) -> COFidValidationResult:
    base_kwargs = {
        "input_cofid": cofid,
        "input_cif": input_cif,
        "checked_cif": checked_cif,
        "expected_monomers": expected.monomers,
        "expected_linkage": expected.linkage,
        "expected_topology": expected.topology,
        "decomposition": decomposition,
        "optimized_cif": optimized_cif,
        "lammps_result": lammps_result,
        "metadata": {"comparison": "monomers_and_linkage_only"},
    }
    if not decomposition.ok:
        return COFidValidationResult(
            status="failed",
            reason=decomposition.reason or "decomposition did not produce a COFid",
            **base_kwargs,
        )

    try:
        recovered = parse_cofid(str(decomposition.cofid))
    except ValueError as exc:
        return COFidValidationResult(
            status="failed",
            recovered_cofid=decomposition.cofid,
            reason=f"recovered COFid could not be parsed: {exc}",
            **base_kwargs,
        )

    monomers_match = expected.monomers == recovered.monomers
    linkage_match = expected.linkage == recovered.linkage
    ok = monomers_match and linkage_match
    return COFidValidationResult(
        status="match" if ok else "mismatch",
        recovered_cofid=decomposition.cofid,
        recovered_monomers=recovered.monomers,
        recovered_linkage=recovered.linkage,
        recovered_topology=recovered.topology,
        monomers_match=monomers_match,
        linkage_match=linkage_match,
        reason=None if ok else _mismatch_reason(monomers_match=monomers_match, linkage_match=linkage_match),
        **base_kwargs,
    )


def _mismatch_reason(*, monomers_match: bool, linkage_match: bool) -> str:
    mismatches = []
    if not monomers_match:
        mismatches.append("monomers")
    if not linkage_match:
        mismatches.append("linkage")
    return "COFid mismatch in " + " and ".join(mismatches)


def _monomer_to_dict(monomer: COFidMonomer) -> dict[str, object]:
    return {
        "connectivity": monomer.connectivity,
        "reactive_group": monomer.reactive_group,
        "canonical_smiles": monomer.canonical_smiles,
    }


__all__ = [
    "COFidValidationResult",
    "validate_cif_against_cofid",
    "validate_lammps_optimized_cif_against_cofid",
]
