from __future__ import annotations

from dataclasses import dataclass, field
from math import atan2, inf, isfinite
from typing import Callable, Mapping

from ..geometry import Frame, centroid, covariance_eigenpairs, cross, dot, norm, normalize, planar_arrangement_mismatch, scale, sub
from ..model import MonomerSpec, ReactiveMotif
from .motif_registry import MotifKindDefinition, MotifKindRegistry, default_motif_kind_registry

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDepictor
except ImportError:  # pragma: no cover - handled at call sites
    Chem = None
    AllChem = None
    rdDepictor = None


# Heuristic — pending calibration: a molecule at or above this atom count
# whose conformer needed a fallback (random-coordinate or 2D) embedding skips
# force-field minimization entirely, bounding the minimization cost for giant
# precursors on the degraded embedding rungs.
_FALLBACK_FORCEFIELD_ATOM_LIMIT = 180
_NONMETAL_ATOMIC_NUMBERS = frozenset(
    {
        1,  # H
        2,  # He
        5,  # B
        6,  # C
        7,  # N
        8,  # O
        9,  # F
        10,  # Ne
        14,  # Si
        15,  # P
        16,  # S
        17,  # Cl
        18,  # Ar
        33,  # As
        34,  # Se
        35,  # Br
        36,  # Kr
        52,  # Te
        53,  # I
        54,  # Xe
        85,  # At
        86,  # Rn
    }
)


# Heuristic — pending calibration: weights combining the components of
# geometry.planar_arrangement_mismatch into a single conformer shape score for
# shape-aware conformer selection. The components mix units (radians,
# dimensionless relative spread, angstrom RMS); equal weights intentionally
# let any large defect dominate, since the ring-forming node placement cannot
# repair any of them by rigid motion.
_SHAPE_SCORE_ANGULAR_WEIGHT = 1.0
_SHAPE_SCORE_RADIAL_WEIGHT = 1.0
_SHAPE_SCORE_PLANARITY_WEIGHT = 1.0

# Shape-aware conformer selection applies only to precursors placed as rigid
# multi-motif nodes (ring-forming node placement); ditopic precursors are
# edge-fitted, where any conformer span works. Single owner for the builder's
# internal gate and the ring-forming CLI pre-flight check.
SHAPE_SELECTION_MIN_MOTIFS = 3

# Heuristic — pending calibration: upper motif-count bound for shape-aware
# conformer selection in binary-bridge (node-linker) builds. A/B evidence on
# the default monomer library (2026-09-29, out/_ab_shape_selection): trigonal
# (3-motif) nodes are always planar targets and flexible tritopic monomers
# improve (clash-free hcb builds vs hydrogen-clash warnings with energy-only
# selection), while tetrahedral 4-connected monomers on dia get mildly worse
# because no planar conformer exists and the regular-polygon target is wrong
# for them. Planar 4-connected nets (sql/kgm) had no library example to
# validate, so they stay on energy selection until evidence exists.
# Ring-forming builds are exempt: they only support 2D topologies, where the
# planar regular target is always correct, and gate on
# SHAPE_SELECTION_MIN_MOTIFS alone.
SHAPE_SELECTION_VALIDATED_MAX_MOTIFS = 3

# Heuristic — pending calibration: conformer-ensemble floor for shape-aware
# conformer selection. The default 4-8-conformer ensembles frequently contain
# no topology-compatible (regular-planar) conformer for flexible multi-arm
# precursors; 16 gave a near-regular conformer for every triazine/boroxine
# precursor and default-library tritopic monomer tried. Owned here because
# the conformer builder is what consumes the ensemble size; callers
# (ring-forming CLI, batch builds) reference this constant, never retype it.
SHAPE_SELECTION_MIN_CONFORMERS = 16

# Heuristic — pending calibration: upper bound on the heavy-atom RMS distance
# to the best-fit molecular plane for a ditopic/monotopic conformer to carry a
# molecular-plane prior. Calibration evidence (2026-09-29, MMFF/UFF-minimized
# ETKDG conformers): planar aromatic ditopic linkers sit at <= 0.05 angstrom
# RMS (p-phenylenediamine 0.04, terephthalaldehyde ~0.00, 1,4-phenylene
# diboronic acid 0.003), while genuinely nonplanar or twisted conformers sit
# at >= 0.18 angstrom (gauche ethylenediamine 0.18, chair
# trans-1,4-diaminocyclohexane 0.57, twisted biphenyl/stilbene linkers >= 0.4).
_MONOMER_PLANARITY_RMS_TOLERANCE_ANGSTROM = 0.1

# Heuristic — pending calibration: a heavy-atom point set is treated as
# collinear (no unique best-fit plane) when the ratio of the second to the
# first principal RMS width falls below this value. A truly linear molecule
# (dicyanoacetylene) measures ~0.0; the narrowest planar linkers measured
# (phenyl-diacetylene diamine) stay above 0.16.
_MONOMER_PLANE_COLLINEAR_RATIO = 0.1


class AromaticityRestoreError(RuntimeError):
    """RDKit MMFF property generation cleared a molecule's aromaticity.

    ``AllChem.MMFFGetMoleculeProperties`` kekulizes the molecule it is handed
    in place and re-marks aromaticity with MMFF's own perception, which leaves
    cumulene-bridged rings with no aromatic atoms at all. ``Chem.SanitizeMol``
    normally restores them; this error means it did not, so the monomer must
    not be built from a molecule with silently altered aromaticity. Present in
    every RDKit release checked, so it cannot be avoided by upgrading.
    """


@dataclass(frozen=True)
class _MonomerPlaneFit:
    """Molecular-plane assessment for the monomer canonical frame.

    ``status`` is one of:

    - ``"anchor-plane"`` — three or more motif anchors define the connection
      plane directly (the historical multi-motif path);
    - ``"planar"`` — the heavy atoms fit a plane within
      ``_MONOMER_PLANARITY_RMS_TOLERANCE_ANGSTROM``;
    - ``"nonplanar"`` — a least-variance axis exists but the conformer is
      genuinely non-planar, so no plane prior is claimed;
    - ``"collinear"`` — the heavy atoms are essentially one-dimensional, so
      no unique plane exists;
    - ``"degenerate"`` — fewer than two distinct heavy-atom positions.

    ``normal`` is the chemically meaningful plane normal (world frame) and is
    only set for ``anchor-plane``/``planar``. ``canonical_axis`` /
    ``principal_axis`` are covariance axes used purely to canonicalize local
    coordinates when no plane exists; they carry no plane claim.
    """

    status: str
    normal: tuple[float, float, float] | None
    canonical_axis: tuple[float, float, float] | None
    principal_axis: tuple[float, float, float] | None
    rms_deviation: float | None


@dataclass(frozen=True)
class _DetectedMotif:
    reactive_atom_id: int
    anchor_atom_id: int
    atom_ids: tuple[int, ...]
    origin: tuple[float, float, float]
    anchor: tuple[float, float, float]
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class _ConformerEmbeddingResult:
    conformer_ids: tuple[int, ...]
    method: str
    attempts: tuple[Mapping[str, object], ...]

    @property
    def used_fallback(self) -> bool:
        return self.method != "etkdg-v3"


@dataclass(frozen=True)
class _ConformerSelectionResult:
    conformer_id: int
    energy: float | None
    forcefield: str
    optimization_status: str
    diagnostics: tuple[str, ...] = ()
    # Energies of every evaluated conformer (kcal/mol) from the force field
    # that produced this selection; used as the tie-breaker for shape-aware
    # conformer selection. Empty when minimization was skipped.
    conformer_energies: Mapping[int, float] = field(default_factory=dict)


RDKitMatchHandler = Callable[[object, object, tuple[int, ...], MotifKindDefinition], _DetectedMotif | None]
RDKitPostprocessHandler = Callable[[object, tuple[_DetectedMotif, ...], MotifKindDefinition], tuple[_DetectedMotif, ...]]


class RDKitMotifBuilder:
    def __init__(
        self,
        *,
        motif_registry: MotifKindRegistry | None = None,
        match_handlers: Mapping[str, RDKitMatchHandler] | None = None,
        postprocess_handlers: Mapping[str, RDKitPostprocessHandler] | None = None,
    ) -> None:
        self.motif_registry = motif_registry or default_motif_kind_registry()
        self._match_handlers: dict[str, RDKitMatchHandler] = dict(match_handlers or {})
        self._postprocess_handlers: dict[str, RDKitPostprocessHandler] = dict(postprocess_handlers or {})

    def register_match_handler(self, motif_kind: str, handler: RDKitMatchHandler) -> None:
        self._match_handlers[motif_kind] = handler

    def register_postprocess_handler(self, motif_kind: str, handler: RDKitPostprocessHandler) -> None:
        self._postprocess_handlers[motif_kind] = handler

    def supported_motif_kinds(self) -> tuple[str, ...]:
        return tuple(sorted(self._match_handlers))

    def build_monomer(
        self,
        monomer_id: str,
        name: str,
        smiles: str,
        motif_kind: str,
        *,
        num_conformers: int = 8,
        random_seed: int = 0xC0F,
        optimization_max_iterations: int = 500,
        optimization_attempts: int = 1,
        select_conformer_by_motif_shape: bool = False,
    ) -> MonomerSpec:
        if Chem is None or AllChem is None:
            raise RuntimeError("RDKit is required for build_rdkit_monomer()")
        definition = self.motif_registry.get(motif_kind)
        if definition.rdkit_smarts is None:
            raise ValueError(f"motif kind {motif_kind!r} has no RDKit SMARTS configuration")

        base = Chem.MolFromSmiles(smiles)
        if base is None:
            raise ValueError(f"RDKit could not parse SMILES for {monomer_id!r}")
        molecule = Chem.AddHs(base)
        try:
            embedding = _embed_conformers(
                molecule,
                num_conformers=num_conformers,
                random_seed=random_seed,
            )
        except ValueError as exc:
            raise ValueError(
                f"RDKit conformer generation failed for monomer {monomer_id!r} "
                f"from SMILES {smiles!r}: {exc}"
            ) from exc
        optimization_skip_reason = None
        if embedding.method == "rdkit-2d":
            optimization_skip_reason = "planar coordinate fallback is not force-field minimized"
        elif embedding.used_fallback and molecule.GetNumAtoms() >= _FALLBACK_FORCEFIELD_ATOM_LIMIT:
            optimization_skip_reason = (
                f"fallback conformer has {molecule.GetNumAtoms()} atoms; "
                f"force-field minimization is capped below {_FALLBACK_FORCEFIELD_ATOM_LIMIT} atoms"
            )
        try:
            selection = _optimize_conformers(
                molecule,
                embedding.conformer_ids,
                skip_reason=optimization_skip_reason,
                max_iterations=optimization_max_iterations,
                max_attempts=optimization_attempts,
            )
        except AromaticityRestoreError as exc:
            raise AromaticityRestoreError(f"monomer {monomer_id!r}: {exc}") from exc
        conformer_id = selection.conformer_id
        conformer = molecule.GetConformer(conformer_id)

        detected = self._detect_motifs(molecule, conformer, definition)
        if not detected:
            raise ValueError(f"no {motif_kind!r} motifs detected in {monomer_id!r}")

        conformer_selection = "energy"
        shape_score: float | None = None
        if select_conformer_by_motif_shape and len(detected) >= SHAPE_SELECTION_MIN_MOTIFS:
            shape_choice = self._select_conformer_by_motif_shape(
                molecule,
                embedding.conformer_ids,
                definition,
                reference_motif_count=len(detected),
                energies=selection.conformer_energies,
            )
            if shape_choice is not None:
                conformer_id, shape_score = shape_choice
                if conformer_id != selection.conformer_id:
                    conformer = molecule.GetConformer(conformer_id)
                    detected = self._detect_motifs(molecule, conformer, definition)
                conformer_selection = "motif_shape"

        atom_positions, motifs, plane_fit = _build_geometry(detected, molecule, conformer, definition)
        atom_symbols = tuple(atom.GetSymbol() for atom in molecule.GetAtoms())
        bonds = tuple(
            (
                int(bond.GetBeginAtomIdx()),
                int(bond.GetEndAtomIdx()),
                float(bond.GetBondTypeAsDouble()),
            )
            for bond in molecule.GetBonds()
        )
        return MonomerSpec(
            id=monomer_id,
            name=name,
            motifs=motifs,
            conformer_ids=(f"rdkit-conf-{conformer_id}",),
            atom_symbols=atom_symbols,
            atom_positions=atom_positions,
            bonds=bonds,
            metadata={
                "source_smiles": smiles,
                "geometry_mode": {
                    "etkdg-v3": "rdkit-etkdg",
                    "etkdg-v3-random": "rdkit-etkdg-random-coordinates",
                    "etkdg-v2-random": "rdkit-etkdg-v2-random-coordinates",
                    "rdkit-2d": "rdkit-2d-coordinate-fallback",
                }[embedding.method],
                "embedding_method": embedding.method,
                "embedding_fallback": embedding.used_fallback,
                "embedding_attempts": [dict(attempt) for attempt in embedding.attempts],
                "motif_detection": f"rdkit-smarts:{definition.rdkit_smarts}",
                "motif_kind": motif_kind,
                "n_atoms": molecule.GetNumAtoms(),
                "n_heavy_atoms": base.GetNumAtoms(),
                "n_conformers": len(embedding.conformer_ids),
                "selected_conformer_id": conformer_id,
                "conformer_selection": conformer_selection,
                **(
                    {"selected_conformer_shape_score": shape_score}
                    if shape_score is not None
                    else {}
                ),
                "forcefield": selection.forcefield,
                "forcefield_optimization_status": selection.optimization_status,
                "forcefield_diagnostics": selection.diagnostics,
                "selected_conformer_energy": selection.conformer_energies.get(conformer_id, selection.energy),
                "plane_normal": plane_fit.normal,
                # Plane-prior provenance: "anchor-plane" (>=3 motif anchors
                # define the connection plane), "planar" (heavy-atom best-fit
                # plane within tolerance), or the uncertain statuses
                # "nonplanar" / "collinear" / "degenerate", for which
                # plane_normal is None and motif frame normals are zero —
                # consumers must skip plane-dependent terms rather than fall
                # back to a fabricated plane (impact-review claim T1-8).
                "plane_status": plane_fit.status,
                "plane_rms_deviation": plane_fit.rms_deviation,
            },
        )

    def _select_conformer_by_motif_shape(
        self,
        molecule,
        conformer_ids: tuple[int, ...],
        definition: MotifKindDefinition,
        *,
        reference_motif_count: int,
        energies: Mapping[int, float],
    ) -> tuple[int, float] | None:
        """Pick the embedded conformer whose motif origins best form a regular
        planar polygon — the shape the ring-forming node placement assumes.

        Returns ``(conformer_id, shape_score)`` minimizing
        ``(shape_score, energy, conformer_id)``, or ``None`` when no conformer
        yields a measurable arrangement (motif count drift or degenerate
        geometry); callers then keep the energy-selected conformer.
        """
        best: tuple[tuple[float, float, int], int, float] | None = None
        for conf_id in conformer_ids:
            conformer = molecule.GetConformer(conf_id)
            detected = self._detect_motifs(molecule, conformer, definition)
            if len(detected) != reference_motif_count:
                continue
            mismatch = planar_arrangement_mismatch(tuple(motif.origin for motif in detected))
            if mismatch is None:
                continue
            score = (
                _SHAPE_SCORE_ANGULAR_WEIGHT * mismatch.angular_max_radians
                + _SHAPE_SCORE_RADIAL_WEIGHT * mismatch.radial_spread
                + _SHAPE_SCORE_PLANARITY_WEIGHT * mismatch.planarity_rms
            )
            key = (score, energies.get(conf_id, inf), conf_id)
            if best is None or key < best[0]:
                best = (key, conf_id, score)
        if best is None:
            return None
        return best[1], best[2]

    def _detect_motifs(self, molecule, conformer, definition: MotifKindDefinition) -> tuple[_DetectedMotif, ...]:
        assert definition.rdkit_smarts is not None
        pattern = Chem.MolFromSmarts(definition.rdkit_smarts)
        if pattern is None:
            raise ValueError(f"invalid SMARTS for motif kind {definition.kind!r}")
        try:
            handler = self._match_handlers[definition.kind]
        except KeyError as exc:
            raise ValueError(f"motif kind {definition.kind!r} has no RDKit match handler") from exc
        matches = molecule.GetSubstructMatches(pattern, uniquify=False)
        detected_by_atoms: dict[tuple[int, ...], _DetectedMotif] = {}
        for match in matches:
            detected = handler(molecule, conformer, tuple(int(atom_id) for atom_id in match), definition)
            if detected is None:
                continue
            detected_by_atoms.setdefault(tuple(sorted(detected.atom_ids)), detected)
        return tuple(
            self._postprocess_detected_motifs(
                molecule,
                tuple(
                    sorted(
                        detected_by_atoms.values(),
                        key=lambda item: (item.reactive_atom_id, item.anchor_atom_id, item.atom_ids),
                    )
                ),
                definition,
            )
        )

    def _postprocess_detected_motifs(
        self,
        molecule,
        detected: tuple[_DetectedMotif, ...],
        definition: MotifKindDefinition,
    ) -> tuple[_DetectedMotif, ...]:
        handler = self._postprocess_handlers.get(definition.kind)
        if handler is None:
            return detected
        return handler(molecule, detected, definition)

    @classmethod
    def builtin(cls, *, motif_registry: MotifKindRegistry | None = None) -> "RDKitMotifBuilder":
        builder = cls(motif_registry=motif_registry)
        builder.register_match_handler("amine", _interpret_primary_amine_match)
        builder.register_match_handler("aldehyde", _interpret_aldehyde_match)
        builder.register_match_handler("hydrazine", _interpret_hydrazine_match)
        builder.register_match_handler("hydrazide", _interpret_hydrazide_match)
        builder.register_match_handler("boronic_acid", _interpret_boronic_acid_match)
        builder.register_match_handler("catechol", _interpret_catechol_match)
        builder.register_match_handler("keto_aldehyde", _interpret_keto_aldehyde_match)
        builder.register_match_handler("activated_methylene", _interpret_activated_methylene_match)
        builder.register_match_handler("nitrile", _interpret_nitrile_match)
        builder.register_postprocess_handler("keto_aldehyde", _postprocess_keto_aldehyde_matches)
        return builder


def build_rdkit_monomer(
    monomer_id: str,
    name: str,
    smiles: str,
    motif_kind: str,
    *,
    num_conformers: int = 8,
    random_seed: int = 0xC0F,
    optimization_max_iterations: int = 500,
    optimization_attempts: int = 1,
    select_conformer_by_motif_shape: bool = False,
    motif_registry: MotifKindRegistry | None = None,
    builder: RDKitMotifBuilder | None = None,
) -> MonomerSpec:
    effective_builder = builder or RDKitMotifBuilder.builtin(motif_registry=motif_registry)
    return effective_builder.build_monomer(
        monomer_id,
        name,
        smiles,
        motif_kind,
        num_conformers=num_conformers,
        random_seed=random_seed,
        optimization_max_iterations=optimization_max_iterations,
        optimization_attempts=optimization_attempts,
        select_conformer_by_motif_shape=select_conformer_by_motif_shape,
    )


def detect_rdkit_motif_count(
    smiles: str,
    motif_kind: str,
    *,
    motif_registry: MotifKindRegistry | None = None,
    builder: RDKitMotifBuilder | None = None,
) -> int:
    """Return the number of buildable motifs without generating a conformer.

    Decomposition uses this lightweight path to verify that a repaired precursor
    really exposes the role and connectivity written to its COFid.  Motif match
    handlers need a conformer for their geometry payloads, but motif acceptance
    itself is graph-based, so a zero-coordinate conformer is sufficient here.
    """

    if Chem is None:
        raise RuntimeError("RDKit is required for motif detection")
    effective_builder = builder or RDKitMotifBuilder.builtin(motif_registry=motif_registry)
    definition = effective_builder.motif_registry.get(motif_kind)
    if definition.rdkit_smarts is None:
        raise ValueError(f"motif kind {motif_kind!r} has no RDKit SMARTS configuration")
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise ValueError(f"RDKit could not parse SMILES {smiles!r}")
    molecule = Chem.AddHs(base)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    molecule.AddConformer(conformer, assignId=True)
    return len(effective_builder._detect_motifs(molecule, conformer, definition))


def monomer_geometry_degradation_warnings(monomer: MonomerSpec) -> tuple[str, ...]:
    """Describe any geometry fallback rungs recorded in monomer metadata.

    The conformer builder degrades through documented fallback ladders —
    embedding: ETKDGv3 -> random-coordinate retries -> 2D planar depiction;
    force-field optimization: MMFF -> UFF -> unconverged -> unminimized — and
    records the rung in ``metadata`` (``embedding_method``,
    ``embedding_fallback``, ``forcefield``, ``forcefield_optimization_status``,
    ``forcefield_diagnostics``). Nothing downstream consumed those flags, so a
    2D-planar or unminimized monomer was indistinguishable from an
    MMFF-optimized one. This surfaces each off-top rung as a
    ``monomer_geometry_degraded:<monomer_id>:<detail>`` warning string; a
    top-rung monomer yields an empty tuple. Metadata-bearing monomers built by
    other means simply produce no warnings when the keys are absent.
    """
    metadata = monomer.metadata
    details: list[str] = []
    embedding_method = metadata.get("embedding_method")
    if embedding_method == "rdkit-2d":
        details.append("conformer embedding fell back to a 2D planar depiction (rdkit-2d)")
    elif metadata.get("embedding_fallback"):
        details.append(f"conformer embedding fell back to random-coordinate embedding ({embedding_method})")
    forcefield = metadata.get("forcefield")
    status = metadata.get("forcefield_optimization_status")
    if status == "unconverged" and forcefield == "none":
        details.append(
            "no supported force field produced a minimized conformer; using the unminimized embedded conformer"
        )
    elif status == "unconverged":
        details.append(
            f"force-field optimization did not converge; using the lowest-energy unconverged {forcefield} conformer"
        )
    elif status == "skipped" or forcefield == "none":
        diagnostics = tuple(str(item) for item in metadata.get("forcefield_diagnostics", ()) or ())
        reason = f" ({diagnostics[0]})" if diagnostics else ""
        details.append(f"conformer is not force-field minimized{reason}")
    elif forcefield == "UFF":
        details.append("force-field optimization fell back to UFF because MMFF parameters were unavailable")
    return tuple(f"monomer_geometry_degraded:{monomer.id}:{detail}" for detail in details)


def _embed_conformers(
    molecule,
    *,
    num_conformers: int,
    random_seed: int,
) -> _ConformerEmbeddingResult:
    attempts: list[Mapping[str, object]] = []
    requested = max(1, num_conformers)
    for method, parameter_factory, use_random_coords in (
        ("etkdg-v3", AllChem.ETKDGv3, False),
        ("etkdg-v3-random", AllChem.ETKDGv3, True),
        ("etkdg-v2-random", AllChem.ETKDGv2, True),
    ):
        molecule.RemoveAllConformers()
        params = parameter_factory()
        params.randomSeed = random_seed
        # Cited: 0.2 A is RDKit's own conventional pruneRmsThresh, used in the
        # ETKDG documentation and EmbedMultipleConfs examples.
        params.pruneRmsThresh = 0.2
        if method != "etkdg-v2-random":
            params.useSmallRingTorsions = True
        params.useRandomCoords = use_random_coords
        if use_random_coords:
            # Bound the retry cost for giant precursors (heuristic — pending
            # calibration: 50 iterations per attempt). ETKDGv2 is the next
            # fallback when v3 exhausts these attempts.
            params.maxIterations = 50
        try:
            raw_ids = AllChem.EmbedMultipleConfs(
                molecule,
                numConfs=requested,
                params=params,
            )
            conformer_ids = tuple(int(conf_id) for conf_id in raw_ids)
        except RuntimeError as exc:
            attempts.append(
                {
                    "method": method,
                    "status": "error",
                    "error": f"{type(exc).__name__}: {str(exc).splitlines()[0]}",
                }
            )
            continue
        geometry_error = _conformer_coordinate_error(molecule, conformer_ids)
        if conformer_ids and geometry_error is None:
            attempts.append({"method": method, "status": "success", "n_conformers": len(conformer_ids)})
            return _ConformerEmbeddingResult(conformer_ids, method, tuple(attempts))
        attempts.append(
            {
                "method": method,
                "status": "failed",
                "error": geometry_error or "RDKit returned no conformers",
            }
        )

    planar_fallback_error = _planar_coordinate_fallback_error(molecule)
    if planar_fallback_error is not None:
        attempts.append(
            {
                "method": "rdkit-2d",
                "status": "skipped",
                "error": planar_fallback_error,
            }
        )
    else:
        molecule.RemoveAllConformers()
        try:
            conformer_id = int(rdDepictor.Compute2DCoords(molecule))
            conformer_ids = (conformer_id,)
            geometry_error = _conformer_coordinate_error(molecule, conformer_ids)
        except (RuntimeError, ValueError) as exc:
            attempts.append(
                {
                    "method": "rdkit-2d",
                    "status": "error",
                    "error": f"{type(exc).__name__}: {str(exc).splitlines()[0]}",
                }
            )
        else:
            if geometry_error is None:
                attempts.append({"method": "rdkit-2d", "status": "success", "n_conformers": 1})
                return _ConformerEmbeddingResult(conformer_ids, "rdkit-2d", tuple(attempts))
            attempts.append({"method": "rdkit-2d", "status": "failed", "error": geometry_error})

    attempt_summary = "; ".join(
        f"{attempt['method']}: {attempt.get('error', attempt['status'])}"
        for attempt in attempts
    )
    raise ValueError(f"RDKit conformer embedding failed after all fallbacks ({attempt_summary})")


def _planar_coordinate_fallback_error(molecule) -> str | None:
    """Return why a last-resort 2D conformer would be unsafe, if applicable."""
    atoms = tuple(molecule.GetAtoms())
    has_coordination_or_charge = any(
        atom.GetAtomicNum() not in _NONMETAL_ATOMIC_NUMBERS
        or atom.GetFormalCharge() != 0
        for atom in atoms
    )
    if not has_coordination_or_charge:
        return "reserved for metal-containing or formally charged precursors"

    allowed_hybridizations = {
        Chem.HybridizationType.SP,
        Chem.HybridizationType.SP2,
    }
    nonplanar_atoms = tuple(
        atom.GetIdx()
        for atom in atoms
        if atom.GetAtomicNum() > 1
        and atom.GetAtomicNum() in _NONMETAL_ATOMIC_NUMBERS
        and not atom.GetIsAromatic()
        and atom.GetHybridization() not in allowed_hybridizations
    )
    if nonplanar_atoms:
        preview = ", ".join(str(index) for index in nonplanar_atoms[:5])
        suffix = ", ..." if len(nonplanar_atoms) > 5 else ""
        return f"non-planar heavy-atom hybridization at atom indices {preview}{suffix}"
    return None


def _conformer_coordinate_error(molecule, conformer_ids: tuple[int, ...]) -> str | None:
    if not conformer_ids:
        return "no conformers were generated"
    for conformer_id in conformer_ids:
        conformer = molecule.GetConformer(conformer_id)
        positions = tuple(conformer.GetAtomPosition(index) for index in range(molecule.GetNumAtoms()))
        if any(not all(isfinite(value) for value in (point.x, point.y, point.z)) for point in positions):
            return f"conformer {conformer_id} contains non-finite coordinates"
        if len(positions) > 1:
            extent = max(
                max(point.x for point in positions) - min(point.x for point in positions),
                max(point.y for point in positions) - min(point.y for point in positions),
                max(point.z for point in positions) - min(point.z for point in positions),
            )
            if extent < 1.0e-6:
                return f"conformer {conformer_id} has collapsed coordinates"
    return None


def _minimize_conformer(field, *, max_iterations: int, max_attempts: int) -> int:
    for _ in range(max_attempts):
        status = field.Minimize(maxIts=max_iterations)
        if status != 1:
            return status
    return status


def _aromaticity_snapshot(molecule) -> tuple[frozenset[int], frozenset[int]]:
    """Return the aromatic atom-index and aromatic bond-index sets."""
    aromatic_atoms = frozenset(
        atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetIsAromatic()
    )
    aromatic_bonds = frozenset(
        bond.GetIdx()
        for bond in molecule.GetBonds()
        if bond.GetBondType() == Chem.BondType.AROMATIC
    )
    return aromatic_atoms, aromatic_bonds


def _restore_mmff_aromaticity(
    molecule,
    snapshot: tuple[frozenset[int], frozenset[int]],
) -> str | None:
    """Repair aromaticity lost to RDKit's in-place MMFF kekulization.

    Returns ``None`` when the aromaticity state is unchanged, a diagnostic
    entry describing a successful repair, and raises
    :class:`AromaticityRestoreError` when ``Chem.SanitizeMol`` fails or does
    not restore the snapshotted sets. The comparison is set equality, never a
    count: MMFF can re-mark the same number of atoms or bonds differently,
    which a scalar count would not notice.
    """
    if _aromaticity_snapshot(molecule) == snapshot:
        return None
    atoms_before = len(snapshot[0])
    try:
        Chem.SanitizeMol(molecule)
    except Exception as exc:
        raise AromaticityRestoreError(
            "RDKit MMFF property generation cleared aromaticity in place and "
            f"Chem.SanitizeMol could not restore it: {type(exc).__name__}: {exc}"
        ) from exc
    after = _aromaticity_snapshot(molecule)
    if after != snapshot:
        raise AromaticityRestoreError(
            "RDKit MMFF property generation cleared aromaticity in place and "
            f"Chem.SanitizeMol did not restore it ({atoms_before} aromatic atoms "
            f"before, {len(after[0])} after)"
        )
    return (
        "MMFF property generation cleared RDKit aromatic flags in place; "
        f"re-sanitized and restored {atoms_before} aromatic atoms"
    )


def _optimize_conformers(
    molecule,
    conformer_ids: tuple[int, ...],
    *,
    skip_reason: str | None = None,
    max_iterations: int = 500,
    max_attempts: int = 1,
) -> _ConformerSelectionResult:
    if any(not isinstance(v, int) or isinstance(v, bool) or v <= 0 for v in (max_iterations, max_attempts)):
        raise ValueError("Conformer minimization iterations and attempts must be positive integers.")
    if skip_reason is not None:
        return _ConformerSelectionResult(
            conformer_id=conformer_ids[0],
            energy=None,
            forcefield="none",
            optimization_status="skipped",
            diagnostics=(skip_reason,),
        )

    diagnostics: list[str] = []
    best_unconverged: tuple[int, float, str] | None = None
    try:
        has_mmff_parameters = bool(AllChem.MMFFHasAllMoleculeParams(molecule))
    except (RuntimeError, ValueError) as exc:
        diagnostics.append(f"MMFF parameter check failed: {type(exc).__name__}: {str(exc).splitlines()[0]}")
        has_mmff_parameters = False
    if has_mmff_parameters:
        # MMFF property generation kekulizes the molecule in place and can drop
        # the aromatic perception the motif detection downstream depends on, so
        # the state is snapshotted here and repaired once MMFF is done with the
        # molecule, never in between force-field construction calls.
        aromaticity_snapshot = _aromaticity_snapshot(molecule)
        try:
            props = AllChem.MMFFGetMoleculeProperties(molecule)
        except (RuntimeError, ValueError) as exc:
            diagnostics.append(
                f"MMFF property generation failed: {type(exc).__name__}: {str(exc).splitlines()[0]}"
            )
            props = None
        best = None
        mmff_energies: dict[int, float] = {}
        if props is not None:
            for conf_id in conformer_ids:
                try:
                    field = AllChem.MMFFGetMoleculeForceField(molecule, props, confId=conf_id)
                    if field is None:
                        continue
                    status = _minimize_conformer(field, max_iterations=max_iterations, max_attempts=max_attempts)
                    energy = float(field.CalcEnergy())
                except (RuntimeError, ValueError) as exc:
                    diagnostics.append(
                        f"MMFF conformer {conf_id} failed: {type(exc).__name__}: {str(exc).splitlines()[0]}"
                    )
                    continue
                if isfinite(energy):
                    mmff_energies[conf_id] = energy
                if status != 0:
                    diagnostics.append(f"Conformer {conf_id} did not converge (status {status})")
                    if isfinite(energy) and (best_unconverged is None or energy < best_unconverged[1]):
                        best_unconverged = (conf_id, energy, "MMFF")
                if status == 0 and isfinite(energy) and (best is None or energy < best[1]):
                    best = (conf_id, energy)
        aromaticity_diagnostic = _restore_mmff_aromaticity(molecule, aromaticity_snapshot)
        if aromaticity_diagnostic is not None:
            diagnostics.append(aromaticity_diagnostic)
        if best is not None:
            return _ConformerSelectionResult(
                best[0],
                best[1],
                "MMFF",
                "optimized",
                tuple(diagnostics),
                conformer_energies=mmff_energies,
            )

    best = None
    uff_energies: dict[int, float] = {}
    for conf_id in conformer_ids:
        try:
            field = AllChem.UFFGetMoleculeForceField(molecule, confId=conf_id)
            if field is None:
                continue
            status = _minimize_conformer(field, max_iterations=max_iterations, max_attempts=max_attempts)
            energy = float(field.CalcEnergy())
        except (RuntimeError, ValueError) as exc:
            diagnostics.append(
                f"UFF conformer {conf_id} failed: {type(exc).__name__}: {str(exc).splitlines()[0]}"
            )
            continue
        if isfinite(energy):
            uff_energies[conf_id] = energy
        if status != 0:
            diagnostics.append(f"Conformer {conf_id} did not converge (status {status})")
            if isfinite(energy) and (best_unconverged is None or energy < best_unconverged[1]):
                best_unconverged = (conf_id, energy, "UFF")
        if status == 0 and isfinite(energy) and (best is None or energy < best[1]):
            best = (conf_id, energy)
    if best is not None:
        return _ConformerSelectionResult(
            best[0],
            best[1],
            "UFF",
            "optimized",
            tuple(diagnostics),
            conformer_energies=uff_energies,
        )
    if best_unconverged is not None:
        diagnostics.append(
            "no force-field minimization fully converged; proceeding with the "
            f"lowest-energy unconverged {best_unconverged[2]} conformer"
        )
        return _ConformerSelectionResult(
            best_unconverged[0],
            best_unconverged[1],
            best_unconverged[2],
            "unconverged",
            tuple(diagnostics),
            conformer_energies=mmff_energies if best_unconverged[2] == "MMFF" else uff_energies,
        )
    diagnostics.append(
        "no supported force field produced a finite minimized conformer; "
        "proceeding with the unminimized embedded conformer"
    )
    return _ConformerSelectionResult(
        conformer_ids[0],
        None,
        "none",
        "unconverged",
        tuple(diagnostics),
    )


def _build_geometry(
    detected: tuple[_DetectedMotif, ...],
    molecule,
    conformer,
    definition: MotifKindDefinition,
) -> tuple[tuple[tuple[float, float, float], ...], tuple[ReactiveMotif, ...], _MonomerPlaneFit]:
    points = tuple(_conformer_point(conformer, atom.GetIdx()) for atom in molecule.GetAtoms())
    center = centroid(points)
    plane_fit = _plane_normal(detected, molecule, conformer)
    # The z axis of the canonical local frame is the molecular plane normal
    # when one exists; otherwise the covariance least-variance axis keeps the
    # coordinates deterministic and molecule-derived without claiming a plane.
    # Only a truly degenerate point set falls back to a world axis.
    z_axis = plane_fit.canonical_axis if plane_fit.canonical_axis is not None else (0.0, 0.0, 1.0)
    motif_normal = (0.0, 0.0, 1.0) if plane_fit.normal is not None else (0.0, 0.0, 0.0)
    provisional_primary = _project_onto_plane(sub(detected[0].origin, detected[0].anchor), z_axis)
    if _is_near_zero(provisional_primary):
        provisional_primary = plane_fit.principal_axis if plane_fit.principal_axis is not None else (1.0, 0.0, 0.0)
    x_axis = normalize(provisional_primary)
    y_axis = normalize(cross(z_axis, x_axis))

    def transform(point: tuple[float, float, float]) -> tuple[float, float, float]:
        shifted = sub(point, center)
        return (
            dot(shifted, x_axis),
            dot(shifted, y_axis),
            dot(shifted, z_axis),
        )

    atom_positions = tuple(transform(point) for point in points)
    motif_rows: list[tuple[float, ReactiveMotif]] = []
    for detected_motif in detected:
        origin = transform(detected_motif.origin)
        anchor = transform(detected_motif.anchor)
        primary = _project_onto_plane(sub(origin, anchor), (0.0, 0.0, 1.0))
        if _is_near_zero(primary):
            primary = (1.0, 0.0, 0.0)
        frame = Frame(origin=origin, primary=normalize(primary), normal=motif_normal)
        angle = atan2(origin[1], origin[0])
        motif_rows.append(
            (
                angle,
                ReactiveMotif(
                    id=f"{definition.id_prefix}{len(motif_rows) + 1}",
                    kind=definition.kind,
                    atom_ids=detected_motif.atom_ids,
                    frame=frame,
                    allowed_reaction_templates=definition.allowed_reaction_templates,
                    metadata={
                        **detected_motif.metadata,
                        "reactive_atom_id": detected_motif.reactive_atom_id,
                        "anchor_atom_id": detected_motif.anchor_atom_id,
                    },
                ),
            )
        )

    motifs = tuple(
        ReactiveMotif(
            id=f"{definition.id_prefix}{index}",
            kind=motif.kind,
            atom_ids=motif.atom_ids,
            frame=motif.frame,
            valence=motif.valence,
            symmetry_order=motif.symmetry_order,
            planarity_class=motif.planarity_class,
            allowed_reaction_templates=motif.allowed_reaction_templates,
            metadata=motif.metadata,
        )
        for index, (_, motif) in enumerate(sorted(motif_rows, key=lambda item: item[0]), start=1)
    )
    return atom_positions, motifs, plane_fit


def _interpret_primary_amine_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    reactive_atom_id, anchor_atom_id = int(match[0]), int(match[1])
    atom_ids = [reactive_atom_id, anchor_atom_id]
    atom_ids.extend(
        neighbor.GetIdx()
        for neighbor in molecule.GetAtomWithIdx(reactive_atom_id).GetNeighbors()
        if neighbor.GetAtomicNum() == 1
    )
    atom_ids = sorted(set(atom_ids))
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
    )


def _interpret_aldehyde_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    reactive_atom_id, oxygen_atom_id, anchor_atom_id = (int(match[0]), int(match[1]), int(match[2]))
    atom_ids = [reactive_atom_id, oxygen_atom_id, anchor_atom_id]
    atom_ids.extend(
        neighbor.GetIdx()
        for neighbor in molecule.GetAtomWithIdx(reactive_atom_id).GetNeighbors()
        if neighbor.GetAtomicNum() == 1
    )
    atom_ids = sorted(set(atom_ids))
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
    )


def _interpret_hydrazine_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif | None:
    del definition
    reactive_atom_id, anchor_atom_id = (int(match[0]), int(match[1]))
    hydrogen_atom_ids = _attached_hydrogen_ids(molecule, reactive_atom_id)
    if len(hydrogen_atom_ids) < 2:
        return None
    atom_ids = sorted({reactive_atom_id, anchor_atom_id, *hydrogen_atom_ids})
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
        metadata={
            "internal_nitrogen_atom_id": anchor_atom_id,
            "hydrogen_atom_ids": hydrogen_atom_ids,
        },
    )


def _interpret_hydrazide_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    terminal_nitrogen_atom_id, internal_nitrogen_atom_id, carbonyl_carbon_atom_id, carbonyl_oxygen_atom_id, anchor_atom_id = (
        int(match[0]),
        int(match[1]),
        int(match[2]),
        int(match[3]),
        int(match[4]),
    )
    hydrogen_atom_ids = _attached_hydrogen_ids(molecule, terminal_nitrogen_atom_id)
    atom_ids = sorted(
        {
            terminal_nitrogen_atom_id,
            internal_nitrogen_atom_id,
            carbonyl_carbon_atom_id,
            carbonyl_oxygen_atom_id,
            anchor_atom_id,
            *hydrogen_atom_ids,
        }
    )
    return _DetectedMotif(
        reactive_atom_id=terminal_nitrogen_atom_id,
        anchor_atom_id=internal_nitrogen_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, terminal_nitrogen_atom_id),
        anchor=_conformer_point(conformer, internal_nitrogen_atom_id),
        metadata={
            "internal_nitrogen_atom_id": internal_nitrogen_atom_id,
            "carbonyl_carbon_atom_id": carbonyl_carbon_atom_id,
            "carbonyl_oxygen_atom_id": carbonyl_oxygen_atom_id,
            "hydrogen_atom_ids": hydrogen_atom_ids,
        },
    )


def _interpret_boronic_acid_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    anchor_atom_id, boron_atom_id, oxygen_atom_id_1, oxygen_atom_id_2 = (
        int(match[0]),
        int(match[1]),
        int(match[2]),
        int(match[3]),
    )
    hydrogen_atom_ids = tuple(
        hydrogen_atom_id
        for oxygen_atom_id in (oxygen_atom_id_1, oxygen_atom_id_2)
        for hydrogen_atom_id in _attached_hydrogen_ids(molecule, oxygen_atom_id)
    )
    atom_ids = sorted(
        {
            anchor_atom_id,
            boron_atom_id,
            oxygen_atom_id_1,
            oxygen_atom_id_2,
            *hydrogen_atom_ids,
        }
    )
    return _DetectedMotif(
        reactive_atom_id=boron_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, boron_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
        metadata={
            "oxygen_atom_ids": (oxygen_atom_id_1, oxygen_atom_id_2),
            "hydrogen_atom_ids": hydrogen_atom_ids,
        },
    )


def _interpret_nitrile_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    anchor_atom_id, carbon_atom_id, nitrogen_atom_id = (
        int(match[0]),
        int(match[1]),
        int(match[2]),
    )
    return _DetectedMotif(
        reactive_atom_id=carbon_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=(anchor_atom_id, carbon_atom_id, nitrogen_atom_id),
        origin=_conformer_point(conformer, carbon_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
        metadata={"nitrogen_atom_id": nitrogen_atom_id},
    )


def _interpret_catechol_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif | None:
    del definition
    oxygen_atom_id_1, carbon_atom_id_1 = (int(match[0]), int(match[1]))
    neighbor_pairs = []
    for aromatic_neighbor in molecule.GetAtomWithIdx(carbon_atom_id_1).GetNeighbors():
        carbon_atom_id_2 = int(aromatic_neighbor.GetIdx())
        if carbon_atom_id_2 == oxygen_atom_id_1:
            continue
        if aromatic_neighbor.GetAtomicNum() != 6 or not aromatic_neighbor.GetIsAromatic():
            continue
        oxygen_atom_id_2 = _aromatic_hydroxyl_substituent(molecule, carbon_atom_id_2, excluded_atom_id=carbon_atom_id_1)
        if oxygen_atom_id_2 is None:
            continue
        neighbor_pairs.append((carbon_atom_id_2, oxygen_atom_id_2))
    if not neighbor_pairs:
        return None
    carbon_atom_id_2, oxygen_atom_id_2 = min(neighbor_pairs)
    if carbon_atom_id_1 > carbon_atom_id_2:
        return None
    hydrogen_atom_ids_1 = _attached_hydrogen_ids(molecule, oxygen_atom_id_1)
    hydrogen_atom_ids_2 = _attached_hydrogen_ids(molecule, oxygen_atom_id_2)
    if len(hydrogen_atom_ids_1) != 1 or len(hydrogen_atom_ids_2) != 1:
        return None
    hydrogen_atom_ids = hydrogen_atom_ids_1 + hydrogen_atom_ids_2
    origin = centroid(
        (
            _conformer_point(conformer, oxygen_atom_id_1),
            _conformer_point(conformer, oxygen_atom_id_2),
        )
    )
    anchor = centroid(
        (
            _conformer_point(conformer, carbon_atom_id_1),
            _conformer_point(conformer, carbon_atom_id_2),
        )
    )
    atom_ids = sorted(
        {
            oxygen_atom_id_1,
            oxygen_atom_id_2,
            carbon_atom_id_1,
            carbon_atom_id_2,
            *hydrogen_atom_ids,
        }
    )
    return _DetectedMotif(
        reactive_atom_id=min(oxygen_atom_id_1, oxygen_atom_id_2),
        anchor_atom_id=min(carbon_atom_id_1, carbon_atom_id_2),
        atom_ids=tuple(atom_ids),
        origin=origin,
        anchor=anchor,
        metadata={
            "reactive_atom_ids": (oxygen_atom_id_1, oxygen_atom_id_2),
            "anchor_atom_ids": (carbon_atom_id_1, carbon_atom_id_2),
            "hydrogen_atom_ids": hydrogen_atom_ids,
        },
    )


def _interpret_keto_aldehyde_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif | None:
    del definition
    reactive_atom_id, oxygen_atom_id, anchor_atom_id = (int(match[0]), int(match[1]), int(match[2]))
    ortho_hydroxyl = _find_ortho_hydroxyl(molecule, anchor_atom_id)
    if ortho_hydroxyl is not None:
        hydroxyl_oxygen_atom_id, hydroxyl_hydrogen_atom_id, hydroxyl_anchor_atom_id = ortho_hydroxyl
        atom_ids = [
            reactive_atom_id,
            oxygen_atom_id,
            anchor_atom_id,
            hydroxyl_oxygen_atom_id,
            hydroxyl_hydrogen_atom_id,
        ]
        atom_ids.extend(_attached_hydrogen_ids(molecule, reactive_atom_id))
        atom_ids = sorted(set(atom_ids))
        return _DetectedMotif(
            reactive_atom_id=reactive_atom_id,
            anchor_atom_id=anchor_atom_id,
            atom_ids=tuple(atom_ids),
            origin=_conformer_point(conformer, reactive_atom_id),
            anchor=_conformer_point(conformer, anchor_atom_id),
            metadata={
                "precursor_route": "ortho_hydroxyl_tautomerization",
                "aldehyde_oxygen_atom_id": oxygen_atom_id,
                "ortho_hydroxyl_oxygen_atom_id": hydroxyl_oxygen_atom_id,
                "ortho_hydroxyl_hydrogen_atom_id": hydroxyl_hydrogen_atom_id,
                "ortho_hydroxyl_anchor_atom_id": hydroxyl_anchor_atom_id,
            },
        )

    beta_keto_carbonyl = _find_beta_ketoaldehyde_carbonyl(
        molecule,
        alpha_carbon_atom_id=anchor_atom_id,
        aldehyde_carbon_atom_id=reactive_atom_id,
    )
    if beta_keto_carbonyl is None:
        return None
    carbonyl_carbon_atom_id, carbonyl_oxygen_atom_id = beta_keto_carbonyl
    atom_ids = {
        reactive_atom_id,
        oxygen_atom_id,
        anchor_atom_id,
        carbonyl_carbon_atom_id,
        carbonyl_oxygen_atom_id,
        *_attached_hydrogen_ids(molecule, reactive_atom_id),
        *_attached_hydrogen_ids(molecule, anchor_atom_id),
    }
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(sorted(atom_ids)),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
        metadata={
            "precursor_route": "beta_ketoenol_michael_addition",
            "aldehyde_oxygen_atom_id": oxygen_atom_id,
            "beta_ketoenol_alpha_carbon_atom_id": anchor_atom_id,
            "beta_keto_carbonyl_carbon_atom_id": carbonyl_carbon_atom_id,
            "beta_keto_carbonyl_oxygen_atom_id": carbonyl_oxygen_atom_id,
        },
    )


def _postprocess_keto_aldehyde_matches(
    molecule,
    detected: tuple[_DetectedMotif, ...],
    definition: MotifKindDefinition,
) -> tuple[_DetectedMotif, ...]:
    del definition
    conventional_indices = tuple(
        index
        for index, motif in enumerate(detected)
        if motif.metadata.get("precursor_route") == "ortho_hydroxyl_tautomerization"
    )
    if not conventional_indices:
        return detected
    conventional = tuple(detected[index] for index in conventional_indices)
    resolved = _assign_unique_keto_aldehyde_hydroxyls(molecule, conventional)
    resolved_by_index = dict(zip(conventional_indices, resolved))
    return tuple(resolved_by_index.get(index, motif) for index, motif in enumerate(detected))


def _find_beta_ketoaldehyde_carbonyl(
    molecule,
    *,
    alpha_carbon_atom_id: int,
    aldehyde_carbon_atom_id: int,
) -> tuple[int, int] | None:
    alpha_carbon = molecule.GetAtomWithIdx(alpha_carbon_atom_id)
    if (
        alpha_carbon.GetAtomicNum() != 6
        or alpha_carbon.GetIsAromatic()
        or len(_attached_hydrogen_ids(molecule, alpha_carbon_atom_id)) != 2
    ):
        return None
    options: list[tuple[int, int]] = []
    for carbonyl_carbon in alpha_carbon.GetNeighbors():
        carbonyl_carbon_atom_id = int(carbonyl_carbon.GetIdx())
        if carbonyl_carbon_atom_id == aldehyde_carbon_atom_id:
            continue
        if carbonyl_carbon.GetAtomicNum() != 6:
            continue
        connecting_bond = molecule.GetBondBetweenAtoms(
            alpha_carbon_atom_id,
            carbonyl_carbon_atom_id,
        )
        if connecting_bond is None or abs(float(connecting_bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
            continue
        oxygen_ids = tuple(
            int(bond.GetOtherAtom(carbonyl_carbon).GetIdx())
            for bond in carbonyl_carbon.GetBonds()
            if (
                bond.GetOtherAtom(carbonyl_carbon).GetAtomicNum() == 8
                and abs(float(bond.GetBondTypeAsDouble()) - 2.0) <= 1.0e-6
            )
        )
        carbon_neighbors = tuple(
            neighbor
            for neighbor in carbonyl_carbon.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != alpha_carbon_atom_id
        )
        if len(oxygen_ids) == 1 and carbon_neighbors:
            options.append((carbonyl_carbon_atom_id, oxygen_ids[0]))
    return min(options) if options else None


def _interpret_activated_methylene_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif | None:
    del definition
    reactive_atom_id = int(match[0])
    hydrogen_atom_ids = _attached_hydrogen_ids(molecule, reactive_atom_id)
    if len(hydrogen_atom_ids) < 2:
        return None
    activating_neighbor_ids = tuple(
        neighbor.GetIdx()
        for neighbor in molecule.GetAtomWithIdx(reactive_atom_id).GetNeighbors()
        if neighbor.GetAtomicNum() != 1 and _is_electron_withdrawing_anchor(neighbor, reactive_atom_id)
    )
    if not activating_neighbor_ids:
        return None
    anchor_atom_id = min(activating_neighbor_ids)
    anchor = centroid(tuple(_conformer_point(conformer, atom_id) for atom_id in activating_neighbor_ids))
    atom_ids = sorted({reactive_atom_id, *hydrogen_atom_ids, *activating_neighbor_ids})
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=anchor,
        metadata={
            "hydrogen_atom_ids": hydrogen_atom_ids,
            "activator_atom_ids": activating_neighbor_ids,
        },
    )


def _plane_normal(
    detected: tuple[_DetectedMotif, ...],
    molecule,
    conformer,
) -> _MonomerPlaneFit:
    """Assess the monomer's plane prior from its own geometry.

    With three or more motif anchors the connection points define the plane
    directly (the ring-forming node path). With fewer than three anchors two
    connection points only define a line, so the plane is fitted from the
    conformer's heavy atoms instead — and genuinely non-planar, collinear, or
    degenerate conformers report ``normal=None`` (uncertain) rather than
    borrowing the arbitrary world z axis of the embedding (impact-review
    claim T1-8).
    """
    anchor_points = [item.anchor for item in detected]
    if len(anchor_points) >= 3:
        first, second, third = anchor_points[:3]
        normal = cross(sub(second, first), sub(third, first))
        if not _is_near_zero(normal):
            normalized = normalize(normal)
            if normalized[2] < 0.0:
                normalized = (-normalized[0], -normalized[1], -normalized[2])
            return _MonomerPlaneFit(
                status="anchor-plane",
                normal=normalized,
                canonical_axis=normalized,
                principal_axis=None,
                rms_deviation=None,
            )
        # Collinear anchors define no plane; fall through to the molecular fit.
    heavy_points = tuple(
        _conformer_point(conformer, atom.GetIdx())
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1
    )
    return _fit_molecular_plane(heavy_points)


def _fit_molecular_plane(points: tuple[tuple[float, float, float], ...]) -> _MonomerPlaneFit:
    """Best-fit plane through a heavy-atom point set, with degeneracy handling.

    The normal's sign is canonicalized against the cross product of two
    geometry-derived in-plane reference displacements, so a rigid rotation of
    the input coordinates rotates the reported normal covariantly instead of
    flipping it through a world-axis convention.
    """
    pairs = covariance_eigenpairs(points)
    if pairs is None or pairs.eigenvalues[2] <= 1e-16:
        return _MonomerPlaneFit(
            status="degenerate",
            normal=None,
            canonical_axis=None,
            principal_axis=None,
            rms_deviation=None,
        )
    least_axis, _, principal_axis = pairs.eigenvectors
    rms = max(pairs.eigenvalues[0], 0.0) ** 0.5
    if max(pairs.eigenvalues[1], 0.0) ** 0.5 < _MONOMER_PLANE_COLLINEAR_RATIO * pairs.eigenvalues[2] ** 0.5:
        # One-dimensional point set: any plane containing the line fits
        # equally well, so no normal is claimed. The least-variance axis is
        # still returned for deterministic coordinate canonicalization.
        return _MonomerPlaneFit(
            status="collinear",
            normal=None,
            canonical_axis=least_axis,
            principal_axis=principal_axis,
            rms_deviation=rms,
        )
    normal = _covariant_axis_sign(points, pairs.center, least_axis)
    principal_axis = _covariant_axis_sign(points, pairs.center, principal_axis)
    if rms > _MONOMER_PLANARITY_RMS_TOLERANCE_ANGSTROM:
        return _MonomerPlaneFit(
            status="nonplanar",
            normal=None,
            canonical_axis=normal,
            principal_axis=principal_axis,
            rms_deviation=rms,
        )
    return _MonomerPlaneFit(
        status="planar",
        normal=normal,
        canonical_axis=normal,
        principal_axis=principal_axis,
        rms_deviation=rms,
    )


def _covariant_axis_sign(
    points: tuple[tuple[float, float, float], ...],
    center: tuple[float, float, float],
    axis: tuple[float, float, float],
) -> tuple[float, float, float]:
    """Choose the sign of *axis* from the point set itself.

    The reference is the cross product of the farthest displacement from the
    centroid and the displacement most orthogonal to it — both transform
    covariantly under rigid motions, so the chosen sign does not depend on
    world axes. For (near-)collinear sets no such reference exists and the
    eigensolver's deterministic sign is kept.
    """
    offsets = tuple(sub(point, center) for point in points)
    farthest = max(offsets, key=norm)
    most_orthogonal = max(offsets, key=lambda offset: norm(cross(farthest, offset)))
    reference = cross(farthest, most_orthogonal)
    if norm(reference) < 1e-12:
        return axis
    return axis if dot(axis, reference) >= 0.0 else scale(axis, -1.0)


def _project_onto_plane(
    vector: tuple[float, float, float],
    plane_normal: tuple[float, float, float],
) -> tuple[float, float, float]:
    scale = dot(vector, plane_normal)
    return (
        vector[0] - scale * plane_normal[0],
        vector[1] - scale * plane_normal[1],
        vector[2] - scale * plane_normal[2],
    )


def _is_near_zero(vector: tuple[float, float, float]) -> bool:
    return abs(vector[0]) + abs(vector[1]) + abs(vector[2]) < 1e-8


def _conformer_point(conformer, atom_id: int) -> tuple[float, float, float]:
    position = conformer.GetAtomPosition(atom_id)
    return (float(position.x), float(position.y), float(position.z))


def _attached_hydrogen_ids(molecule, atom_id: int) -> tuple[int, ...]:
    return tuple(
        int(neighbor.GetIdx())
        for neighbor in molecule.GetAtomWithIdx(atom_id).GetNeighbors()
        if neighbor.GetAtomicNum() == 1
    )


def _find_ortho_hydroxyl(molecule, anchor_atom_id: int) -> tuple[int, int, int] | None:
    options = _find_ortho_hydroxyl_options(molecule, anchor_atom_id)
    return options[0] if options else None


def _find_ortho_hydroxyl_options(molecule, anchor_atom_id: int) -> tuple[tuple[int, int, int], ...]:
    anchor_atom = molecule.GetAtomWithIdx(anchor_atom_id)
    options: list[tuple[int, int, int]] = []
    for aromatic_neighbor in anchor_atom.GetNeighbors():
        if aromatic_neighbor.GetAtomicNum() != 6 or not aromatic_neighbor.GetIsAromatic():
            continue
        oxygen_atom_id = _aromatic_hydroxyl_substituent(molecule, aromatic_neighbor.GetIdx(), excluded_atom_id=anchor_atom_id)
        if oxygen_atom_id is None:
            continue
        hydrogen_atom_ids = _attached_hydrogen_ids(molecule, oxygen_atom_id)
        if len(hydrogen_atom_ids) == 1:
            options.append((oxygen_atom_id, hydrogen_atom_ids[0], aromatic_neighbor.GetIdx()))
    return tuple(sorted(set(options)))


def _assign_unique_keto_aldehyde_hydroxyls(molecule, detected: tuple[_DetectedMotif, ...]) -> tuple[_DetectedMotif, ...]:
    if len(detected) <= 1:
        return detected

    options_by_index = [
        _find_ortho_hydroxyl_options(molecule, motif.anchor_atom_id)
        for motif in detected
    ]
    if any(not options for options in options_by_index):
        return detected

    ordered_indices = tuple(sorted(range(len(detected)), key=lambda idx: (len(options_by_index[idx]), detected[idx].reactive_atom_id)))
    best_assignment: dict[int, tuple[int, int, int]] | None = None
    best_score: tuple[int, tuple[int, ...]] | None = None

    def _search(position: int, used_oxygen_ids: set[int], assignment: dict[int, tuple[int, int, int]]) -> None:
        nonlocal best_assignment, best_score
        if position == len(ordered_indices):
            chosen_oxygen_ids = tuple(assignment[idx][0] for idx in range(len(detected)))
            score = (len(set(chosen_oxygen_ids)), tuple(chosen_oxygen_ids))
            if best_score is None or score > best_score:
                best_score = score
                best_assignment = dict(assignment)
            return

        detected_index = ordered_indices[position]
        options = sorted(options_by_index[detected_index], key=lambda option: (option[0] in used_oxygen_ids, option[0], option[2]))
        for option in options:
            assignment[detected_index] = option
            next_used = set(used_oxygen_ids)
            next_used.add(option[0])
            _search(position + 1, next_used, assignment)
            del assignment[detected_index]

    _search(0, set(), {})
    if best_assignment is None:
        return detected

    resolved: list[_DetectedMotif] = []
    for index, motif in enumerate(detected):
        oxygen_atom_id, hydrogen_atom_id, hydroxyl_anchor_atom_id = best_assignment[index]
        original_oxygen_atom_id = motif.metadata.get("ortho_hydroxyl_oxygen_atom_id")
        original_hydrogen_atom_id = motif.metadata.get("ortho_hydroxyl_hydrogen_atom_id")
        atom_ids = set(motif.atom_ids)
        if isinstance(original_oxygen_atom_id, int):
            atom_ids.discard(original_oxygen_atom_id)
        if isinstance(original_hydrogen_atom_id, int):
            atom_ids.discard(original_hydrogen_atom_id)
        atom_ids.update((oxygen_atom_id, hydrogen_atom_id, hydroxyl_anchor_atom_id))
        metadata = dict(motif.metadata)
        metadata["ortho_hydroxyl_oxygen_atom_id"] = oxygen_atom_id
        metadata["ortho_hydroxyl_hydrogen_atom_id"] = hydrogen_atom_id
        metadata["ortho_hydroxyl_anchor_atom_id"] = hydroxyl_anchor_atom_id
        resolved.append(
            _DetectedMotif(
                reactive_atom_id=motif.reactive_atom_id,
                anchor_atom_id=motif.anchor_atom_id,
                atom_ids=tuple(sorted(atom_ids)),
                origin=motif.origin,
                anchor=motif.anchor,
                metadata=metadata,
            )
        )
    return tuple(resolved)


def _is_electron_withdrawing_anchor(atom, excluded_neighbor_id: int) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    if atom.GetIsAromatic():
        aromatic_nitrogen_neighbors = sum(
            1
            for neighbor in atom.GetNeighbors()
            if neighbor.GetIdx() != excluded_neighbor_id and neighbor.GetAtomicNum() == 7 and neighbor.GetIsAromatic()
        )
        if aromatic_nitrogen_neighbors >= 2:
            return True
        if _aromatic_ring_has_conjugated_activation(atom, excluded_neighbor_id):
            return True
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetIdx() == excluded_neighbor_id:
            continue
        bond_order = float(bond.GetBondTypeAsDouble())
        if bond_order >= 2.0 and other.GetAtomicNum() in (7, 8, 16):
            return True
        if bond_order >= 3.0 and other.GetAtomicNum() == 7:
            return True
    return False


def _aromatic_ring_has_conjugated_activation(atom, excluded_neighbor_id: int) -> bool:
    """Recognize methyl donors activated through an aromatic ring.

    Vinylene COFs commonly use methylated aza-heterocycles,
    cyano-substituted aromatic donors, or phenyl methyl groups conjugated to
    an aza-heterocycle.  The withdrawing atom need not be a direct neighbor
    of the methyl-bearing ring carbon, but it must belong to the same aromatic
    ring or to a directly attached five- or six-membered aromatic heterocycle.
    """

    molecule = atom.GetOwningMol()
    Chem.GetSymmSSSR(molecule)
    atom_idx = int(atom.GetIdx())
    for raw_ring in molecule.GetRingInfo().AtomRings():
        ring = {int(candidate) for candidate in raw_ring}
        if atom_idx not in ring or len(ring) not in {5, 6}:
            continue
        if not all(molecule.GetAtomWithIdx(candidate).GetIsAromatic() for candidate in ring):
            continue
        if any(
            molecule.GetAtomWithIdx(candidate).GetAtomicNum() in {7, 8, 16}
            for candidate in ring
        ):
            return True
        for candidate in ring:
            ring_atom = molecule.GetAtomWithIdx(candidate)
            for substituent in ring_atom.GetNeighbors():
                substituent_idx = int(substituent.GetIdx())
                if substituent_idx in ring or substituent_idx == excluded_neighbor_id:
                    continue
                if substituent.GetAtomicNum() != 6:
                    continue
                if any(
                    float(bond.GetBondTypeAsDouble()) >= 3.0
                    and bond.GetOtherAtom(substituent).GetAtomicNum() == 7
                    for bond in substituent.GetBonds()
                ):
                    return True
                if substituent.GetIsAromatic() and _belongs_to_aromatic_heterocycle(
                    molecule,
                    substituent_idx,
                    excluded_ring=ring,
                ):
                    return True
    return False


def _belongs_to_aromatic_heterocycle(
    molecule,
    atom_id: int,
    *,
    excluded_ring: set[int],
) -> bool:
    """Return whether an aromatic atom belongs to an attached heterocycle."""

    for raw_ring in molecule.GetRingInfo().AtomRings():
        ring = {int(candidate) for candidate in raw_ring}
        if ring == excluded_ring or atom_id not in ring or len(ring) not in {5, 6}:
            continue
        if not all(molecule.GetAtomWithIdx(candidate).GetIsAromatic() for candidate in ring):
            continue
        if any(
            molecule.GetAtomWithIdx(candidate).GetAtomicNum() in {7, 8, 16}
            for candidate in ring
        ):
            return True
    return False


def _aromatic_hydroxyl_substituent(molecule, carbon_atom_id: int, *, excluded_atom_id: int) -> int | None:
    for substituent in molecule.GetAtomWithIdx(carbon_atom_id).GetNeighbors():
        if substituent.GetIdx() == excluded_atom_id or substituent.GetAtomicNum() != 8:
            continue
        if len(_attached_hydrogen_ids(molecule, substituent.GetIdx())) == 1:
            return int(substituent.GetIdx())
    return None
