# Plan: fixing the 2D stacking mechanism

Status: revised 2026-06-XX (independently verified against codebase).
Original proposal 2026-05-24. Severity tags map to the audit findings
(C1–C4 = critical, T1–T6 = traps/conventions).

> **Verification note.** Every flaw claim in the original plan was tested
> against the live source. All critical flaws are real and confirmed below.
> Four inaccuracies were found — all path-level or description-level errors,
> none affecting the fixes themselves — and are corrected in-place.

## 1. Goals and non-goals

**Goals**

- No stacking output can contain interpenetrating or sterically impossible
  layers without a loud, recorded diagnostic (C1).
- Named registries (`AA`, `AB`, `slipped`) mean the same geometry in every
  build path and every cell setting (C2).
- Every emitted distance is self-explaining: derivable from the artifact,
  provenance-tagged, and free of silent fallbacks (audit turn 3).
- Ranking is invariant to stacking expansion (C4).
- Validation can actually see interlayer clashes (C3).

**Non-goals (keep current scope honesty)**

- No engine-level stacking exploration (`COFProject.stacking_mode` stays
  `"disabled"`), no stacking-aware scoring/energy terms
  (`stacking_considered: False` stays), no multilayer count > 2, no
  turbostratic/rotated registries, no ASE.
- Seed assembly remains honest about not being physical relaxation. The
  stacking fix produces correct *starting models*; relaxation stays in the
  LAMMPS layer.

## 2. Design constraints (from AGENTS.md)

- Warn over hard stop for non-fatal problems; never abort a batch record.
- Batch-facing data goes into `batch_models.py` dataclasses, not ad hoc dicts.
- Scoring and optimization stay separate.
- Same-change updates: `CHANGELOG.md`, `docs/CURRENT_SCOPE.md`,
  `ARCHITECTURE.md` / `agent-docs/AGENT_CODEBASE_MAP.md` if seams move,
  `skills/cofkit-navigator/SKILL.md` if CLI/output artifacts change.
- Tests mirror `tests/test_*.py`; no LAMMPS requirement in the suite.

## 3. Workstreams

### W1 — Distance semantics and the span pipeline (fixes C1, C2-convention)

**W1.1 Define the quantity we actually compute.**
`interlayer_distance` is used as a *clearance between the extreme nuclear
planes* of adjacent layers (`src/cofkit/stacking.py:93`:
`c2c = interlayer_distance + layer_z_span`), not the literature
mean-plane repeat. Make this explicit:

- Docstring on `LayerRegistry` (`src/cofkit/stacking.py:15-19`) stating the
  clearance semantics, the derivation `c2c = clearance + span`,
  `c = 2·c2c` for the bilayer, and that the resulting interlayer repeat
  (`c/2`) is what to compare to experimental interlayer distances.
- Metadata rename with compatibility alias (see W3 schema):
  `interlayer_clearance` becomes the canonical key;
  `interlayer_distance` remains as a deprecated mirror for one release.
- Registry defaults (`3.4/3.5/3.6`) gain an explicit provenance comment
  ("heuristic nuclear-plane clearances, registry-dependent; not fitted, not
  chemistry-aware") and become overridable per call
  (`LayerRegistry` already a dataclass; add a `clearance=` CLI override
  later — optional, W8).

**W1.2 One shared span measurer.**
There are two independent builders that measure span, plus one reader that
silently falls back to 0.0:

- `src/cofkit/batch.py:1858-1914` — two functions (`_decorated_bex_layer_spacing`
  and `_realized_decorated_bex_layer_spacing`) that measure along the bare
  cartesian z axis. Crucially, **neither** includes `pose.translation` in the
  z accumulation — both only rotate local positions (`matmul_vec(...)[2]`).
- `src/cofkit/build_workflows/ring_forming.py:565-574` — measures along the
  bare cartesian z axis but **does** include `pose.translation` via
  `add(matmul_vec(pose.rotation_matrix, position), pose.translation)[2]`.
- `src/cofkit/stacking.py:277-286` (`_candidate_layer_z_span`) — not an
  independent measurer; only reads whatever was deposited in
  `embedding["layer_z_span"]` and silently returns `0.0` on failure.

The root cause of C1 is that `_candidate_layer_z_span` reads metadata that
may be absent, rather than measuring directly. Introduce a shared measurer,
e.g. in `src/cofkit/geometry.py` or `src/cofkit/stacking.py`:

```python
@dataclass(frozen=True)
class LayerSpanReport:
    span: float
    mode: str          # "atomistic_product" | "precursor_coordinates"
    axis: str          # "c_hat"
    n_atoms: int
    z_translation_included: bool

def measure_layer_z_span(state, monomer_specs, instance_to_monomer,
                         realization=None) -> LayerSpanReport
```

Rules: measure over *realized product atoms when available* (matches the
ring-forming `atomistic_product` behavior), along `normalize(cell[2])`
instead of bare cartesian z, **including `pose.translation`**. Refactor
`_decorated_bex_layer_spacing`, `_realized_decorated_bex_layer_spacing`
(`src/cofkit/batch.py`), and the `_annotate_ring_embedding` helper in
`src/cofkit/build_workflows/ring_forming.py:549-599` to call it.

**W1.3 Remove the silent fallback (the root cause of C1).**
`enumerate_candidate_stackings` gains a `monomer_specs` parameter (both call
sites have specs in hand: `src/cofkit/batch.py:4597-4604` with `first`/`second`,
`src/cofkit/build_workflows/ring_forming.py:102`). `_candidate_layer_z_span`
is replaced by a fresh measurement; the metadata value is used only as a
cross-check with a warning on mismatch. If measurement is genuinely impossible
(coarse pseudo-atom monomers), the span falls back to `precursor_coordinates`
mode; if even that fails, stacking is skipped with
`stacking_skipped:span-unavailable` recorded — never a bare `0.0`.
`max(0.0, span)` clamping (`src/cofkit/stacking.py:286`) becomes a warning.

**W1.4 Unify the "3.4 Å" conventions.**
Four overlapping usages exist today:
- `RingFormationConfig.layer_spacing = 3.4` in
  `src/cofkit/build_workflows/ring_forming.py:39` (single-layer flat repeat)
- `default_ring_layer_spacing = 3.4` in `src/cofkit/engine.py:35`
- stacking clearance default `3.4` in `src/cofkit/stacking.py:35`
- vacuum default `default_layer_spacing = 8.0` in `src/cofkit/engine.py:34`

Decide and document:

- Single-layer exports describe a slab in vacuum (keep c = span + vacuum gap);
  add `c_axis_semantics: "vacuum_slab"` to embedding metadata.
- Stacked exports describe the bulk bilayer (clearance semantics).
- `RingFormationConfig.layer_spacing` becomes `vacuum_gap`-style or is
  documented as "implicit AA repeat"; recommended: switch ring-forming
  single-layer c to `span + 8.0`-style vacuum slab for consistency and note
  the change in CHANGELOG (this alters existing single-layer ring-forming
  outputs — call it out explicitly).

### W2 — Registry correctness and cell-basis normalization (fixes C2)

**W2.1 One shared cell classifier.**
`src/cofkit/batch.py:5003-5014` (`_single_node_cell_kind`) uses cosine
tolerance `abs(cosine - 0.5) < 1e-3` to detect hexagonal cells. A fitted
cell with γ ≈ 60.09° gives cosine ≈ 0.50019, which **fails this test** and
falls through to `"oblique"`, breaking the registry selection. The function
in `src/cofkit/single_node_topologies.py:336` (`_metric_family`) and
`src/cofkit/indexed_topology_layouts.py:171` (`_metric_family`) each have
their own independent classification logic. Move to one function (e.g.
`geometry.classify_2d_cell`):

- angle tolerance ±0.5°, relative length tolerance 1% (fitted cells measure
  γ = 60.09° today and fall through to `"oblique"`).
- returns `(kind, setting)` where `setting ∈ {"60deg", "120deg"}` for
  hexagonal (the registry shift depends on it).

Update `src/cofkit/embedding.py` (`_cell_kind`, `_single_node_bipartite_cell_kind`),
`src/cofkit/batch.py` (`_topology_cell_kind`), and stacking to use it.

**W2.2 Basis-aware registry shifts.**
The honeycomb Bernal shift (vertex → ring/pore center, length `a/√3`) is
`(1/3, 2/3)` in the 120° setting and `(1/3, 1/3)` in the 60° setting. The
registry table (`src/cofkit/stacking.py:32-43`) stores only the 60° form
(`(1/3, 1/3)`) while the ring-forming/indexed builders emit 120° cells
(measured: applied shift |s| = a/3 = 5.04 Å instead of a/√3 = 8.73 Å for
|a| = 15.1 Å). Fix:

- Store per-setting shifts (`AB_hex_60 = (1/3, 1/3)`,
  `AB_hex_120 = (1/3, 2/3)`), selected by `setting`.
- `slipped = (1/2, 0)` and `AB_square = (1/2, 1/2)` are setting-independent;
  document their meaning (half-vector slip / hollow-site registry).
- Optional (W8): post-shift symmetry check that the shifted layer's net
  vertices eclipse expected high-symmetry sites; warn otherwise.

**W2.3 Update the tests and docs that enshrine the wrong value.**
`tests/test_ring_forming.py:205-215` asserts `(1/3, 1/3)` unconditionally
(the assertion is at line 214); `docs/building.md:106` documents it.
Replace with setting-aware assertions (see W7).

### W3 — Self-explaining artifacts (audit turn 3)

**W3.1 Metadata schema (extend, keep aliases).**

```python
metadata["stacking"] = {
    "id": "AA",
    "layer_count": 2,
    "interlayer_clearance": 3.4,              # NEW canonical (alias: interlayer_distance)
    "registry_shift_fractional": (0.0, 0.0),  # NEW canonical (alias: lateral_shift_fractional)
    "cell_classification": {"kind": "hexagonal", "setting": "60deg", "gamma_deg": 60.09},
    "layer_z_span": 3.586,
    "layer_z_span_mode": "atomistic_product", # NEW (ring-forming already computes this; stop dropping it)
    "layer_z_span_axis": "c_hat",             # NEW
    "center_to_center_distance": 6.986,
    "derivation": "c2c = interlayer_clearance + layer_z_span; c = layer_count * c2c",
    "min_interlayer_contact": 4.29,           # NEW (self-check result)
    "comment_suffix": "stacking=AA",
    "source_candidate_id": "...",
    "warnings": [],                           # NEW
}
```

Route this through a `StackingSummary` dataclass in
`src/cofkit/batch_models.py` (per AGENTS.md) and render it in `summary.md` /
console (`src/cofkit/batch.py:4711-4713`, `src/cofkit/cli_build.py:641-665`).

**W3.2 Derivation block in the CIF.**
Next to the existing `# COFid: ... stacking=AA` line (`src/cofkit/cofid.py`,
`src/cofkit/cif.py:_render_cif`), emit:

```
# stacking-geometry: registry=AA basis=hexagonal:60deg shift_frac=0,0 \
#   interlayer_clearance=3.4 layer_z_span=3.586 (atomistic_product, c_hat) \
#   c2c=6.986 c=2*c2c min_interlayer_contact=4.29
```

so a CIF alone answers "why is c this long?" (`_STACKING_COMMENT_PREFIX`
keeps `read_cofid_comment_from_cif` parsing intact — verify
`src/cofkit/lammps.py`/`src/cofkit/graspa.py` comment round-tripping still
works; those parse the first line only, and the stacking block goes below it).

### W4 — Validation and safety rails (fixes C1 residual, C3)

**W4.1 A clash criterion that can see clashes.**
`CoarseValidationThresholds.min_nonbonded_heavy_distance = 1.05`
(`src/cofkit/validation.py:34`) is used as **both** the acceptance threshold
and the `cutoff` argument to `NeighborSearch`
(`src/cofkit/validation.py:426-427`). With a 1.05 Å search radius the search
never finds a 2 Å interpenetration contact. Change to a vdW-aware ratio: flag
when `d / (r_vdw(i) + r_vdw(j)) < 0.75` for nonbonded heavy pairs, with an
independent search radius (≈ 3.5 Å) and H contacts reported separately
(`min_nonbonded_h_distance`). Keep the plain-distance threshold as a coarse
backstop (e.g. 2.2 Å heavy-heavy). Diagnostic only at first (warn +
`diagnostics` entry) per the warn-over-stop principle; promote to
classification failure one release later.

**W4.2 Stacking self-check at construction.**
In `_apply_layer_registry`, after building the bilayer: verify
`c2c > layer_z_span` (warn otherwise) and compute `min_interlayer_contact`
over the two layers with minimum-image distances (cheap; monomer-level then
atom-level refinement). Sub-vdW contacts produce
`warning:` on stderr plus `stacking.warnings` entries and a
`stacking_clash` flag. Acceptance target: shipped evidence structures
(`out/stacking_audit_evidence/SHIPPED_binary-bridge_*.cif`) must be flagged
or rejected by this check.

**W4.3 No more silent eligibility skips.**
`_is_eligible_2d_candidate` swallows all exceptions
(`src/cofkit/stacking.py:258-265`); `enumerate_candidate_stackings` silently
returns the input. When `registry_ids` were requested but stacking is skipped,
attach `stacking_skipped:<reason>` to `candidate.flags`, record the reason in
metadata, and warn on stderr (batch: per-record diagnostic, no abort).

Both call sites of `enumerate_candidate_stackings` are:
- `src/cofkit/batch.py:4597-4604`
- `src/cofkit/build_workflows/ring_forming.py:102`

The `monomer_specs` parameter added in W1.3 should be keyword-only; `None`
triggers the W4.3 skip-with-warning path.

**W4.4 Axis guard.** Warn if `cell[2]` is not orthogonal to the ab plane
(the offsets ride along `normalize(cell[2])` while spans were measured on
cartesian z — W1.2 unifies this, but the geometry only makes sense for
c ⊥ ab).

### W5 — Ranking integrity (fixes C4)

**W5.1 Residual bookkeeping.**
`_duplicate_ring_geometry_metrics` doubles `total_residual`
(`src/cofkit/stacking.py:322-327`) but the duplicated `bridge_event_metrics`
leave `bridge_geometry_residual` (`src/cofkit/scoring.py:81`) undoubled, so
`residual_ranking_key` (`src/cofkit/model.py:157-164`) computes a halved
per-event mean for stacked candidates (measured: 0.0029 vs true 0.00577),
because `n_bridge_events` doubles while the numerator stays the same. Fix by
doubling `bridge_geometry_residual` alongside the metrics (mirror the ring
branch), and add a parity test: the same candidate pre/post stacking must
yield the same `candidate_ranking_key` mean residual.

**W5.2 Keep enumeration order honest.**
Registry variants tie on residual and sort lexicographically
(`…__AA < …__AB < …__slipped`). Document in `docs/building.md` and the CLI
help that registry order is enumerative, not energetic. Optionally sort
registry variants by `registry_id` explicitly and stop letting them occupy
"topology" rank slots: compute `topology_rank`/`topology_count` on the
pre-expansion candidate list (`src/cofkit/batch.py:4607-4639`), and add
`stacking_variant_index` instead.

### W6 — Latent trap hardening (T1–T6)

1. **Ring-center reinterpretation (T1).** `ring_center_fractional` is reused
   against the *new* cell with `ring_center_offset_cartesian` compensation
   (`src/cofkit/reaction_realization.py:3026-3044`,
   `src/cofkit/ring_geometry.py:105-115`) — exact only when its c-component
   `w = 0`. On stacking: assert `|w| < 1e-9` (clear error naming the event)
   or, better, rewrite `ring_center_fractional` into the new cell basis from
   a stored `ring_center_cartesian`. Also make `ring_geometry.py` consume the
   offset metadata through one shared helper (it already calls
   `_event_ring_center_offset`; keep the coupling documented).
2. **Torsion keys (T2).** `state.torsions` is copied without L0/L1 key
   suffixing (`src/cofkit/stacking.py:239`). Suffix keys per layer (or
   document as passthrough-only until a consumer exists).
3. **Span/axis consistency (T3, T4).** Covered by W1.2 (measure along
   `c_hat`, include pose translations) plus the W4.4 guard.
4. **`_safe_normalize` fallback (T5).** The `(0,0,1)` fallback for a
   degenerate c (`src/cofkit/stacking.py:341-344`) should warn.
5. **Cosine tolerance (T6).** Covered by W2.1.
6. **Eligibility by `net_plan` presence (T6-adjacent).** Seed-assembly
   candidates without `net_plan` skip stacking silently — covered by W4.3.

### W7 — Tests, docs, changelog

New `tests/test_stacking.py` (plus targeted additions to existing files,
mirroring current layout):

| # | Test | Guards |
|---|------|--------|
| 1 | TAPB/TFB hcb `--stacking AA`: `c == 2*(clearance + span)` with span > 3 Å, min interlayer contact > 2.5 Å | C1 regression |
| 2 | Candidate without `embedding.layer_z_span` gets a *measured* span; never `0.0` | C1 root cause |
| 3 | Hexagonal AB shift maps a vertex onto the ring/pore center in both 60° and 120° cells (`|s_cart| == a/√3`) | C2 |
| 4 | `classify_2d_cell` on fitted γ=60.09 cell returns `("hexagonal", "60deg")` | C2/T6 |
| 5 | Ranking parity: `candidate_ranking_key` invariant under stacking expansion | C4 |
| 6 | Validator flags a 2.0 Å heavy-heavy nonbonded pair (fixture: trimmed `SHIPPED_binary-bridge_AB.cif`) | C3 |
| 7 | Stacking self-check flags the original shipped interpenetrated geometry | C1/C3 |
| 8 | Ring event with `ring_center_fractional` w≠0 raises/warns (T1) | T1 |
| 9 | Metadata + CIF carry the derivation block; `read_cofid_comment_from_cif` and lammps/graspa comment round-trips still pass | W3 |
| 10 | `stacking_skipped:<reason>` + stderr warning when eligibility fails with `registry_ids` set | W4.3 |
| 11 | Existing suites updated: `test_ring_forming` AB assertions become setting-aware (line 214); `test_batch`/`test_cif` gain the new metadata keys | W2.3/W3 |

Docs in the same change: `docs/building.md` (semantics table: clearance vs
repeat, registry definitions per setting, "registry order is not energetic"),
`docs/python-api.md` (new `monomer_specs` parameter),
`docs/CURRENT_SCOPE.md`, `agent-docs/AGENT_CODEBASE_MAP.md` (new
measurer/classifier seams in `src/cofkit/build_workflows/ring_forming.py`
and `src/cofkit/geometry.py`),
`skills/cofkit-navigator/SKILL.md` (new metadata keys, CIF comment line,
warning behavior), `CHANGELOG.md` (noting the two observable changes:
ring-forming single-layer c convention if W1.4 is accepted, and richer
stacking metadata/comments).

### W8 — Follow-ups (tracked, explicitly out of this change)

- `src/cofkit/lammps.py` data prep: unwrap molecules via the spanning tree
  instead of keeping wrapped CIF coordinates (kills the "Inconsistent image
  flags" / "bond extent > half box" warnings seen in the optimization runs);
  offer a layered-structure preset (tighter box-relax, damped c response).
- Rerun the base/fixed/shipped comparison converged; adjudicate the gallery
  pairing (0.89/7.93 Å) as physics vs artifact.
- Optional `--stacking-clearance` CLI override; chemistry-aware clearance
  defaults (element-pair vdW based) instead of the 3.4/3.5/3.6 constants.
- Optional rigid-layer registry/clearance scan ("stacking_mode: scan")
  through the existing optimizer — only if/when engine scope is revisited.

## 4. Sequencing

| Milestone | Contents | Risk if shipped alone |
|-----------|----------|-----------------------|
| **M1 — stop the bleeding** | W1.2, W1.3, W4.2, W5.1 + tests 1, 2, 5, 7 | Removes impossible outputs and the ranking bias; metadata names still ambiguous |
| **M2 — registry truth** | W2.1, W2.2, W2.3 + tests 3, 4 | AB geometry changes (intentional); downstream comparisons against old outputs shift |
| **M3 — self-explanation + validator** | W3.1, W3.2, W4.1 + tests 6, 9 | Validator becomes stricter (warn-level first) |
| **M4 — traps + docs** | W1.1, W1.4, W4.3, W4.4, W5.2, W6 + tests 8, 10, 11 + full doc pass | — |

M1 is the only *urgent* milestone; M2 changes observed AB geometry and must
be called out in CHANGELOG as a correction, not a regression.

## 5. Compatibility and risk register

| Change | Blast radius | Mitigation |
|--------|--------------|------------|
| `enumerate_candidate_stackings(..., monomer_specs=...)` | Python API (`src/cofkit/batch.py`, `src/cofkit/build_workflows/ring_forming.py`, external users) | Keyword arg; `None` triggers the W4.3 skip-with-warning path; CHANGELOG migration note |
| `bridge_geometry_residual` doubled | Any consumer reading the raw total | Documented as per-layer duplication (consistent with `ring_geometry.total_residual`); ranking unaffected within a layer set |
| Metadata key renames | `summary.md`/JSONL consumers, skills docs | Old keys kept as aliases one release; `StackingSummary.to_dict()` emits both |
| Ring-forming single-layer c convention (W1.4) | Existing single-layer ring-forming CIFs change c | Gate behind a config flag defaulting to the new behavior; CHANGELOG entry; test updates |
| Stricter validation (W4.1) | Batch records currently "valid" may flip to warn | Warn-only first release (per AGENTS warn-over-stop), promote later |
| AB geometry correction (M2) | Any golden outputs / user comparisons | Treat as bug fix; regenerate goldens; note in CHANGELOG |
| Existing tests asserting `(1/3, 1/3)` | `tests/test_ring_forming.py:214` | Update in the same change (W2.3) |

## 6. Definition of done

1. Regenerated evidence pack (the TAPB/TFB AA/AB/slipped trio) shows min
   interlayer contact ≥ 2.5 Å with `c = 2·(clearance + span)`, and the
   original shipped interpenetrated CIFs are flagged by `cofkit validate`.
2. No code path can emit `layer_z_span: 0.0` without a measurement behind it
   (grep gate in tests).
3. Hexagonal AB shift passes the vertex→pore-center assertion in both
   settings.
4. `candidate_ranking_key` parity test green.
5. CIFs carry the `# stacking-geometry:` derivation line; metadata carries
   the full W3.1 schema.
6. `uv run pytest -q` green; `uv lock --check`, ruff gate, `uv build --wheel`
   green; all doc/skill/changelog updates landed in the same change.
