# Stacking fix plan — implementation status

Snapshot against `STACKING_FIX_PLAN.md` as of 2026-06.

Legend: ✅ done · ⚠️ partial · ❌ not started

---

## W1 — Distance semantics and span pipeline

| Item | Status | Notes |
|---|---|---|
| W1.1 `LayerRegistry` clearance semantics docstring | ✅ | `stacking.py:17-38`; `interlayer_clearance` canonical key (deprecated `interlayer_distance` metadata alias removed 2026-09-25, pre-release) |
| W1.2 `measure_layer_z_span` shared measurer in `geometry.py` | ✅ | `LayerSpanReport` dataclass + `measure_layer_z_span()`; explicit `axis`/`axis_label` (2026-09-25: axis is now the layer normal via `layer_normal_axis` = normalize(a × b), not c_hat), includes `pose.translation` |
| W1.2 `_decorated_bex_layer_spacing` in `batch.py` uses shared measurer | ✅ | Done 2026-09-25 (audit A7): delegates to `measure_layer_z_span` along the fitted layer normal; precursor estimate has zero translations by construction (poses do not exist yet) |
| W1.2 `_realized_decorated_bex_layer_spacing` in `batch.py` uses shared measurer | ✅ | Done 2026-09-25 (audit A7): delegates to `measure_layer_z_span` with the realized product atoms; `pose.translation` now included — the C1 root cause is fixed |
| W1.2 `_annotate_ring_embedding` in `ring_forming.py` uses shared measurer | ✅ | Done 2026-09-25 (audit A7): inline loop deleted; delegates to `measure_layer_z_span` along the fitted cell's layer normal and reports `layer_z_span_mode`/`layer_z_span_axis` from the report |
| W1.3 `enumerate_candidate_stackings` `monomer_specs` kwarg | ✅ | Keyword-only; both call sites updated: `batch.py:4600` passes `{first.id: first, second.id: second}`; `ring_forming.py:104` passes `{monomer.id: monomer}` |
| W1.3 Silent `0.0` fallback replaced by skip-with-warning | ✅ | `stacking.py:503-550` (`_measure_span`): fresh measurement; falls back to `precursor_coordinates` mode; if unavailable records `stacking_skipped:span-unavailable` |
| W1.4 `c_axis_semantics: "vacuum_slab"` for single-layer exports | ❌ | Not started; `RingFormationConfig.layer_spacing` and `default_ring_layer_spacing` still undocumented; no `c_axis_semantics` key emitted |

---

## W2 — Registry correctness and cell-basis normalisation

| Item | Status | Notes |
|---|---|---|
| W2.1 `classify_2d_cell` in `geometry.py` | ✅ | ±0.5° angle tolerance, 1% relative length tolerance; returns `(kind, setting)` with `"60deg"` / `"120deg"` disambiguation for hexagonal |
| W2.1 `batch.py` `_single_node_cell_kind` migrated | ✅ | `batch.py:5003` now calls `classify_2d_cell`; old cosine-tolerance logic removed |
| W2.1 `embedding.py` `_cell_kind_from_vectors` helper added | ✅ | Now the single classification route in embedding (2026-09-25, audit A8): `_cell_kind` (topology-id-derived) and `_single_node_bipartite_cell_kind` (60°-only cosine) deleted; both metadata sites classify the built cell vectors, warning + `"oblique"` on failure |
| W2.1 `single_node_topologies.py` `_metric_family` migrated | ✅ | Done 2026-09-25 (audit A8): delegates to the new `geometry.classify_2d_cell_parameters` adapter |
| W2.1 `indexed_topology_layouts.py` `_metric_family` migrated | ✅ | Done 2026-09-25 (audit A8): 2D branch delegates to `classify_2d_cell_parameters`; the 3D branch (cubic/orthorhombic/triclinic) is unchanged |
| W2.2 Hexagonal AB shift cell-setting-aware | ✅ | `StackingExplorer` takes `cell_setting` param; 60° → `(1/3, 1/3)`, 120° → `(1/3, 2/3)` |
| W2.3 `tests/test_ring_forming.py:214` AB assertion updated | ❌ | Still asserts `(1/3, 1/3)` unconditionally; needs to be setting-aware |

---

## W3 — Metadata schema and CIF annotation

| Item | Status | Notes |
|---|---|---|
| W3.1 `StackingSummary` dataclass in `batch_models.py` | ❌ dropped | Removed 2026-09-25 (audit D2/#23): never routed — `batch.py` reads the raw metadata dict |
| W3.1 Full metadata schema in `_stacking_metadata()` | ✅ | `stacking.py` `_stacking_metadata`; canonical keys only (deprecated aliases removed 2026-09-25, pre-release); `cell_classification`, `layer_z_span_mode`, `layer_z_span_axis`, `derivation`, `warnings` all emitted |
| W3.2 CIF `# stacking-geometry:` derivation comment | ✅ | `cif.py:548-596` `_stacking_geometry_comment()` reads from `metadata["stacking"]`; emits single-line comment with registry, basis, shift, clearance, span, c2c, contact |
| W3.1 `StackingSummary` routed through `batch.py` summary/console output | ❌ | Obsolete — dataclass dropped (audit D2); `batch.py` keeps reading the raw metadata dict |

---

## W4 — Validation and safety rails

| Item | Status | Notes |
|---|---|---|
| W4.1 Wider search radius for clash detection | ✅ | `nonbonded_heavy_search_radius: float = 3.5` in `CoarseValidationThresholds`; sized to cover the largest flaggable vdW sum in the Bondi table (Si...Si 3.15 Å) plus margin |
| W4.1 vdW-ratio criterion `d / (r_i + r_j) < 0.75` | ✅ | Done 2026-09-25 (audit #11): wired in `validation.py:_nonbonded_contact_scan` over the shared Bondi table `cofkit.vdw` (not DREIDING radii) with periodic-image-aware 1-2/1-3/1-4 exclusions (`_bond_graph_exclusions`), an explicit 1-4 policy (excluded from the ratio check, heavy-heavy pairs still floored at 2.2 Å), a severe-overlap check for excluded pairs (`excluded_pair_severe_overlap`), and a separate hydrogen channel (`min_nonbonded_hydrogen_*` metrics + `hydrogen_atom_clash` warning) |
| W4.2 `c2c > layer_z_span` self-check | ✅ | `stacking.py:204-213` warns + records when `center_to_center_distance <= layer_z_span` |
| W4.2 `min_interlayer_contact` computed and flagged | ✅ replaced | 2026-09-25 (audit #12): the monomer-center estimate (could never fire) is replaced by `_measure_interlayer_contact` — atomistic contact over the realized layer atoms (sharing the span path's `ReactionRealizer` realization; honest `atomistic_product`/`precursor_coordinates`/`mixed` mode) across both `(0,0,±1)` c galleries via `vdw.min_periodic_pair_contact`; `stacking_clash` uses the shared vdW-ratio + 2.2 Å-floor criterion; metadata records atoms/image/vdW ratio/H involvement/cutoff; the CIF `# stacking-geometry:` comment labels the value as an atomistic contact; measurement failures warn + degrade |
| W4.3 `stacking_skipped:<reason>` flag + stderr warning | ✅ | `stacking.py:128-143`; eligibility failures attach flag to candidate and warn on stderr |
| W4.4 c-axis orthogonality guard | ✅ | `_c_axis_is_orthogonal()` at `stacking.py:660-673`; warns when c is not ⊥ to ab within 2° |

---

## W5 — Ranking integrity

| Item | Status | Notes |
|---|---|---|
| W5.1 `bridge_geometry_residual` doubled on stacking | ✅ replaced | 2026-09-25 (audit A9): the doubling hack and its `try/except: pass` are gone; the sum-over-events contract is documented on `scoring.BridgeGeometryReport` and `stacking.py` recomputes each aggregate from the duplicated per-event metrics in one place (`_recompute_total_residual`, warns on non-numeric data) |
| W5.1 Ranking parity test | ✅ | `tests/test_stacking.py` (2026-09-25): normalized residual + unreacted-motif tie-breaker invariant under stacking expansion (bridge and ring-geometry paths), expanded-vs-unexpanded ordering falls to the id tie-breaker — per the review, not the whole ranking-key tuple |
| W5.2 Registry sort order documented | ❌ | `topology_rank`/`stacking_variant_index` separation not done; not documented in `docs/building.md` |

---

## W6 — Latent trap hardening

| Item | Status | Notes |
|---|---|---|
| T1 `ring_center_fractional` w≠0 warning | ✅ | `stacking.py:269-284`; warns + records in `stacking_warnings` |
| T2 Torsion key passthrough documented | ✅ | `stacking.py:302-303`; comment documents passthrough-only, no per-layer suffixing yet |
| T3/T4 Span/axis consistency | ✅ | Covered by W1.2 — `measure_layer_z_span` measures along the layer normal (`layer_normal_axis`), not c_hat or bare z (updated 2026-09-25, audit A7) |
| T5 Degenerate c-axis fallback warns | ✅ | `_safe_normalize` at `stacking.py:686-698` |
| T6 Fitted γ = 60.09° cell classified correctly | ✅ | `classify_2d_cell` ±0.5° tolerance handles it |

---

## W7 — Tests, docs, changelog

### Tests

| # | Description | Status |
|---|---|---|
| 1 | TAPB/TFB hcb AA: `c == 2*(clearance+span)`, span > 3 Å, min contact > 2.5 Å | ⚠️ | Span-pipeline half covered 2026-09-25 by `tests/test_geometry.py` (translation-inclusive spans, tilted-cell layer-normal invariance) and the decorated-bex tests in `tests/test_batch.py`; the full TAPB/TFB AA integration assertion is still unwritten |
| 2 | Span-never-zero gate: candidate without `embedding.layer_z_span` gets measured span | ⚠️ | 2026-09-25: `tests/test_geometry.py` distinguishes a measured 0.0 (planar layer, `n_atoms > 0`) from `mode="unavailable"`; the candidate-level gate test is still unwritten |
| 3 | Hexagonal AB shift: vertex→pore-center in both 60° and 120° cells | ❌ |
| 4 | `classify_2d_cell` on fitted γ=60.09° cell returns `("hexagonal", "60deg")` | ✅ | `tests/test_geometry.py::CellClassificationTests` (2026-09-25), plus a cross-module consistency test through both migrated `_metric_family` paths |
| 5 | Ranking parity: `candidate_ranking_key` invariant under stacking expansion | ✅ | `tests/test_stacking.py` (2026-09-25, audit A9) — normalized residual + tie-breaker invariance and expanded/unexpanded ordering; the id component intentionally changes |
| 6 | Validator flags 2.0 Å heavy-heavy nonbonded pair | ✅ | `tests/test_validation.py::test_validator_flags_two_angstrom_nonbonded_heavy_pair` (2026-09-25, audit #11), plus benzene / angle-chain / cis-1-4 negative controls, severe-overlap tests, and hydrogen-channel tests |
| 7 | Stacking self-check flags original shipped interpenetrated geometry | ✅ | `tests/test_stacking.py::test_stacking_self_check_flags_interpenetrated_bilayer` (synthetic sub-vdW bilayer, 2026-09-25, audit #12) plus a cross-c-boundary gallery regression at the shared-contact level |
| 8 | `ring_center_fractional` w≠0 warns (T1) | ❌ |
| 9 | Metadata + CIF carry derivation block; lammps/graspa round-trips pass | ❌ |
| 10 | `stacking_skipped:<reason>` + stderr warning when eligibility fails | ❌ |
| 11 | `test_ring_forming.py:214` AB assertion made setting-aware | ✅ | Asserts the 120°-setting `(1/3, 2/3)` shift (W2.3) |
| — | `resolve_stacking_pattern` and `compute_interlayer_offset` unit tests | ❌ dropped | Dead API deleted with its tests 2026-09-25 (audit D1/#21) |

### Docs and changelog

| File | Status |
|---|---|
| `CHANGELOG.md` | ❌ |
| `docs/CURRENT_SCOPE.md` | ❌ |
| `ARCHITECTURE.md` | ❌ |
| `agent-docs/AGENT_CODEBASE_MAP.md` | ❌ |
| `skills/cofkit-navigator/SKILL.md` | ❌ |
| `docs/building.md` (semantics table, registry-order note) | ❌ |
| `docs/python-api.md` (`monomer_specs` parameter) | ❌ |

---

## Definition-of-done checklist

From `STACKING_FIX_PLAN.md` §6:

| Criterion | Status |
|---|---|
| 1. Regenerated evidence pack shows min contact ≥ 2.5 Å; shipped interpenetrated CIFs flagged by `cofkit validate` | ❌ — self-check infra is in place but evidence pack not regenerated |
| 2. No code path can emit `layer_z_span: 0.0` without a measurement (grep gate in tests) | ❌ — `_decorated_bex_layer_spacing` can still return 0.0 without measurement |
| 3. Hexagonal AB shift passes vertex→pore-center assertion in both 60° and 120° settings | ❌ — test not written |
| 4. `candidate_ranking_key` parity test green | ✅ — `tests/test_stacking.py` (2026-09-25, audit A9) |
| 5. CIFs carry `# stacking-geometry:` derivation line; metadata carries full W3.1 schema | ✅ |
| 6. `uv run pytest -q` green; ruff gate; `uv build --wheel` green; all doc/skill/changelog updates landed | ⚠️ — existing tests pass; new integration tests not yet written; docs/changelog not updated |

---

## Highest-priority remaining work

1. ~~**W1.2 `batch.py` span measurers**~~ — **done 2026-09-25** (audit A7):
   both helpers delegate to `measure_layer_z_span` along the layer normal and
   include `pose.translation`; the C1 root cause is closed.

2. ~~**W4.1 vdW-ratio loop**~~ — **done 2026-09-25** (audit #11): the ratio
   comparison is wired into `validation.py:_nonbonded_contact_scan` over the
   shared `cofkit.vdw` Bondi table with graph-aware exclusions; W4.2's dead
   monomer-center check is likewise replaced by the atomistic
   `_measure_interlayer_contact` (audit #12).

3. **W7 integration tests** — tests 1–11 from the plan are all absent. Tests 1, 2, 5, and 7 directly guard the four critical flaws.

4. ~~**W2.1 remaining classifiers**~~ — **done 2026-09-25** (audit A8):
   `single_node_topologies.py:336` and `indexed_topology_layouts.py:171` both
   delegate to `classify_2d_cell` via `classify_2d_cell_parameters`.

5. **Docs and changelog** — all seven doc files listed above are unmodified.

---

## Post-audit update — 2026-09-25 (commits cf9776d, 408fbec, 6c966e8)

The builder edge-case audit (`agent-docs/BUILDER_EDGE_CASE_AUDIT.md`) re-bucketed
this plan's remaining items. Current accounting:

**Implemented and kept:** W1.1 docstring; W1.2 shared measurer (call sites
migrated 2026-09-25, audit A7 — see below); W1.3 `monomer_specs` + skip-with-warning (caveat: the
`stacking_skipped:span-unavailable` flag promised in the `stacking.py:123`
docstring is never attached — only eligibility skips tag flags); W2.1 canonical
`classify_2d_cell` (all copies migrated 2026-09-25, audit A8); W2.2 setting-aware AB shift;
**W2.3 now ✅** (test asserts the 120°-setting `(1/3, 2/3)` shift);
W3.1 canonical metadata schema; W3.2 CIF derivation comment; W4.1 wider search
radius only; W4.2/W4.3/W4.4 warn-only shells; W5.1 residual handling (the
doubling symptom patch was replaced by the sum-over-events contract fix on
2026-09-25, audit A9); W6 T1/T2/T5/T6.
Hardening beyond the plan: per-candidate exception isolation in batch stacking
expansion and narrowed exception swallowing with stderr warnings (408fbec).

**Strategically cancelled** (Tier 3 drops, 6c966e8): `StackingSummary` and its
routing (W3.1's ❌ closed by deletion — the raw dict stays the channel); the
"kept for one release" deprecated aliases (never shipped — removed outright);
`tests/test_stacking.py` (the only ✅ test row above — it covered the deleted
dead API); the `stacking_considered` / `stacking_penalty` stamps (removed with
the legacy scorer).

**Awaiting implementation:** ~~W1.2 root cause~~ **landed 2026-09-25**
(audit A7: both `batch.py` span helpers and `ring_forming._annotate_ring_embedding`
now delegate to `measure_layer_z_span` along the layer normal; `pose.translation`
included; honest `layer_z_span_axis`/`layer_z_span_mode` provenance);
~~W2.1's two remaining classifiers~~ **landed 2026-09-25** (audit A8: both
`_metric_family` copies delegate to `classify_2d_cell` via the new
`classify_2d_cell_parameters` adapter; embedding and ring-forming classify the
built/fitted cell, keeping topology-declared families as provenance only);
audit A10 primitive consolidation **landed 2026-09-25** (`safe_normalize` /
`orthogonal_component` / `angle_degrees` / `distance` in `geometry.py`;
`scoring.py`'s `_safe_normalize` copy landed with A9, delegating to
`geometry.safe_normalize` with its `(1,0,0)` fallback preserved);
audit A9 **landed 2026-09-25** (sum-over-events residual contract on
`scoring.BridgeGeometryReport`; stacking recomputes aggregates from duplicated
per-event metrics via `_recompute_total_residual` instead of the doubling
hack; ranking parity tests in `tests/test_stacking.py`); the finding-#16
soft-relax keep-conditions **landed 2026-09-25** (non-converged passes keep
the original structure; pair-list decisions evaluated at the minimum-distance
image; clash cutoff shared with `CoarseValidationThresholds`).
Still awaiting: W1.4 `c_axis_semantics`; W5.2; and W7's remaining integration tests plus user-facing docs.
~~W4.1 vdW wiring (audit #11)~~ **landed 2026-09-25** (ratio criterion + 2.2 Å
floor wired in `validation.py:_nonbonded_contact_scan` over the new shared
`cofkit.vdw` Bondi table with periodic-image-aware 1-2/1-3/1-4 exclusions, an
excluded-pair severe-overlap check, and a separate hydrogen metric channel;
negative controls for ideal benzene, ordinary angle neighbors, and cis 1-4
contacts in `tests/test_validation.py`) and ~~a real atomistic contact metric
to replace the mathematically dead W4.2 check (audit #12)~~ **landed
2026-09-25** (`stacking._measure_interlayer_contact` measures realized layer
atoms — sharing the span realization — across both `(0,0,±1)` c galleries via
`cofkit.vdw.min_periodic_pair_contact`; `stacking_clash` aligned to the shared
criterion; atom labels/image/H-involvement/cutoff recorded; CIF comment labels
the value atomistic; per-candidate failure isolation).
Definition-of-done criteria 1, 3 and 6 remain unmet; criterion 2's grep gate
is now covered behaviorally by the `tests/test_geometry.py` span tests.

**Second post-audit update — 2026-09-25 (later same day):** audit #13 **landed**
(owner chose orthogonalize-on-stacking): `_apply_layer_registry` builds the
stacked c along the layer normal with length `2·c2c`, so the W4.4 warn-only
shell is superseded — the derivation is true by construction, tilt is recorded
as provenance (`c_axis_basis` / `c_axis_orthogonalized` / `base_c_tilt_degrees`)
and in the CIF comment, and the degenerate in-plane case falls back with a
warning. Audit #14 **landed** (quinoid ring re-bond-ordering in keto-enamine
realization via constrained perfect matching). Audit #15 consistency fix
**landed** (per-template retraction in mixed-linkage builds); full calibration
remains an explicitly documented future project. These close the remaining
c-axis correctness half of review point 2; W1.4 `c_axis_semantics` (the
vacuum-slab vs implicit-AA-repeat distinction for single-layer exports) is
still open and is now the only remaining c-semantics item.
