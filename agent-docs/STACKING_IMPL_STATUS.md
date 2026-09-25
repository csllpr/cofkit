# Stacking fix plan — implementation status

Snapshot against `STACKING_FIX_PLAN.md` as of 2026-06.

Legend: ✅ done · ⚠️ partial · ❌ not started

---

## W1 — Distance semantics and span pipeline

| Item | Status | Notes |
|---|---|---|
| W1.1 `LayerRegistry` clearance semantics docstring | ✅ | `stacking.py:17-38`; `interlayer_clearance` canonical key (deprecated `interlayer_distance` metadata alias removed 2026-09-25, pre-release) |
| W1.2 `measure_layer_z_span` shared measurer in `geometry.py` | ✅ | `LayerSpanReport` dataclass + `measure_layer_z_span()` added; projects onto `c_hat`, includes `pose.translation` |
| W1.2 `_decorated_bex_layer_spacing` in `batch.py` uses shared measurer | ❌ | `batch.py:1869` still does `matmul_vec(rotation, position)[2]` with no `pose.translation` — root cause of C1 |
| W1.2 `_realized_decorated_bex_layer_spacing` in `batch.py` uses shared measurer | ❌ | `batch.py:1908` still does `matmul_vec(pose.rotation_matrix, local_position)[2]` with no `pose.translation` |
| W1.2 `_annotate_ring_embedding` in `ring_forming.py` uses shared measurer | ⚠️ | Already includes `pose.translation` (`ring_forming.py:574`); not refactored to call `measure_layer_z_span` but produces correct results |
| W1.3 `enumerate_candidate_stackings` `monomer_specs` kwarg | ✅ | Keyword-only; both call sites updated: `batch.py:4600` passes `{first.id: first, second.id: second}`; `ring_forming.py:104` passes `{monomer.id: monomer}` |
| W1.3 Silent `0.0` fallback replaced by skip-with-warning | ✅ | `stacking.py:503-550` (`_measure_span`): fresh measurement; falls back to `precursor_coordinates` mode; if unavailable records `stacking_skipped:span-unavailable` |
| W1.4 `c_axis_semantics: "vacuum_slab"` for single-layer exports | ❌ | Not started; `RingFormationConfig.layer_spacing` and `default_ring_layer_spacing` still undocumented; no `c_axis_semantics` key emitted |

---

## W2 — Registry correctness and cell-basis normalisation

| Item | Status | Notes |
|---|---|---|
| W2.1 `classify_2d_cell` in `geometry.py` | ✅ | ±0.5° angle tolerance, 1% relative length tolerance; returns `(kind, setting)` with `"60deg"` / `"120deg"` disambiguation for hexagonal |
| W2.1 `batch.py` `_single_node_cell_kind` migrated | ✅ | `batch.py:5003` now calls `classify_2d_cell`; old cosine-tolerance logic removed |
| W2.1 `embedding.py` `_cell_kind_from_vectors` helper added | ✅ | `embedding.py` has `_cell_kind_from_vectors` calling `classify_2d_cell` |
| W2.1 `single_node_topologies.py` `_metric_family` migrated | ❌ | Still uses its own independent classifier at `single_node_topologies.py:336` |
| W2.1 `indexed_topology_layouts.py` `_metric_family` migrated | ❌ | Still uses its own independent classifier at `indexed_topology_layouts.py:171` |
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
| W4.1 Wider search radius for clash detection | ✅ | `nonbonded_heavy_search_radius: float = 3.5` in `CoarseValidationThresholds`; `_min_nonbonded_heavy_distance_below_cutoff` uses it |
| W4.1 vdW-ratio criterion `d / (r_i + r_j) < 0.75` | ⚠️ | Threshold field `min_nonbonded_heavy_vdw_ratio = 0.75` added; plain-distance backstop `hard_min_nonbonded_heavy_plain_distance = 2.2` added; but the actual per-pair vdW-radius lookup and ratio comparison are not yet wired into the loop — only the wider search radius was applied |
| W4.2 `c2c > layer_z_span` self-check | ✅ | `stacking.py:204-213` warns + records when `center_to_center_distance <= layer_z_span` |
| W4.2 `min_interlayer_contact` computed and flagged | ✅ | `_min_interlayer_contact()` at `stacking.py:552-589`; `stacking_clash` flag on `candidate.flags` when contact < 2.0 Å |
| W4.3 `stacking_skipped:<reason>` flag + stderr warning | ✅ | `stacking.py:128-143`; eligibility failures attach flag to candidate and warn on stderr |
| W4.4 c-axis orthogonality guard | ✅ | `_c_axis_is_orthogonal()` at `stacking.py:660-673`; warns when c is not ⊥ to ab within 2° |

---

## W5 — Ranking integrity

| Item | Status | Notes |
|---|---|---|
| W5.1 `bridge_geometry_residual` doubled on stacking | ✅ | `stacking.py:337-343` doubles `score_metadata["bridge_geometry_residual"]` alongside `total_residual` |
| W5.1 Ranking parity test | ❌ | No test asserting `candidate_ranking_key` is invariant under stacking expansion |
| W5.2 Registry sort order documented | ❌ | `topology_rank`/`stacking_variant_index` separation not done; not documented in `docs/building.md` |

---

## W6 — Latent trap hardening

| Item | Status | Notes |
|---|---|---|
| T1 `ring_center_fractional` w≠0 warning | ✅ | `stacking.py:269-284`; warns + records in `stacking_warnings` |
| T2 Torsion key passthrough documented | ✅ | `stacking.py:302-303`; comment documents passthrough-only, no per-layer suffixing yet |
| T3/T4 Span/axis consistency | ✅ | Covered by W1.2 — `measure_layer_z_span` projects onto `c_hat` not bare z |
| T5 Degenerate c-axis fallback warns | ✅ | `_safe_normalize` at `stacking.py:686-698` |
| T6 Fitted γ = 60.09° cell classified correctly | ✅ | `classify_2d_cell` ±0.5° tolerance handles it |

---

## W7 — Tests, docs, changelog

### Tests

| # | Description | Status |
|---|---|---|
| 1 | TAPB/TFB hcb AA: `c == 2*(clearance+span)`, span > 3 Å, min contact > 2.5 Å | ❌ |
| 2 | Span-never-zero gate: candidate without `embedding.layer_z_span` gets measured span | ❌ |
| 3 | Hexagonal AB shift: vertex→pore-center in both 60° and 120° cells | ❌ |
| 4 | `classify_2d_cell` on fitted γ=60.09° cell returns `("hexagonal", "60deg")` | ❌ |
| 5 | Ranking parity: `candidate_ranking_key` invariant under stacking expansion | ❌ |
| 6 | Validator flags 2.0 Å heavy-heavy nonbonded pair | ❌ |
| 7 | Stacking self-check flags original shipped interpenetrated geometry | ❌ |
| 8 | `ring_center_fractional` w≠0 warns (T1) | ❌ |
| 9 | Metadata + CIF carry derivation block; lammps/graspa round-trips pass | ❌ |
| 10 | `stacking_skipped:<reason>` + stderr warning when eligibility fails | ❌ |
| 11 | `test_ring_forming.py:214` AB assertion made setting-aware | ❌ |
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
| 4. `candidate_ranking_key` parity test green | ❌ |
| 5. CIFs carry `# stacking-geometry:` derivation line; metadata carries full W3.1 schema | ✅ |
| 6. `uv run pytest -q` green; ruff gate; `uv build --wheel` green; all doc/skill/changelog updates landed | ⚠️ — existing tests pass; new integration tests not yet written; docs/changelog not updated |

---

## Highest-priority remaining work

1. **W1.2 `batch.py` span measurers** — `_decorated_bex_layer_spacing` (`batch.py:1869`) and `_realized_decorated_bex_layer_spacing` (`batch.py:1908`) do not add `pose.translation[2]`. This is the stated root cause of C1. Fix: add `+ pose.translation[2]` (or delegate to `measure_layer_z_span`).

2. **W4.1 vdW-ratio loop** — `_min_nonbonded_heavy_distance_below_cutoff` has the wider search radius but does not yet compute `d / (r_vdw_i + r_vdw_j)` per pair. Threshold fields are present; the per-element radius lookup needs wiring in.

3. **W7 integration tests** — tests 1–11 from the plan are all absent. Tests 1, 2, 5, and 7 directly guard the four critical flaws.

4. **W2.1 remaining classifiers** — `single_node_topologies.py:336` and `indexed_topology_layouts.py:171` still have independent `_metric_family` implementations; both should delegate to `classify_2d_cell`.

5. **Docs and changelog** — all seven doc files listed above are unmodified.

---

## Post-audit update — 2026-09-25 (commits cf9776d, 408fbec, 6c966e8)

The builder edge-case audit (`agent-docs/BUILDER_EDGE_CASE_AUDIT.md`) re-bucketed
this plan's remaining items. Current accounting:

**Implemented and kept:** W1.1 docstring; W1.2 shared measurer (helper only —
call sites unmigrated); W1.3 `monomer_specs` + skip-with-warning (caveat: the
`stacking_skipped:span-unavailable` flag promised in the `stacking.py:123`
docstring is never attached — only eligibility skips tag flags); W2.1 canonical
`classify_2d_cell` (two modules unmigrated); W2.2 setting-aware AB shift;
**W2.3 now ✅** (test asserts the 120°-setting `(1/3, 2/3)` shift);
W3.1 canonical metadata schema; W3.2 CIF derivation comment; W4.1 wider search
radius only; W4.2/W4.3/W4.4 warn-only shells; W5.1 residual doubling (kept as a
known symptom patch pending the sum-vs-mean contract, audit A9); W6 T1/T2/T5/T6.
Hardening beyond the plan: per-candidate exception isolation in batch stacking
expansion and narrowed exception swallowing with stderr warnings (408fbec).

**Strategically cancelled** (Tier 3 drops, 6c966e8): `StackingSummary` and its
routing (W3.1's ❌ closed by deletion — the raw dict stays the channel); the
"kept for one release" deprecated aliases (never shipped — removed outright);
`tests/test_stacking.py` (the only ✅ test row above — it covered the deleted
dead API); the `stacking_considered` / `stacking_penalty` stamps (removed with
the legacy scorer).

**Awaiting implementation:** W1.2 root cause (`batch.py:1860-1877` still omits
`pose.translation` — audit A7, highest priority); W1.4 `c_axis_semantics`;
W2.1's two remaining classifiers (audit A8); W4.1 vdW wiring (audit #11 —
threshold fields at `validation.py:37-41` still declared-but-unwired, docstring
still claims the comparison); a real atomistic contact metric to replace the
mathematically dead W4.2 check (audit #12); the proper W5.1 fix + parity test
(audit A9); W5.2; and W7's 11 integration tests plus user-facing docs.
Definition-of-done criteria 1-4 and 6 remain unmet.
