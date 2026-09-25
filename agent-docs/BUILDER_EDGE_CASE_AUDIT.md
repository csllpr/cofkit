# Builder edge-case audit — fix-or-drop decision document

Audited 2026-09-25. Method: four parallel read-only scouts over the builder
subsystems (reaction realization / chemistry, stacking+embedding+geometry,
build workflows+planner+engine, validation+soft-relax+scoring+optimizer),
cross-verified against `STACKING_FIX_PLAN.md` / `STACKING_FIX_REVIEW.md` /
`STACKING_IMPL_STATUS.md`, with live experiments on real `out/` CIFs for the
validator findings.

**Purpose:** the builder accumulates edge-case patches that do not fix, or only
partially fix, the underlying problems. This document lists each cluster,
scores confidence in a genuine fix, and marks fix vs drop so each can be
actioned deliberately. Tier 1 (fix, high confidence) and Tier 3 (drop) have
action lists at the bottom.

## Meta-patterns (why partial fixes accumulate)

1. A canonical helper gets written and the migration stalls halfway —
   `geometry.py:classify_2d_cell` and `geometry.py:measure_layer_z_span` each
   still have 3–4 divergent local copies in live use.
2. Safety rails are declared but inert: threshold fields nothing reads (W4.1),
   a self-check that mathematically cannot fire (W4.2), validator checks that
   parse labels wrong and silently never run.
3. Warn-then-ship-anyway (tilted-c cells) or swallow-without-warning, violating
   the repo's own "warn, don't silently degrade" principle.
4. Tests cover dead code while live machinery is untested; the working tree
   carried a red test (`tests/test_ring_forming.py:212`).
5. Same-day fix-commit pairs (`5a11e02`+`44f72fc`, `7e26fba`+`17f9d8f`,
   `4fc403e`+`2e24179`) indicate symptom whack-a-mole.

---

## Tier 1 — Fix, high confidence (root cause known, bounded effort)

### 1. Layer z-span computed four divergent ways (~90%)
`batch.py:1869` and `batch.py:1908` omit `pose.translation` — documented root
cause (C1) of shipped interpenetrated stacks, still unfixed while downstream
patches (soft-relax, W4.2, W5.1) accumulate on top. `ring_forming.py:574`
uses bare `[2]`; `geometry.py:158 measure_layer_z_span` is canonical.
Metadata stamps `layer_z_span_axis: "c_hat"` even when the number came from a
bare-z fallback (false provenance). Fix: delegate all sites to
`measure_layer_z_span`; per the review, measure along the layer normal
`a×b`, not c_hat.

### 2. Cell classification split-brain (~85%)
4–5 classifiers with different tolerances and contradictory hexagonal
conventions: topology defs are γ=120°, embedding builds γ=60°, and
`embedding.py:561` only detects the 60° form (the T6 bug).
`geometry.py:246 classify_2d_cell` is canonical; local copies at
`single_node_topologies.py:336`, `indexed_topology_layouts.py:171`,
`embedding.py:363,561`. The W2.1 "done" in embedding is hollow:
`_cell_kind_from_vectors` (`embedding.py:368`) has zero call sites.
Fix: delete local classifiers, route through `classify_2d_cell`, pick one
canonical cell setting per net at build time.

### 3. Validator `_instance_id` parses labels wrong (~95%)
`validation.py:499-501` uses `partition("_")`; instance ids contain
underscores, so on a real `bex_d2h` CIF all 183 bonds collapse into one
pseudo-instance. Consequences: the disconnected-instance-graph check can
never fire, and the boronate B–O bond-window check (advertised in
`docs/CURRENT_SCOPE.md:43`) never runs on real output.
`decompose.py:3563-3566` parses the same labels correctly with
`rsplit("_", 1)`. Empirically demonstrated.

### 4. Silent fit failures with metadata claiming success (~95%)
Imine chain-closure fit (`reaction_realization.py:2220-2274`) returns `None`
on degenerate inputs → bridge exported collinear, but notes unconditionally
claim the fit was applied. Azine fit silently skipped for ≠2-event
hydrazines (`:1648-1650`), for ≠2 reactive atoms (`:1732`), and per-endpoint
failures (`:1840-1846`). The boronate path already does conditional notes
correctly (`:2721-2728`).

### 5. Batch stacking expansion lacks per-candidate exception isolation (~95%)
`batch.py:4597-4605` expands all candidates × registries in one tuple
comprehension; a single exception discards every sibling variant of the
pair. Violates the repo's batch principle; flagged by the stacking review.

### 6. Cross-boundary bonded pairs bypass the validator's exclusion (~90%)
`validation.py:455-462`: the `image_index == 0` gate misses bonds crossing a
cell boundary → false `heavy_atom_clash` positives and polluted
minimum-distance metric (12/400 scanned real CIFs show bond-range
"minima"). Same gap copy-pasted at `soft_relax.py:409-412`.
Fix: exclude on unordered label pairs (labels are unique in P1).

### 7. W5.1 residual doubling compensates a symptom (~85%)
`stacking.py:337-343,651-656` doubles residuals so stacked candidates match
unstacked ranking — scoring and stacking disagree on sum-vs-mean
normalization. Wrapped in `try/except: pass`; `stacking_considered: False`
written onto stacked candidates (`stacking.py:348`). Fix: define the
sum-vs-mean contract in `scoring.py`, recompute aggregate in one place.

### 8. Duplicate primitives with divergent semantics (~90%)
`_safe_normalize` in six modules (`scoring.py:226`, `optimizer.py`,
`batch.py:5327`, `embedding.py:652`, `indexed_topology_layouts.py:276`,
`stacking.py:686`) with two different fabricated fallback directions;
angle-from-three-points ×4 with different degenerate returns;
`_orthogonal_component` twice in one class
(`reaction_realization.py:1568,3852`). Mechanical consolidation into
`geometry.py`.

### 9. Validation ordering hides clashes (~90%)
`validation.py:194-204` with `skip_cif_checks_when_metadata_invalid=True`:
when metadata checks fail, the whole CIF-level block (including the clash
scan) is skipped and the record is demoted to `needs_optimization` —
interpenetrated structures are routed to "repairable" by construction.
Fix: run clash/degeneracy checks unconditionally.

### 10. Swallowed-exception cluster (~95%)
Bare/narrow-less `except` with silent downgrade: `stacking.py:434-435`
(cell classification), `stacking.py:535-541` (span realization),
`ring_forming.py:556-565` (embedding annotation),
`validation.py:491-497` (topology dimensionality → None → wrong
degeneracy check). Plus `validation.py:523-531` `_distance_residual`
returns 0.0 for missing data (missing masquerades as perfect).

---

## Tier 2 — Fix, medium confidence (effort or calibration risk real)

| # | Issue | Confidence |
|---|-------|-----------|
| 11 | W4.1 vdW-ratio clash criterion declared but never wired (`validation.py:34-40` vs `:431-465`); naive wiring false-positives 1–3 angle neighbors → needs periodic-image-aware 1–2/1–3 exclusions first | ~75% |
| 12 | W4.2 self-check mathematically dead: `_min_interlayer_contact` (`stacking.py:552-587`) measures monomer-center translations (≥3.4 Å by construction) yet is exported into CIFs as an atomistic contact; four inconsistent clash/contact definitions across validation/stacking/soft-relax | ~80% |
| 13 | Tilted-c cells: warn and ship wrong geometry anyway (`stacking.py:660-673`); derivation false exactly when the warning fires; 2° guard leaves a silent window | ~75% (lands with #1) |
| 14 | Keto-enamine tautomerization half-realized: ring never re-bond-ordered → valence-5 aromatic carbon in exported CIFs (`reaction_realization.py:2558-2621`); vinylene same shape (`:2804`); inverse algorithm exists in `decompose_bond_orders.py` | ~60% (~95% for the honest-degradation variant: stop writing `_ccdc_geom_bond_type`) |
| 15 | Two independent imine-collinearity compensators (`linkage_geometry.py:7-8` retraction 0.11/0.08 + realization 721-step fit toward uncited 127.2°/127.8° targets); retraction silently absent for mixed-linkage builds; constants don't match literature ~116–121° | ~55% |
| 16 | Soft-relax is a symptom patch over #1 (`soft_relax.py` whole module): cannot fix fixed-cell problems, ignores `converged=False`, "nearest image" is "first image", second independent clash threshold | ~70% as labeled stopgap; drop long-term once #1 lands |
| 17 | Optimizer proposals fight the objective: rotation moves the attachment endpoints it aligns (`optimizer.py:264-322`); sum-vs-mean inconsistency forced the W5.1 patch | ~65% |
| 18 | Conformer quality degradation ladder invisible downstream (`chem/rdkit.py:336-413,522-642`): unminimized monomer indistinguishable from MMFF-optimized | ~75% |
| 19 | Topology-inference heuristics: finite-patch bipartite coloring (`single_node_topologies.py:620-674`), handcrafted rotate/mirror-trigonal (`:288-304`), two-space-group ceiling with silent omission (`:482-493,:98-101`) | ~60% |

---

## Tier 3 — Drop candidates (removal is the fix)

| # | Dead weight | Confidence removal is safe |
|---|-------------|---------------------------|
| 20 | Benzothiazole dead code: ~390 lines never called (`reaction_realization.py:884-1272`); `local_relaxation_applied` hard-coded False; internal-only per AGENTS.md | ~100% (dead functions; the live grid-search path is a separate judgment call) |
| 21 | Dead public stacking API `resolve_stacking_pattern`/`compute_interlayer_offset` (`stacking.py:701-777`): zero call sites, contradicts live Bernal-shift implementation, promises "ABC" the enumerator can't build — and it's the only stacking code with tests | ~95% |
| 22 | Legacy scorer path (`scoring.py:43-98`): deprecated but still invoked unconditionally for metadata; placeholder fields that can never be nonzero | ~90% (keep `bridge_geometry_report`) |
| 23 | Dead planner ring branch (`planner.py:56-80`, divergent default policy vs `engine.py:158-165`); vestigial `stacking_mode` (`engine.py:78-79,323-324`); unrouted `StackingSummary` (`batch_models.py:106-151`); dead `_candidate_layer_z_span` (`stacking.py:478-488`); unreachable negative-span clamp (`stacking.py:193-200`); `geometry.py:194` self-import | ~95% |
| 24 | Deprecated aliases "kept for one release" with calendar versioning and no removal gate | Decision, not effort |

---

## Action list — Tier 1 (fix now)

- [x] A0. Fix red test `tests/test_ring_forming.py:212` — make AB-shift
      assertion setting-aware (W2.3).
- [x] A1. (#3) `validation.py:_instance_id` — parse labels with
      `rsplit("_", 1)`; regression test on instance ids containing
      underscores (assert non-zero inter-instance edges and realized
      bridge-bond count).
- [x] A2. (#4) Imine/azine fits — conditional notes + per-instance skip
      diagnostics in metadata instead of unconditional success claims.
- [x] A3. (#5) `batch.py` stacking expansion — per-candidate try/except
      recording `f"{type(exc).__name__}: {exc}"` per the batch convention.
- [x] A4. (#6) Cross-boundary bonded-pair exclusion in
      `validation.py` and `soft_relax.py` — exclude on unordered label
      pairs; regression test with a boundary-crossing bond.
- [x] A5. (#9) Validation ordering — run clash/degeneracy checks
      unconditionally even when metadata checks fail.
- [x] A6. (#10) Narrow swallowed exceptions; emit standard
      `warning:` stderr lines; fix `_distance_residual` missing-data 0.0.
- [ ] A7. (#1) Span pipeline — route all four z-span sites through
      `measure_layer_z_span`; fix `batch.py:1869,1908` translation
      omission; layer-normal axis; honest `layer_z_span_axis` provenance.
- [ ] A8. (#2) Cell classifiers — migrate all copies to
      `classify_2d_cell`; decide canonical cell setting per net; wire or
      delete `embedding.py:_cell_kind_from_vectors`.
- [ ] A9. (#7) Define sum-vs-mean residual contract in `scoring.py`;
      remove doubling compensation.
- [ ] A10. (#8) Consolidate `_safe_normalize`/angle/distance primitives
      into `geometry.py`.

## Action list — Tier 3 (drop)

- [x] D1. (#21) Delete `resolve_stacking_pattern`/`compute_interlayer_offset`
      and their tests (or replace tests with live-machinery tests).
      Done 2026-09-25: both functions, `__all__` entries, and
      `tests/test_stacking.py` (which covered only this dead API) deleted.
- [x] D2. (#23) Delete dead planner ring branch, vestigial `stacking_mode`,
      dead `_candidate_layer_z_span`, unreachable span clamp,
      `geometry.py:194` self-import; route or delete `StackingSummary`.
      Done 2026-09-25: all deleted (`StackingSummary` deleted, not routed).
      The `stacking_disabled` candidate flag is still emitted unconditionally
      because `stacking.py` consumes it when expanding stackings; constant
      `"stacking_mode"` metadata stamps in `ring_forming.py` / `stacking.py`
      kept as informational. No CLI flag ever exposed `stacking_mode`.
- [x] D3. (#22) Stop invoking legacy `CandidateScorer.score` for metadata;
      drop `total`/breakdown fields.
      Done 2026-09-25: `score()`/`ScoreResult`/magic weights deleted; live
      metadata now comes from `CandidateScorer.scoring_metadata`;
      `enable_legacy_scoring` plumbing removed from all configs. The public
      `--legacy-scoring` CLI flags were kept as accepted deprecated no-ops
      with a stderr warning (repo deprecation convention).
- [x] D4. (#20) Delete uncalled benzothiazole fit/relax functions and the
      hard-coded-False metric.
      Done 2026-09-25: both functions plus the now-orphaned
      `_project_onto_plane` helper and the `local_relaxation_applied`
      metric/parameter removed.
- [x] D5. (#24) Decide removal release for deprecated stacking aliases.
      Done 2026-09-25: git archaeology showed the alias keys were introduced
      in the unreleased checkpoint `cf9776d` (no tag contains it; the
      `2026.9.12` release predates it), so the aliases never shipped and were
      deleted outright; consumers switched to the canonical keys. Note:
      `LayerRegistry.interlayer_distance` is the dataclass's actual field
      name (not an alias property), so it was kept.

## Cross-cutting requirement

Whichever items are fixed: write guard tests first. The current suite passes
over every failure in this document — that is how the patches accumulated.
The stacking plan's 11 unwritten tests (W7 in `STACKING_IMPL_STATUS.md`) are
the concrete checklist for the stacking cluster.
