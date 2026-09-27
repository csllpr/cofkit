# Magic-Number & Single-Owner Fix Plan

Date: 2026-09-27. Companion to `BUILDER_EDGE_CASE_AUDIT.md`; supersedes nothing
there — this plan covers the *wiring* disease, not individual algorithmic
edge cases.

## Diagnosis

A two-agent sweep of `src/cofkit` (2026-09-27) found **69 distinct sites
(~110+ literals)** that are both non-intuitive and non-referenced, plus **9
groups** whose comments openly admit uncalibrated-heuristic status.

The count itself is not the real problem. The boronate-ester `range(721)` /
`720.0` scan surviving the #15 rework (imine/azine scan deleted, its twin
left behind at `reaction_realization.py:1295`) demonstrates the structural
defect: **the same quantity exists in several places, owned by no one**, so
fixes land on one copy and silently miss the others. Concrete instances
found today:

- c-axis vacuum slab `8.0` — five literal owners (`embedding.py:41`,
  `model.py:127`, `engine.py:33`, `ring_forming.py:57`, `cli_build.py:509`),
  no shared constant. The W1.4 "unification" unified the *value*, not the
  *owner*.
- `cli_analyze.py:47–58` retypes all nine `CoarseValidationThresholds`
  values as argparse defaults — guaranteed to drift from
  `validation.py:33–41`.
- default bridge target distance: `1.5` (`reactions.py:121,192`) vs `1.4`
  (`embedding.py:43`) — same quantity, two values, chosen by entry path.
- covalent-radius fallback: `0.7` (`chem/detector.py:66`) vs `0.75`
  (`decompose.py:3524`), with two incompatible bond-tolerance models
  (multiplicative ×1.3 vs additive +0.30/+0.45 + B–O special case).
- RASPA/EQeq/cutoff defaults retyped across `graspa.py`, `hybrid_mdmc.py`,
  and ~20 argparse sites in `cli_calculate.py` (12.8 Å ×16, EQeq params ×5,
  seed 246813 ×4, Ewald 1e-6 ×11).
- ring bond lengths 1.38/1.35 in three modules; angular residual weight
  30.0 in two; `_FRACTIONAL_WRAP_TOLERANCE` defined twice.

## Governing principle (proposed convention)

> Every tunable numeric has **exactly one owner**: a named constant in the
> module that semantically owns the quantity, or a field on the owning
> config dataclass. Every owner carries a provenance comment of one of
> three kinds: **derived** (formula shown), **cited** (literature/source),
> or **heuristic — pending calibration** (explicitly labeled). CLI argparse
> defaults must reference the owner, never retype a literal. The same
> quantity must never exist with two values.

## Workstreams

### W1 — Single-owner wiring fixes (mechanical, behavior-preserving)

Each item: introduce/choose one owner, rewire all consumers, delete
duplicate literals. No value changes except where flagged for decision.

| # | Item | Sites | Notes |
|---|------|-------|-------|
| W1.1 | Cell-default constants (`DEFAULT_MONOLAYER_C_ANGSTROM = 8.0`, `DEFAULT_LATERAL_SPAN_ANGSTROM = 30.0`) | embedding.py, model.py, engine.py, ring_forming.py, cli_build.py | finishes W1.4 of the stacking plan properly |
| W1.2 | `cli_analyze.py` argparse defaults ← `CoarseValidationThresholds` fields | cli_analyze.py:47–58 | removes 9 retyped literals |
| W1.3 | `cli_calculate.py` argparse defaults ← `graspa.py`/`lammps.py` config dataclass defaults | ~20 sites (12.8 cutoff, EQeq, seed, Ewald, temperatures) | biggest literal-count reduction in the repo |
| W1.4 | One covalent-radius table + one fallback | chem/detector.py, decompose.py:3513–3524 | **value conflict 0.7 vs 0.75 — owner decision (W3)**; also reconcile or document the two tolerance models |
| W1.5 | Bridge target distances single home | reactions.py profiles vs embedding.py:43 | **value conflict 1.4 vs 1.5 — owner decision (W3)** |
| W1.6 | Ring bond lengths (1.38/1.35) + angular weight (30.0) single home | reactions.py, ring_geometry.py, reaction_realization.py, stacking.py:875 | |
| W1.7 | `_FRACTIONAL_WRAP_TOLERANCE` shared | topology_analysis.py:26, topology_symmetry.py:7 | trivial |

### W2 — Finish the #15 rework: kill the boronate 721-scan

`reaction_realization.py:1295–1296` still runs a 721-step angular scan with
unexplained constants (0.01 O–B–O weight at :1321) for boronate-ester —
the exact pattern the imine/azine rework removed.

- **Route A (preferred, mirrors #15):** extend the `BridgeGeometryPriors`
  closed-form constructor to boronate-ester (B–O 1.38 already at
  `chem/linkage.py:28`; O–B–O prior 120° sp2). Delete the scan. Run the
  same clash-regression screen used for the imine rework (baseline vs
  reworked, `out/_rework_baseline` pattern).
- **Route B:** if investigation shows the scan is load-bearing for
  boronate (e.g. genuinely two-degree-of-freedom fit), extract its
  constants as named, derived values and document why this template
  differs.
- First step is a short investigation to pick the route.

### W3 — Value-conflict decisions (owner input required)

| # | Conflict | Options |
|---|----------|---------|
| W3.1 | covalent-radius fallback 0.7 vs 0.75 | pick one; 0.7 matches the Cordero-table minimum region, 0.75 is decompose's historical value |
| W3.2 | default bridge distance 1.4 vs 1.5 | 1.5 is the reactions.py fallback; 1.4 is embedding's. Realistically templates override both, so this is about which fallback is defensible |
| W3.3 | default temperature 300 vs 298 K | likely intentional per-command (Widom vs adsorption); if so, document, don't unify |
| W3.4 | soft_relax clash floors vs vdw.py | route through vdw.py constants, or document as an intentionally independent coarse model |

### W4 — Provenance documentation pass (admitted heuristics)

For each cluster: name the constants, add derived/cited/heuristic
provenance comments. No behavior change.

- **W4.1** `reaction_realization.py` fitting machinery (~20 sites):
  benzothiazole block (:778–861), azine fit weights (:1770–1843),
  H-reposition (:3628–3905). Benzothiazole is internal-only per
  AGENTS.md — label `un calibrated, internal prototype` and freeze.
- **W4.2** `decompose*.py` scoring weights and bond-order windows
  (~10 sites).
- **W4.3** `validation.py` nine thresholds + `engine.py`/`lammps.py`/
  `soft_relax.py`/`graspa.py` config defaults (~12 sites).
- **W4.4 (optional, needs decision):** calibrate the two flagged
  carry-overs against the CoRE-COF subset in `files_for_reference/`:
  azine N–N 1.408 and stacking clearances 3.4/3.5/3.6. This changes
  values, not just docs — treat as its own decision.

### W5 — Guardrails (prevent regression)

| # | Item |
|---|------|
| W5.1 | Extend the grep-gate pattern from `tests/test_linkage_geometry.py`: forbidden literals outside owner modules (c-axis `8.0` in cell defaults, deleted retraction constants, `range(721)`) |
| W5.2 | Single-owner assertion test: CLI argparse defaults equal the owning dataclass fields (fails on retyped literals) |
| W5.3 | AGENTS.md: add the governing-principle convention above to "Constraints and conventions" |
| W5.4 | `AGENT_CODEBASE_MAP.md`: list the constant-owner modules so future agents land new tunables in the right place |

## Proposed tiering for fix/drop decisions

- **Tier 1 (do now, mechanical, zero behavior change):** W1.1–W1.3,
  W1.6, W1.7, W5.3, W5.4.
- **Tier 2 (finish the rework):** W2 — investigate, then Route A or B.
- **Tier 3 (owner decisions, small):** W3.1–W3.4, then W1.4/W1.5
  implementations once decided.
- **Tier 4 (documentation pass, medium effort, zero behavior change):**
  W4.1–W4.3. W4.4 is a separate calibration decision.
- **Guardrail tests W5.1/W5.2** land with the Tier-1 wiring they protect.

## Verification plan

- Every workstream: `uv run pytest -q` (baseline 603 passed / 4 skipped)
  + `uvx ruff==0.12.12 check src --select E9,F63,F7,F82`.
- W1: assert no behavior change — exported CIFs from a fixture build are
  byte-identical before/after (except W2, which has its own clash screen).
- W2: baseline-vs-reworked clash regression screen, mirroring commit
  `7f893c3`'s evidence in `out/_rework_baseline/`.
- CHANGELOG.md entry per tier; this document carries resolution notes per
  row as work lands.

## Resolution log

**2026-09-27 — Tier 1 landed (commit `d5391d3`).** W1.1: `cofkit.constants`
created (`DEFAULT_MONOLAYER_C_ANGSTROM`, `DEFAULT_LATERAL_SPAN_ANGSTROM`),
all five 8.0 owners and both 30.0 owners rewired. W1.2: all 12
classify-output argparse defaults read from `CoarseValidationThresholds`.
W1.3: ~118 argparse sites across the five `calculate` subcommands rewired
to the owning dataclasses; graspa-internal duplicates hoisted to
`DEFAULT_CUTOFF_ANGSTROM` / `DEFAULT_OVERLAP_CRITERIA` /
`DEFAULT_EWALD_PRECISION` / `_KCAL_PER_MOL_TO_KELVIN`. W1.6: ring bond
lengths owned by `ring_geometry` (`BOROXINE_BO_BOND_LENGTH`,
`TRIAZINE_CN_BOND_LENGTH`), angular weight owned by
`geometry.ANGULAR_RESIDUAL_DOWN_WEIGHT`. W1.7: wrap tolerance single
definition in `topology_symmetry`. W5.3/W5.4: AGENTS.md convention +
codebase-map owner list. W5.2: `tests/test_cli_defaults_wiring.py`.
W1.4/W1.5 remain blocked on W3.1/W3.2 owner decisions. Left in place
deliberately: `--timeout-seconds` 300.0 literals (owned by workflow
function signatures, not config dataclasses — candidate for a later
micro-consolidation), `HybridMdMcSettings` cutoffs (per-command owner,
same pattern as the 300/298 K split), `guest_restart.py:688-689` retyped
fallbacks (out of scope, minor).

**2026-09-27 — Tier 2 landed (Route A).** Investigation verdict: the
boronate fit is genuinely 1-DOF (B rotation α; oxygens exact via
circle-circle intersection) — the scan was not load-bearing. Rework
mirrors the imine/azine constructor: bracketed golden-section on the pure
angle residual (analytic feasible intervals + branch-kink splitting;
brackets proved necessary against a multi-well residual), exact closure at
the priors when feasible, honest best-effort diagnostics otherwise. New
constants: `linkage_geometry.BORONATE_OBO_EXACT_TOLERANCE_DEG` (derived),
`reaction_realization.BORONATE_GOLDEN_SECTION_{BRACKETS_PER_PIECE,ITERATIONS}`
(heuristic-labeled). Placement distance now derived (`t*cos(θ/2)` =
0.8017). Corrections to this plan's Route-A text recorded: boronate's
realized B–O target is 1.44 (1.38 is boroxine's); `reactions.py:434`'s
1.36 is keto-enamine's; `chem/linkage.py:28` is a dead 1.38 duplicate
(follow-up cleanup candidate). Regression screen
(`out/_rework_boronate/EVIDENCE.md`): ∠OBO 113.8→112.4° exact at all six
BDBA×HHTP ring sites, no new clashes or validation regressions, cells
+0.012%, LAMMPS fixed-cell converges ~89 kcal/mol closer to the minimum.
W5.1 grep gate extended to the retired scan constants.

