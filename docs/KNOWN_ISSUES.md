# Known Technical Issues

_Audit date: 2026-08-17. Scope: software-engineering only — scientific/chemical
correctness was explicitly out of scope for this audit._

This document tracks known technical issues in the codebase that have **not**
been fixed yet, ordered by severity. Each entry lists the evidence location,
the impact, and a suggested fix direction.

For fixed issues, see the bottom of this file.

## High

### H1. Topology graph matching silently gives up on large unit cells

- **Location**: `src/cofkit/decompose.py:1832` (`_periodic_graphs_match`)
- **Issue**: when the observed graph has `node_count > 10`, the matcher returns
  `False` instead of reporting "unsupported size". Downstream topology
  detection treats this as a definitive mismatch.
- **Impact**: silent false negatives for large unit cells — a structure whose
  topology actually matches is reported as non-matching, with no indication
  that the comparison was never attempted.
- **Suggested fix**: return an explicit "match not attempted (graph too
  large)" outcome (e.g. a tri-state or a metadata flag on
  `TopologyDetectionCandidate`) instead of a bare `False`.

### H2. Monomer build cache is keyed by record id only

- **Location**: `src/cofkit/batch.py:431-460` (`BatchStructureGenerator.build_monomer`)
- **Issue**: `self._monomer_cache` is keyed solely by `record.id`. Two records
  sharing an id but differing in SMILES or motif kind (possible when merging
  libraries from different sources) collide silently.
- **Impact**: the wrong monomer geometry is reused for the colliding record,
  producing incorrect structures without any error.
- **Suggested fix**: include SMILES and motif kind in the cache key, or verify
  consistency on a hit and raise on conflicts.

### H3. Top-N result retention is O(n²·log n)

- **Location**: `src/cofkit/batch.py` (`_record_batch_pair_results`, the
  `top_results` append + full `sort` + truncate per accepted summary)
- **Issue**: every accepted summary triggers a full list sort followed by
  truncation to `retain_top_results`.
- **Impact**: noticeable hotspot in large screening runs
  (`batch-all-binary-bridges` over 10⁴–10⁵ pairs).
- **Suggested fix**: use `heapq.nlargest` semantics or accumulate unsorted and
  sort once at the end of the run.

### H4. Tests exercise the source tree, not the installed package; no test CI

- **Location**: every file under `tests/` starts with
  `sys.path.insert(0, ... / "src")`; `.github/workflows/` contains only
  `release.yml`.
- **Issue**: the test suite never imports the built/installed distribution, so
  packaging-layout regressions (broken default paths, missing data files,
  broken entry points) are invisible to tests. No CI job runs pytest at all.
- **Impact**: packaging bugs ship unnoticed (the `parents[2]` default-path bug
  fixed in this audit is a direct example).
- **Suggested fix**: add a CI workflow that builds the wheel, installs it into
  a clean environment, and runs the suite (plus a CLI smoke test) from a
  directory outside the repo.

## Medium

### M1. `batch.py` has outgrown a single file and contains parallel duplicate implementations

- **Location**: `src/cofkit/batch.py` (~5,200 lines). Examples:
  - `generate_pair_candidate` vs `generate_monomer_pair_candidate`
    (~75 nearly identical lines each);
  - `_build_expanded_single_node_node_linker_outcome` vs
    `_build_indexed_node_linker_outcome`, and the corresponding `_embed_*`
    pair (~140 lines each, differing mostly in node-site metadata);
  - the outcome/embed builder family repeats the same
    slot/instance/event/net-plan assembly pattern per builder.
- **Impact**: parallel copies have already started drifting apart (see L6);
  every behavior change must be mirrored across copies, which is a regression
  incubator.
- **Suggested fix**: extract the shared outcome/embedding assembly into one
  parameterized implementation; split topology builders into a dedicated
  module package.

### M2. `cif_export_index` metadata is meaningless in parallel batch mode

- **Location**: `src/cofkit/batch.py` (`_collect_parallel_pair_results` runs
  tasks with the default `cif_export_start_index=0`; only the serial branch
  threads the running `cifs_written` counter through)
- **Impact**: `manifest.jsonl` records repeat export indices in parallel runs;
  downstream consumers cannot rely on the field.
- **Suggested fix**: either assign indices after collection (single-threaded,
  deterministic) or document the field as serial-only.

### M3. Broad `except Exception` blocks convert programming errors into normal results

- **Location**: 19 sites, notably `src/cofkit/decompose.py:474,604,956,1421,1760`,
  `src/cofkit/decompose_events.py:243,1376`, `src/cofkit/batch.py:452,1154,1240`.
- **Issue**: the result-object pattern is intentional, but `TypeError` /
  `KeyError` / `AttributeError` from genuine bugs are also swallowed into
  `"status": "error"` results, usually without logging.
- **Impact**: real defects surface as routine per-item failures and are hard
  to diagnose from aggregate output.
- **Suggested fix**: catch expected exception types explicitly; in the
  catch-all path, record/log the full traceback (e.g. `logging.exception`)
  before converting to an error result.

### M4. Generated simulation scripts embed unescaped file names

- **Location**: `src/cofkit/lammps.py:3725,3737,3866,3873,3901,3912`
  (`read_data {data_file}`, `dump ... {dump_file}`, `write_dump ...`);
  gRASPA `simulation.input` rendering embeds `framework_name` similarly.
- **Issue**: base names are interpolated into LAMMPS/RASPA scripts without
  quoting or validation. A CIF file name containing spaces or shell/script
  metacharacters silently misparses.
- **Impact**: confusing failures (or wrong atom counts) for perfectly legal
  input file names.
- **Suggested fix**: validate/sanitize file names up front (whitelist
  `[A-Za-z0-9_.-]`) and fail fast with a clear message, or quote per the
  target tool's syntax.

### M5. Dependency declaration and import guards disagree; facade eagerly imports heavy stacks

- **Location**: `pyproject.toml` lists `gemmi`, `rdkit`, `openbabel`,
  `pandas`, `pymatgen` as hard dependencies, while `src/cofkit/lammps.py:34-60`
  (and others) guard the same imports as optional; `src/cofkit/__init__.py`
  re-exports nearly every submodule, so any `import cofkit.xxx` pulls the full
  dependency chain.
- **Impact**: environments cannot actually run a minimal install despite the
  guards, and lightweight use cases pay the import cost of pymatgen/pandas.
- **Suggested fix**: either move heavy deps to extras and make the facade
  lazy (`__getattr__`-based re-exports), or drop the guards and treat the
  dependencies as required everywhere.

## Low

### L1. Planar motif rotation averages wrapped angles arithmetically

- **Location**: `src/cofkit/batch.py` (`_rotation_for_planar_motifs`, the
  `sum(diffs) / len(diffs)` mean of `_wrap_angle` outputs)
- **Issue**: arithmetic mean of circular quantities is wrong near ±π
  (e.g. +179° and −179° average to 0° instead of 180°).
- **Impact**: occasionally suboptimal monomer orientations.
- **Suggested fix**: use vector averaging (sum of unit phasors) for circular
  means.

### L2. Unbounded in-memory caches

- **Location**: `src/cofkit/batch.py:247-249` (`_spatial_rotation_cache`,
  `_topology_id_cache`)
- **Impact**: slow memory growth in long-lived processes screening large
  libraries.
- **Suggested fix**: bound with `functools.lru_cache` or an explicit size cap.

### L3. SMARTS patterns recompiled per call

- **Location**: `src/cofkit/chem/rdkit.py:117` (`Chem.MolFromSmarts` inside
  `_detect_motifs`)
- **Impact**: repeated pattern compilation for every monomer; measurable in
  large library scans.
- **Suggested fix**: cache compiled patterns per motif kind.

### L4. `_json_safe` does not handle sets or dataclasses

- **Location**: `src/cofkit/batch.py` (`_json_safe`)
- **Impact**: such values degrade to `str(value)` in JSON output.
- **Suggested fix**: add `set`/`frozenset` (sorted list) and
  `dataclasses.asdict` handling.

### L5. Tooling gaps

- No ruff/mypy configuration; scattered `# type: ignore` comments
  (e.g. `src/cofkit/batch.py`); `dev` extra contains only `pytest`; no
  pre-commit hooks; no lint/typecheck CI.
- **Suggested fix**: add ruff (+ mypy in warn mode) config and a minimal
  pre-commit/CI check.

### L6. `net_plan` metadata hardcodes `topology_family: "single-node-2d"` for 3D builds

- **Location**: `src/cofkit/batch.py:2200,2496` (expanded node-linker and
  node-node outcome builders) hardcode the string, while the corresponding
  embed functions (`batch.py:2432,2647`) correctly compute
  `f"single-node-{dimensionality.lower()}"`. (`batch.py:1911,2119` are in the
  2D-only legacy builder and are correct.)
- **Impact**: mislabeled metadata for 3D single-node builds; the two parallel
  implementations already disagree — direct evidence for M1.
- **Suggested fix**: compute from `topology.dimensionality` like the embed
  path does.

### L7. Repeated full library reloads during template discovery

- **Location**: `src/cofkit/batch.py:404-426`
  (`available_binary_bridge_template_ids` reloads the whole library tree per
  candidate template)
- **Impact**: redundant IO and monomer re-inference during
  `batch-all-binary-bridges` discovery.
- **Suggested fix**: load once and reuse across templates.

### L8. CLI error surfacing is inconsistent

- Deep failures propagate as raw tracebacks in some commands while others
  use friendly `SystemExit` messages (e.g. the input-dir checks added in
  `cli_build.py`).
- **Suggested fix**: adopt a top-level handler in `cli.main` that converts
  expected `ValueError`/`FileNotFoundError` into concise stderr messages with
  non-zero exit codes, keeping tracebacks behind a `--debug` flag.

## Fixed in the 2026-08-17 audit round (critical severity)

Kept for reference; these are resolved and covered by regression tests:

1. **`engine.py` unbound `connectivity` NameError** — 3D node-node
   equal-connectivity projects without an explicit topology crashed through
   the documented `COFEngine.run()` API. Fixed by binding
   `connectivity = connectivities[0]` in the node-node branch
   (`src/cofkit/engine.py`). Regression test:
   `tests/test_core.py::EngineTests::test_engine_supports_3d_node_node_project_without_explicit_topology`.
2. **Unguarded `shutil.rmtree(output_dir)`** — `cofkit build default-library`
   deleted any user-specified output directory without confirmation. Fixed
   with a marker-file guard plus an explicit `--force` flag
   (`src/cofkit/cli_build.py`, `_prepare_default_library_output_dir`).
   Regression tests: `tests/test_cli.py::DefaultLibraryOutputGuardTests`.
3. **`parents[2]` default paths broken when installed** — 8 CLI argument
   defaults resolved into the Python installation directory for packaged
   installs. Fixed via `_default_repo_path()` (repo-root in a source
   checkout, cwd-relative when installed) plus input-dir existence checks.
   Regression tests: `tests/test_cli.py::DefaultRepoPathTests`.
4. **Process-pool fallback flaws in batch generation** — `BrokenProcessPool`
   (a `RuntimeError`) was not caught, so a dead worker crashed the whole run;
   and a mid-iteration pool failure caused already-recorded results to be
   recorded twice by the thread-pool fallback. Fixed by catching
   `BrokenProcessPool` and buffering results before recording
   (`src/cofkit/batch.py`, `_collect_parallel_pair_results`). Regression
   tests: `tests/test_batch.py::BatchProcessPoolFallbackTests`.

## Resolved separately (2026-08-17): legacy candidate score disabled by default

The `CandidateScorer` total score (`src/cofkit/scoring.py`) was a heuristic
whose magnitude was dominated by the number of reaction events in the unit
cell — measured correlation between score and event count was **1.0** on a
TAPB+TFB all-topology run, so the ranking systematically preferred larger
unit cells over better bridge geometry (the geometrically near-perfect `hcb`
candidate ranked last). The score is now **disabled by default**:

- `Candidate.score` / `BatchPairSummary.score` are `None` and outputs no
  longer contain `score_breakdown`; candidate and summary metadata carry
  `scoring_mode: "residual"`.
- All ranking/selection (best-candidate selection, topology ordering,
  top-results retention, ensemble ordering) uses the mean per-bridge-event
  geometry residual (lower is better), via
  `cofkit.model.candidate_ranking_key` / `order_candidates`.
- Opt-in restore: `BatchGenerationConfig(enable_legacy_scoring=True)` /
  `COFEngineConfig(enable_legacy_scoring=True)` /
  `RingFormationConfig(enable_legacy_scoring=True)`, or `--legacy-scoring`
  on the build CLIs, which re-attaches the score and restores score-based
  ranking exactly as before.
- Not affected: the coarse validator and the continuous optimizer never used
  the total score — they consume the per-event bridge residuals, which are
  still always computed in `score_metadata`.

Regression tests: `tests/test_batch.py::ScoringModeTests`,
`tests/test_core.py::EngineTests::test_imine_project_legacy_scoring_restores_total_score`.
