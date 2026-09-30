# Agent Codebase Map

This document is for later coding/agent sessions that need to move quickly in `cofkit` without rebuilding the mental model from scratch.

## First orientation

Read these files first, in this order:

1. [README.md](../README.md)
2. [docs/README.md](../docs/README.md)
3. [docs/CURRENT_SCOPE.md](../docs/CURRENT_SCOPE.md)
4. [ARCHITECTURE.md](../ARCHITECTURE.md)
5. [CHANGELOG.md](../CHANGELOG.md)

Then go straight to the module that matches the task.

## Main user entry points

- [src/cofkit/cli.py](../src/cofkit/cli.py)
  - Installed `cofkit` CLI.
  - Routes the top-level `build`, `analyze`, `calculate`, and `validate` namespaces.
- [src/cofkit/cli_build.py](../src/cofkit/cli_build.py)
  - Owns build-facing commands such as `cofkit build single-pair` and the batch workflows.
- [src/cofkit/cli_analyze.py](../src/cofkit/cli_analyze.py)
  - Owns analysis-facing commands such as `cofkit analyze classify-output`, `cofkit analyze decompose`, and `cofkit analyze zeopp`.
- [src/cofkit/cli_calculate.py](../src/cofkit/cli_calculate.py)
  - Owns external calculation commands such as `cofkit calculate lammps-optimize`, `graspa-widom`, `graspa-isotherm`, `graspa-mixture`, and `hybrid-mdmc`.
- [src/cofkit/cli_validate.py](../src/cofkit/cli_validate.py)
  - Owns validation commands such as `cofkit validate simple` and `cofkit validate optimize`.
- [src/cofkit/engine.py](../src/cofkit/engine.py)
  - Direct project-style API via `COFEngine`.
- [src/cofkit/batch.py](../src/cofkit/batch.py)
  - Practical execution layer for topology-guided single-pair and batch generation.
  - Failure-isolation seams: per-candidate realization/CIF-export boundary
    `_export_candidate_cif_guarded` (status `export-failed`, staged artifacts
    removed), last-resort per-pair boundary `_run_batch_pair_task` /
    `_failed_pair_task_result` (status `pair-task-failed`, also applied to
    worker task exceptions in `_collect_parallel_pair_results` — pool-level
    failures keep the discard-and-rerun-in-threads fallback), and the
    run-level `_require_writable_output_root` probe that aborts on an
    unwritable output root before any record is attempted.
  - Durable reporting (A12): `run_binary_bridge_batch` also writes

## Core chemistry seams

- [src/cofkit/reactions.py](../src/cofkit/reactions.py)
  - Reaction templates and linkage profiles, including per-template bridge
    geometry priors (`BridgeGeometryPriors`) that own both the
    realization-time imine/azine bridge constructor and the embedding-time
    motif-origin retraction (`linkage_geometry.derived_origin_retraction_fraction`).
  - This is the first file to touch for a new linkage.
- [src/cofkit/chem/motif_registry.py](../src/cofkit/chem/motif_registry.py)
  - Motif-kind metadata.
  - Add new motif kinds here first.
- [src/cofkit/chem/rdkit.py](../src/cofkit/chem/rdkit.py)
  - Practical SMILES-to-`MonomerSpec` path.
  - Plane-prior seam (A06): `_plane_normal` / `_fit_molecular_plane` derive the
    monomer plane from molecular geometry — the motif-anchor connection plane
    for 3+-motif monomers, a heavy-atom covariance best-fit plane otherwise.
    Non-planar / collinear / degenerate conformers claim no plane
    (`plane_status` metadata, `plane_normal=None`, zero motif frame normals);
    scoring/optimizer consumers must skip plane-dependent terms on zero
    normals and never fabricate a fallback plane. Owners of the fit
    tolerances: `_MONOMER_PLANARITY_RMS_TOLERANCE_ANGSTROM` and
    `_MONOMER_PLANE_COLLINEAR_RATIO`.
  - Conformer selection seam: energy-only by default at the library level;
    `select_conformer_by_motif_shape=True` picks the conformer whose motif
    origins best form a regular planar polygon. Owners of the gates:
    `SHAPE_SELECTION_MIN_MOTIFS` / `SHAPE_SELECTION_VALIDATED_MAX_MOTIFS`
    (binary-bridge builds apply it to 3-connecting monomers only — A/B
    evidence: trigonal improves, tetrahedral worsens; ring-forming applies it
    to all 3+-motif precursors) and `SHAPE_SELECTION_MIN_CONFORMERS`
    (ensemble floor 16).
  - New motif kinds normally need a match handler here.
- [src/cofkit/chem/detector.py](../src/cofkit/chem/detector.py)
  - Lightweight non-RDKit fallback detector.
  - Only some motif kinds are implemented here.
- [src/cofkit/reaction_realization.py](../src/cofkit/reaction_realization.py)
  - Atomistic bond/deletion realization for CIF export.
  - A new linkage is not really complete without this.

## Monomer-library and batch-input seams

- [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py)
  - `MonomerRoleResolver` for autodetection from SMILES.
  - `MonomerRoleResolver.forced_kind_warnings` for non-blocking motif-overlap warnings: a monomer assigned a generic kind (per `_AUTO_DETECT_GENERIC_SUPPRESSION`, currently `aldehyde`) that also builds as the shadowed specific kind (currently `keto_aldehyde`) yields a warning suggesting the specific kind's templates. `infer_record` and `BinaryBridgeLibraryLoader.load_smiles_library` attach these to record metadata as `overlap_warnings`; the CLI echoes them to stderr, and batch runs persist them in the `monomers.jsonl` ledger.
  - Autodetection retains per-kind failure causes (A12): `auto_detect_monomer_candidates_with_causes` returns each skipped kind's `TypeName: message` cause alongside the candidates (expected non-matches and internal failures alike); `infer_record` stores them as `autodetect_failure_causes` in record metadata, and a total detection failure raises a `ValueError` listing the per-kind causes.
  - `BinaryBridgeLibraryLoader` for explicit and autodetected batch libraries.
- [src/cofkit/batch_models.py](../src/cofkit/batch_models.py)
  - Neutral batch-facing data classes.
  - Use these instead of adding new ad hoc summary dictionaries.
  - `BatchPairSummary` carries typed validation accessors (`validation_classification`, `validation_coverage`, `unmeasured_required_checks`) over the serialized `metadata["validation"]` record; `BatchRunSummary.validation_counts` counts validation classifications per run, and `BatchRunSummary.record_failures` aggregates every non-`ok` manifest record as structure id → `TypeName: message` (surfaced in `summary.md` and the console/JSON summaries).
  - `BatchRunSummary` yield accounting (A12): derived `constructed_structures` (status `ok`; the legacy `successful_structures` field remains as an alias — construction success is not a validated yield), `exported_structures` (alias of `cifs_written`), `screened_structures` (validation record attached), and `unvalidated_structures` properties, plus `monomer_records_path` pointing at the durable per-monomer `monomers.jsonl` ledger.

## Topology seams

- [src/cofkit/node_shape.py](../src/cofkit/node_shape.py)
  - Geometric node-shape classification for 4-connecting monomers (square / rectangular / tetrahedral / unknown) from the embedded conformer's connector positions (the four motif-frame origins, fingerprinted by rigid-motion-invariant pairwise angles); the `NODE_SHAPE_*` tolerances are owned here with heuristic provenance. Net node shapes are canonical labels for the curated 4-connected nets (`sql` / `kgm` / `dia`); any other topology id is `unknown`. The graph-automorphism family label (`classify_monomer_graph_shape`) survives as a diagnostic only — it cannot separate square-planar from tetrahedral nodes and feeds no filtering decision.
  - Consumed by `batch.py` topology selection (`_topology_ids_for_pair`, `_topology_unavailable_errors`, decorated-`bex` gating, `_requested_shape_warnings`) and by `planner.py` compatibility checks; unknown shapes are always treated as "no opinion".
  - Shape detection only prunes enumerated topologies. Explicitly requested topologies (`--topology`, `--cofid`, `BatchGenerationConfig.topology_ids` / `single_node_topology_ids`, `NetPlanner` `target_topologies`) are kept and report `shape_warnings` in pair-summary / `NetPlan` / `COFEngine` candidate metadata instead. It can be switched off via `shape_aware_topology_filter` on `BatchGenerationConfig`, the CLI `--no-shape-aware-topology-filter`, `NetPlanner(...)`, and `COFEngineConfig(...)`; enumerations that are passed as `topology_ids` can stay filtered with `BatchGenerationConfig(shape_filter_explicit_topologies=True)`.
  - Dimensionality is a separate, exact check: `BatchGenerationConfig.target_dimensionality` restricts enumerated pools in `_topology_ids_for_pair` (with per-topology reasons from `_topology_unavailable_errors`), and an explicit request whose dimensionality conflicts with the target is unsatisfiable — `BatchStructureGenerator(...)` construction and `NetPlanner.propose` raise `ValueError` (CLI: `error: ...`).
- [src/cofkit/topology_builders.py](../src/cofkit/topology_builders.py)
  - Shared dispatch for supported topology-family builders.
- [src/cofkit/single_node_topologies.py](../src/cofkit/single_node_topologies.py)
  - Space-group-expanded `2D` one-node families. Direction stars are measured from the expanded P1 node sites ("expanded" mode; `fxt` keeps an "explicit" mode), bipartiteness is exact on the quotient graph (parity-expanded BFS over image-shift mod-2 equations), and plane-group operations come from gemmi's 3D space-group tables projected onto the layer (any parseable group, not just `P6/mmm` / `P4/mmm`).
- [src/cofkit/single_node_topologies_3d.py](../src/cofkit/single_node_topologies_3d.py)
  - Supported `3D` one-node families.
- [src/cofkit/indexed_topology_layouts.py](../src/cofkit/indexed_topology_layouts.py)
  - Generic indexed-topology layout reconstruction.
- [src/cofkit/topology_analysis.py](../src/cofkit/topology_analysis.py)
  - Chemistry-facing `two_monomer_*` metadata and lower-level graph diagnostics.
- [src/cofkit/topology_symmetry.py](../src/cofkit/topology_symmetry.py)
  - Generic symmetry expansion for topology analysis.

## Geometry / scoring / validation

- [src/cofkit/geometry.py](../src/cofkit/geometry.py)
  - Shared vector/frame primitives plus the canonical stacking helpers: `measure_layer_z_span` (layer z-span along an explicit axis — prefer `layer_normal_axis`, the `a × b` normal, which is invariant under in-plane periodic images in tilted cells; `LayerSpanReport` carries honest mode/axis provenance), `classify_2d_cell` / `classify_2d_cell_parameters` (the single 2D cell classifiers for built vectors and RCSR parameter cells), `safe_normalize` (explicit per-call-site fallback direction, optional stderr warning), `orthogonal_component`, and `angle_degrees` (caller-chosen degenerate policy). Also owns `covariance_eigenpairs` (Jacobi 3x3 eigensolver returning mean-covariance eigenvalues + eigenvectors, the basis of the `chem/rdkit.py` molecular plane fit), `smallest_covariance_axis` (least-variance-axis convenience wrapper over it), and `planar_arrangement_mismatch` (regular-polygon deviation of a point set — the conformer shape metric used by shape-aware monomer selection). All builder modules delegate here; do not add local copies. Degenerate-vector policy per seam: `stacking.py` warns and falls back, while the thin `_safe_normalize` wrappers in `optimizer.py` / `embedding.py` / `batch.py` raise a descriptive `ValueError` (probe-verified unreachable on live paths; batch loops absorb it per record).
- [src/cofkit/embedding.py](../src/cofkit/embedding.py)
  - Initial periodic placement. `cell_kind` metadata classifies the built cell vectors through `classify_2d_cell`, never the topology id.
- [src/cofkit/optimizer.py](../src/cofkit/optimizer.py)
  - Lightweight continuous refinement.
- [src/cofkit/soft_relax.py](../src/cofkit/soft_relax.py)
  - Experimental in-process clash-repair / strain-relief pass (bond springs + Urey-Bradley 1-3 restraints + ramped soft repulsion) applied to the staged CIF before validation bucketing when `BatchGenerationConfig.soft_relax` / `--soft-relax` is on; not a physical relaxation.
- [src/cofkit/scoring.py](../src/cofkit/scoring.py)
  - Bridge-geometry residual metrics (`bridge_geometry_report`) and the `scoring_metadata` packaging consumed by the optimizer, validator, and candidate ranking. The legacy event-count heuristic score was removed. Plane-dependent terms (planarity residual, normal-misalignment surrogate) are evaluated only for motifs whose monomer carries a plane prior (nonzero frame normal); per-event `plane_prior_coverage` (`both` / `first` / `second` / `none`) records the coverage and `normal_alignment` is null when unevaluated.
- [src/cofkit/vdw.py](../src/cofkit/vdw.py)
  - Shared vdW contact capability: the Bondi radii table (`BONDI_VDW_RADII`, explicit supported-element set; unsupported elements warn on stderr once and use `FALLBACK_VDW_RADIUS` — deliberately not the DREIDING force-field radii), the shared clash criterion (`assess_pair`: `d / (r_i + r_j) < 0.75` or, heavy-heavy only, `d < 2.2 Å`), and `min_periodic_pair_contact` (minimum cross-set contact over every lattice image within the cutoff, including both `(0,0,±1)` c galleries). Consumed by `validation.py` (coarse clash scan) and `stacking.py` (interlayer contact self-check); extend here when a new element needs an exact radius.
- [src/cofkit/validation.py](../src/cofkit/validation.py)
  - `valid` / `warning` / `needs_optimization` / `hard_invalid` / `hard_hard_invalid` / `unvalidated` triage. The contact scan (`_nonbonded_contact_scan`) applies the `vdw.py` criterion with periodic-image-aware 1-2/1-3/1-4 bond-graph exclusions (`_bond_graph_exclusions`), a severe-overlap floor for excluded pairs, and a separate hydrogen-contact metric channel. Bridge-linkage verdicts are derived from the final exported coordinates: `measure_inter_instance_bond_distances` recomputes inter-monomer bond lengths from the atom-site loop plus the explicit `_geom_bond` image codes (never the `_geom_bond_distance` column) and residuals/ratios are checked against the linkage-profile target (or `realized_bridge_bond_distance_windows`); seed `bridge_event_metrics` are recorded as informational `seed_*` metrics only. Every report carries a per-check `coverage` map with statuses `measured` / `no_contacts` (completed search, zero assessable neighbors) / `missing_data` (required check could not run) / `not_applicable` / `skipped`; required-but-unmeasured checks land in `unmeasured_required_checks` and classify the record as `unvalidated` (`is_valid is None`), routed to the `unvalidated` bucket by both `classify_batch_output` and the batch export path.
- [src/cofkit/cif.py](../src/cofkit/cif.py)
  - CIF export, including realized inter-monomer bonds.
- [src/cofkit/decompose_cif.py](../src/cofkit/decompose_cif.py)
  - CIF atom-site and explicit-bond extraction for decomposition without adding ASE.
- [src/cofkit/decompose.py](../src/cofkit/decompose.py)
  - Explicit-bond binary-bridge decomposition into recovered monomers and COFid serialization.
  - Topology auto-detection (`detect_cif_topology`) owns the typed confidence/provenance vocabulary (`TOPOLOGY_CONFIDENCE_*`, `TOPOLOGY_IDENTIFICATION_*`): an embedded COFid comment is promoted to `exact` only when the quotient-graph matcher independently verifies the annotated net (`verified_periodic_graph`); a merely rank/connectivity-compatible comment is an unverified hint (`compatible_unverified`, tie-break confidence `compatible` / `compatible_annotation_unverified`).
  - This logic was adapted from the deCOFpose project at `https://github.com/r-fedorov/deCOFpose`.
- [src/cofkit/decompose_events.py](../src/cofkit/decompose_events.py)
  - Default event/hypothesis decomposition engine: local event detection, bounded per-family hypothesis enumeration (`_MAX_FAMILY_HYPOTHESES`), global validation, and selection. Owns the typed verdict vocabulary (`EVENT_STATUS_*`, including `FAILED_INTERNAL_ERROR` for internal failures — never a chemistry verdict — and `UNSUPPORTED_MIXED_LINKAGE_FAMILY` for multi-family frameworks the one-family COFid contract cannot serialize) and the search-coverage status (`SEARCH_STATUS_*` with explored/theoretical counts; a truncated failed search must not read as exhaustive).
- [src/cofkit/decompose_bond_orders.py](../src/cofkit/decompose_bond_orders.py)
  - Shared graph normalization before event/legacy detection: geometry-constrained, valence-preserving repair of quinoid imine assignments; no new bonds or periodic edge removal. Receives periodic edge multiplicities, checks periodic valences, and retries failing components with parallel-image orders fixed before leaving them unresolved. ReDD-COFFEE regressions and independent CHK boundary labels live in `tests/test_decompose_bond_orders.py` and `tests/fixtures/redd_coffee/`.

## Constant-owner modules

Per the single-owner convention in [AGENTS.md](../AGENTS.md), every
tunable numeric lives in the module that semantically owns it — new
tunables must land here, not at consumer sites:

- [src/cofkit/constants.py](../src/cofkit/constants.py)
  - Cross-cutting assembly defaults (`DEFAULT_MONOLAYER_C_ANGSTROM`, `DEFAULT_LATERAL_SPAN_ANGSTROM`).
- [src/cofkit/ring_geometry.py](../src/cofkit/ring_geometry.py)
  - Ring-template bond lengths for boroxine/triazine ring formation, the ring-arrangement acceptance tolerances (`RingGeometryProfile`), and the exocyclic attachment acceptance criteria (`RING_ATTACHMENT_IDEAL_ANGLE_DEGREES`, derived; `attachment_*` warning/rejection tolerances, heuristic — pending calibration) measured by `ring_attachment_report` and combined in `validate_ring_geometry` (arrangement and attachment verdicts reported separately).
- [src/cofkit/validation.py](../src/cofkit/validation.py)
  - `CoarseValidationThresholds` owns the nine validation thresholds; `cli_analyze.py` argparse defaults reference its fields.
- [src/cofkit/graspa.py](../src/cofkit/graspa.py) and [src/cofkit/lammps.py](../src/cofkit/lammps.py)
  - Config dataclasses own the simulation defaults (cutoffs, EQeq parameters, seeds, temperatures, Ewald tolerances); `cli_calculate.py` argparse defaults reference them.
- [src/cofkit/vdw.py](../src/cofkit/vdw.py)
  - Bondi radii table and clash-assessment constants (pre-existing owner).
- [src/cofkit/node_shape.py](../src/cofkit/node_shape.py)
  - `NODE_SHAPE_*` connector-geometry classification tolerances and the derived ideal tetrahedral angle (`TETRAHEDRAL_NODE_ANGLE_DEGREES`).

## External calculation seams

- [src/cofkit/lammps.py](../src/cofkit/lammps.py)
  - UFF/DREIDING-backed explicit-bond CIF to LAMMPS data/input generation, EQeq charge staging, local minimization/MD orchestration, optional guest-restart atom/force-field merging for MD, and optimized or MD-updated CIF export. The DREIDING path also retypes N/O/F-bound hydrogens as `H__HB` and exports Table V hydrogen-bond terms (`pair_style hybrid/overlay` + `hbond/dreiding/lj` pair coeffs, labeled PairIJ rows) unless `dreiding_hbond=False`. Completed minimizations that miss the force acceptance criteria are exported and marked unconverged (warnings, report/convergence flags, CIF comment) rather than raising; unreadable convergence diagnostics still raise `LammpsExecutionError`.
- [src/cofkit/cif_checks.py](../src/cofkit/cif_checks.py)
  - Structural preconditions shared by the calculation adapters (`read_ordered_structure`) and EQeq charge-assignment validation (`validate_charge_assignment`: atom-count/cell/charge-coverage/net-charge checks plus a verified geometric atom bijection that restores the input labels in the charged CIF). `cif_value_str` is the shared unquoting helper: raw CIF loop tokens keep their delimiters, so any label read from a loop must go through it before being compared with a gemmi-parsed label.
- [src/cofkit/graspa.py](../src/cofkit/graspa.py)
  - EQeq to gRASPA/RASPA2 Widom, single-component isotherm, and mixture workflows; framework mixing-rule generation; simulation.input rendering; result parsing.
- [src/cofkit/hybrid_mdmc.py](../src/cofkit/hybrid_mdmc.py)
  - Cyclic LAMMPS MD plus gRASPA/RASPA2 GCMC workflow. Default framework exchange passes the MD-updated framework CIF between segments; opt-in guest-restart exchange parses final GCMC guest snapshots and passes massive guest atoms into the following LAMMPS MD segment.
  - Partial/failure reporting (A12): `hybrid_mdmc_report.json` is (re)written after every completed cycle (`status: "in_progress"`) and on failure (`status: "failed"`); `HybridMdMcResult.status` / `.failure` (`HybridMdMcFailure`: cycle, failed stage — `lammps_md` / `md_to_gcmc_restart` / `gcmc` / `gcmc_to_lammps_restart` / `cycle_record` — and `TypeName: message` cause) preserve the completed cycles and their warnings. Cycles are dependent simulation states, so a failure stops the run and re-raises after a stderr warning naming the partial report path.
  - Model-contract reporting (A11): every run builds a typed `HybridModelContract` (`_build_hybrid_model_contract`) recorded as `model_contract` on `HybridMdMcResult` and in `hybrid_mdmc_report.json` — kind `approximate_alternating_md_mc` (never a consistently sampled ensemble), the effective interaction settings of each leg (`HybridInteractionSide` md/mc), a `differences` list flagging every divergent aspect (per-run setting mismatches plus the standing structural ones: NVT-MD vs grand-canonical MC, flexible vs rigid framework, element-keyed MC framework LJ rows vs per-type MD assignment, guest presence/rigidity), and `guest_model_diagnostics` (Feynman-Hibbs → classical-LJ conversions and conflicting bundle overrides from `guest_restart.guest_site_model_diagnostics`). A run-level summary warning lists the differing aspects.
- [src/cofkit/guest_restart.py](../src/cofkit/guest_restart.py)
  - GCMC movie/restart snapshot discovery, guest force-field asset synchronization from packaged/bundle RASPA rows into LAMMPS-ready guest sites/templates, binary guest parsing, and zero-mass pseudo-site rejection.
  - Site-to-component lookup is deduplicated per template (`_site_component_candidates`) so multisite guests with repeated site labels (CO2, SO2) parse; a label is ambiguous only when genuinely shared between components.
  - Guest-model conversion/conflict diagnostics (A11): `LammpsGuestSite` keeps the source row's `interaction` kind (`lennard-jones` / `feynman-hibbs-lennard-jones`) and `override_conflicts` (disagreements between a bundle's `lammps.masses` / `charges` / `pair_coeff_rows` overrides and the RASPA-side pseudo-atom/mixing-rule rows); `guest_site_model_diagnostics` turns them into explicit warning strings that ride on the restart state warnings, so a Feynman-Hibbs row staged as classical LJ or a divergent MD-side override is always reported, never silent.
  - Population accounting is explicit in both handoff directions: `LammpsGuestRestartState` carries `skipped_molecules` / `skipped_unknown_site_counts` / `n_skipped_ambiguous_atoms`, and the MD→MC conversion warns with a recovered-vs-input atom summary whenever counts differ. Framework atoms embedded in gRASPA movies are excluded by label and counted, never silently dropped.
  - gRASPA snapshot production is verified against the backend: movies are written during production at `cycle % MoviesEvery == 0` (0-based; default 5000). `GraspaIsothermSettings` / `GraspaMixtureSettings` `movies_every` renders an explicit `MoviesEvery` line (gRASPA only), and hybrid `guest_restart` mode sets it to `max(1, production_cycles - 1)` to guarantee a final-cycle snapshot.
  - Real-backend smoke coverage lives in `tests/engines/test_graspa_engine.py`, gated on `COFKIT_TEST_GRASPA` / `COFKIT_TEST_EQEQ` / `COFKIT_TEST_LMP` like the LAMMPS engine tests.
- [src/cofkit/guest_bundles.py](../src/cofkit/guest_bundles.py)
  - Shared external guest parameter-bundle contract for gRASPA/RASPA2 workflows. Bundles add force-field provenance and compatibility metadata, RASPA molecule definitions, pseudo-atom rows, mixing-rule rows, aliases, rotatability, and a required non-empty `lammps` section for synchronized hybrid MD/MC guest restart.
- [src/cofkit/guest_forcefields.py](../src/cofkit/guest_forcefields.py)
  - Strict loader and public metadata model for packaged guest native parameter families and allowed framework-force-field combinations.

## Typical execution paths

### Single pair from CLI

1. `cofkit build single-pair` in [src/cofkit/cli_build.py](../src/cofkit/cli_build.py)
2. motif-kind autodetection via [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py) if explicit kinds are omitted
3. monomer construction via [src/cofkit/chem/rdkit.py](../src/cofkit/chem/rdkit.py)
4. candidate generation via [src/cofkit/batch.py](../src/cofkit/batch.py)
5. CIF export and validation bucketing via [src/cofkit/cif.py](../src/cofkit/cif.py) and [src/cofkit/validation.py](../src/cofkit/validation.py)

### Full batch from CLI

1. `cofkit build batch-binary-bridge` or `cofkit build batch-all-binary-bridges`
2. library resolution via [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py)
3. pair enumeration and topology selection in [src/cofkit/batch.py](../src/cofkit/batch.py)
4. topology-family dispatch via [src/cofkit/topology_builders.py](../src/cofkit/topology_builders.py)
5. validation-aware CIF writing into `cifs/valid`, `cifs/warning`, `cifs/needs_optimization`, `cifs/unvalidated`, or `cifs/hard_invalid`; `hard_hard_invalid` structures stay manifest-only

### CIF decomposition from CLI

1. `cofkit analyze decompose <CIF>` in [src/cofkit/cli_analyze.py](../src/cofkit/cli_analyze.py), with optional `--topology <TOKEN>` when the topology is known or auto-detection is ambiguous
2. atom and bond extraction via [src/cofkit/decompose_cif.py](../src/cofkit/decompose_cif.py)
3. linkage-specific cutting and monomer repair via [src/cofkit/decompose.py](../src/cofkit/decompose.py)
4. COFid serialization through [src/cofkit/cofid.py](../src/cofkit/cofid.py)

### gRASPA/RASPA2 calculation from CLI

1. `cofkit calculate graspa-widom`, `graspa-isotherm`, or `graspa-mixture` in [src/cofkit/cli_calculate.py](../src/cofkit/cli_calculate.py)
2. optional external guest loading and alias canonicalization via [src/cofkit/guest_bundles.py](../src/cofkit/guest_bundles.py) and [src/cofkit/graspa.py](../src/cofkit/graspa.py)
3. EQeq framework charge staging in [src/cofkit/graspa.py](../src/cofkit/graspa.py)
4. backend-specific `simulation.input` and force-field asset materialization in [src/cofkit/graspa.py](../src/cofkit/graspa.py)
5. parsed JSON/CSV report writing from backend `Output/**/*.data`

### Hybrid MD/MC calculation from CLI

1. `cofkit calculate hybrid-mdmc` in [src/cofkit/cli_calculate.py](../src/cofkit/cli_calculate.py)
2. cycle orchestration in [src/cofkit/hybrid_mdmc.py](../src/cofkit/hybrid_mdmc.py)
3. per-cycle LAMMPS MD framework update through [src/cofkit/lammps.py](../src/cofkit/lammps.py)
4. per-cycle pure-component or mixture GCMC through [src/cofkit/graspa.py](../src/cofkit/graspa.py)
5. in `framework` mode, the next cycle starts from the LAMMPS MD output CIF only
6. in `guest_restart` mode, [src/cofkit/guest_restart.py](../src/cofkit/guest_restart.py) parses the final GCMC `Movies/System_0/result_*.data` snapshot and [src/cofkit/lammps.py](../src/cofkit/lammps.py) injects those guests into the next MD data file

## Current architectural constraints

These are important before editing:

- The practical generation path is still centered on one binary-bridge template per run.
- Batch autodetection still requires one resolved motif kind per monomer.
- RDKit is the practical chemistry path; the fallback detector is narrower.
- Non-binary / ring-forming chemistries are not yet fully wired into topology-guided generation.
- Validation thresholds are still globally configured, not fully per-linkage.

## Where to add tests

- [tests/test_cli.py](../tests/test_cli.py)
  - CLI coverage.
- [tests/test_batch.py](../tests/test_batch.py)
  - Batch and single-pair generation behavior.
- [tests/test_rdkit_monomer.py](../tests/test_rdkit_monomer.py)
  - RDKit motif detection and monomer building.
- [tests/test_reaction_realization.py](../tests/test_reaction_realization.py)
  - Atomistic reaction realization.
- [tests/test_decompose_cif.py](../tests/test_decompose_cif.py)
  - No-ASE CIF extraction and explicit-bond preservation.
- [tests/test_decompose.py](../tests/test_decompose.py)
  - CIF-to-COFid decomposition and generated hcb round trips across buildable binary-bridge linkages.
- [tests/test_lammps.py](../tests/test_lammps.py)
  - LAMMPS data/input generation, force-field parameter paths, optimization/MD orchestration, and CIF preservation behavior.
- [tests/test_soft_relax.py](../tests/test_soft_relax.py)
  - Soft-relax clash repair on synthetic CIFs plus the `--soft-relax` build-pipeline wiring.
- [tests/test_graspa.py](../tests/test_graspa.py)
  - EQeq/gRASPA/RASPA2 workflow staging, CLI parsing, guest bundles, force-field asset generation, and parser behavior.
- [tests/test_hybrid_mdmc.py](../tests/test_hybrid_mdmc.py)
  - Hybrid MD/MC cycle orchestration, framework snapshot handoff behavior, and binary guest-restart propagation.
- [tests/test_guest_restart.py](../tests/test_guest_restart.py)
  - Guest restart force-field synchronization, GCMC movie snapshot parsing, and unsupported massless pseudo-site handling.
- [tests/test_core.py](../tests/test_core.py)
  - End-to-end project / template behavior.

## Practical advice for later sessions

- Start from the registry file, not the executor.
- If a new chemistry exports CIFs incorrectly, inspect [src/cofkit/reaction_realization.py](../src/cofkit/reaction_realization.py) before changing scoring or topology code.
- If autodetected libraries behave strangely, inspect [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py) before [src/cofkit/batch.py](../src/cofkit/batch.py).
- If a topology is present in metadata but never generated, inspect:
  - [src/cofkit/topology_analysis.py](../src/cofkit/topology_analysis.py)
  - [src/cofkit/topology_builders.py](../src/cofkit/topology_builders.py)
  - the default-selector logic in [src/cofkit/batch.py](../src/cofkit/batch.py)
