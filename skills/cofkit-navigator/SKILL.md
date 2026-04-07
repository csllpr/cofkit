---
name: cofkit-navigator
description: Translate natural-language requests about using cofkit into the correct CLI or Python workflow. Use when a user wants to generate or screen COFs from SMILES or monomer-library text files, inspect supported reaction templates, classify finished batch outputs, run Zeo++ pore analysis, run LAMMPS optimization, run the EQeq-to-gRASPA Widom wrapper, rebuild the detector-scanned default library, or choose between cofkit's CLI, BatchStructureGenerator, and COFEngine surfaces.
license: MIT
---

# Cofkit Navigator

## Overview

Use this skill to operate `cofkit` for end-user tasks. Translate the request into the narrowest real `cofkit` workflow, run it, then summarize the generated artifacts instead of answering from memory.

Prefer the installed `cofkit` CLI for human-facing requests. If you are working inside the repository without an editable install, use `PYTHONPATH=src python -m cofkit.cli ...` or the current project interpreter equivalent instead of assuming a checked-in `.venv`. That repo-local launcher still requires the project dependencies to be installed in the interpreter you use.

## Choose The Surface

- `cofkit build list-templates`: capability discovery. Use this first when the user asks what chemistries are available or when you need to confirm whether a template supports topology-guided pair generation.
- `cofkit build single-pair`: one pair of monomers from SMILES. This is the default for requests like "build me a COF from these two monomers". If no `--topology` is supplied, the current CLI default is to enumerate all applicable topologies for the pair.
- `cofkit build batch-binary-bridge`: one binary-bridge template over a monomer-library directory. If no `--topology` is supplied, the current CLI default is to enumerate all applicable topologies for each compatible pair.
- `cofkit build batch-all-binary-bridges`: one monomer-library directory across every currently supported binary-bridge template that the input can satisfy.
- `cofkit analyze classify-output`: post-process an existing batch output tree into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid`.
- `cofkit analyze zeopp`: run Zeo++ pore-property analysis on one CIF, including the default point-probe baseline and optional repeated accessibility-aware probe scans.
- `cofkit calculate lammps-optimize`: run the current UFF-backed LAMMPS local optimization on one explicit-bond `P1` CIF.
- `cofkit calculate graspa-widom`: run the current staged `EQeq -> gRASPA` Widom insertion workflow on one CIF and parse the selected packaged Widom component screen.
- `cofkit build default-library`: rebuild the detector-scanned grouped example library from `examples/batch_test_monomers`.
- `BatchStructureGenerator`: use this when code already has `MonomerSpec` objects or when you need library/template discovery from Python.
- `COFEngine`: use this when the user already knows the target topology and wants the direct project-style API.

## Operating Rules

- Treat the current practical generation path as binary-bridge-first. Not every registered reaction template can run pair generation; confirm with `cofkit build list-templates` and inspect `supports_pair_generation` plus `supports_topology_guided_generation`.
- `single-pair` auto-detects motif kinds by default. If the user already knows the motif kinds, pass them explicitly to reduce ambiguity.
- Unless the user asks to restrict topology selection, leaving `--topology` unset keeps `--all-topologies` enabled by default. In the current CLI that means both `single-pair` and `batch-binary-bridge` enumerate every applicable topology per monomer pair; add `--no-all-topologies` to keep only the best topology per pair.
- Use `--auto-detect-libraries` only for raw generic `.txt` SMILES libraries. Do not use it on `examples/default_monomers_library`; that directory is already grouped detector output.
- Prefer the grouped generic commands over legacy shortcuts. Deprecated flat aliases such as `cofkit batch-imine` still work, but default to `cofkit build batch-binary-bridge --template-id imine_bridge`.
- Read the generated artifacts before reporting success. Single-pair runs write `summary.json`. Batch runs write `manifest.jsonl` and `summary.md`. `batch-all-binary-bridges` also writes `combined_summary.json`. Classification writes `classification_manifest.jsonl`. Zeo++ writes `zeopp_report.json` plus raw outputs and logs. LAMMPS writes `lammps_report.json`, an optimized CIF, the generated data/input files, logs, and a dump trajectory. `graspa-widom` writes `graspa_widom_report.json`, staged `eqeq/` and `widom/` directories, gRASPA logs, raw `Output/*.data`, and `widom/Output/results.csv`. Library rebuilding writes `README.md`, `registry.jsonl`, and `failures.jsonl`.
- Keep current constraints visible: `stacking_mode` stays `"disabled"`, RDKit-backed monomer construction is part of the normal install and is the required path for SMILES-driven workflows, the shipped example library currently exposes `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge`, and the public build workflow is still binary-bridge-first.
- Do not route users to the internal sulfur-enabled benzothiazole conversion prototype through the public CLI. That path remains internal-only until its geometry is more reliable.
- `cofkit analyze zeopp` requires the Zeo++ `network` binary through `COFKIT_ZEOPP_PATH` or `--zeopp-path`. Without any `--probe-radius`, it still writes a point-probe baseline. Probe scans can fail partially; inspect `zeopp_report.json` instead of assuming one success count captures the whole run.
- `cofkit calculate lammps-optimize` requires an explicit-bond `P1` CIF with `_geom_bond_*` and `_ccdc_geom_bond_type`. The public force field is `UFF` only. `COFKIT_LMP_PATH` or `--lmp-path` selects the executable, and if `OMP_NUM_THREADS` is unset, cofkit defaults LAMMPS to half the machine core count.
- `cofkit calculate graspa-widom` requires one CIF, an EQeq executable through `COFKIT_EQEQ_PATH` or `--eqeq-path`, and a gRASPA executable through `COFKIT_GRASPA_PATH` or `--graspa-path`. The current wrapper runs EQeq directly on the CIF, materializes `widom/framework.cif`, derives `UnitCells` from the cell lengths plus the larger cutoff, and expects gRASPA output under `widom/Output/*.data`. On many installations the gRASPA binary is named `nvc_main.x`.
- The current public gRASPA wrapper is intentionally narrow: it runs one bundled Widom-template family, activates packaged probes on demand through repeated `--component NAME` flags or `--all-components`, derives total production cycles from `--widom-moves-per-component` unless `--production-cycles` is set, and parses the resulting Widom / Henry summaries. Packaged probes currently cover `TIP4P`, `CO2`, `H2`, `N2`, `SO2`, `Xe`, and `Kr`. Non-finite uncertainty fields can appear in raw gRASPA output; check `graspa_widom_report.json` warnings and expect those specific values to be recorded as `null`.

## Minimal Input Checklist

Ask only for the missing field that blocks execution:

- one-pair generation: template ID or desired chemistry, two SMILES strings, optional topology or dimensionality, output directory. If topology is omitted, the current default is to evaluate all applicable topologies for the pair.
- library screening: input directory, template ID or "all supported", whether the libraries are already grouped, optional topology or dimensionality, output directory. If topology is omitted, the current default is to evaluate all applicable topologies for each compatible pair.
- output triage: existing batch output directory and optional classification output directory
- Zeo++ pore analysis: CIF path, optional output directory, optional probe radii, and Zeo++ binary path only if `COFKIT_ZEOPP_PATH` is not already configured
- LAMMPS optimization: explicit-bond `P1` CIF path, optional output directory, optional optimization tuning flags, and LAMMPS binary path only if `COFKIT_LMP_PATH` is not already configured
- gRASPA Widom insertion: CIF path, optional output directory, optional EQeq / Widom tuning flags, and EQeq / gRASPA binary paths only if `COFKIT_EQEQ_PATH` and `COFKIT_GRASPA_PATH` are not already configured

## References

- For exact intent-to-command routing and output interpretation, read `references/intent-routing.md`.
- For broader project scope and examples, use `../../README.md` and `../../USER_MANUAL.md`.
- For code-editing or chemistry-extension tasks rather than end-user operation, use `../../agent-docs/AGENT_CODEBASE_MAP.md` and `../../agent-docs/ADDING_LINKAGES_AND_MONOMERS.md`.
