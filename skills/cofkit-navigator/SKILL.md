---
name: cofkit-navigator
description: Translate natural-language requests about using cofkit into the correct CLI or Python workflow. Use when a user wants to generate or screen COFs from SMILES or monomer-library text files, inspect supported reaction templates, classify finished batch outputs, run Zeo++ pore analysis, run LAMMPS optimization, run the EQeq-to-gRASPA/RASPA2 Widom, adsorption-isotherm, or mixture wrappers, rebuild the detector-scanned default library, or choose between cofkit's CLI, BatchStructureGenerator, and COFEngine surfaces.
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
- `cofkit calculate lammps-optimize`: run the current UFF/DREIDING-backed LAMMPS local optimization on one explicit-bond `P1` CIF.
- `cofkit calculate graspa-widom`: run the current staged `EQeq -> gRASPA/RASPA2` Widom insertion workflow on one CIF and parse the selected packaged Widom component screen.
- `cofkit calculate graspa-isotherm`: run the current staged `EQeq -> gRASPA/RASPA2` single-component adsorption workflow on one CIF and parse one or more pressure-point loading summaries.
- `cofkit calculate graspa-mixture`: run the current staged `EQeq -> gRASPA/RASPA2` multi-component adsorption workflow on one CIF and parse component loadings plus pairwise selectivities.
- `cofkit build default-library`: rebuild the detector-scanned grouped example library from `examples/batch_test_monomers`.
- `BatchStructureGenerator`: use this when code already has `MonomerSpec` objects or when you need library/template discovery from Python.
- `COFEngine`: use this when the user already knows the target topology and wants the direct project-style API.

## Operating Rules

- Treat the current practical generation path as binary-bridge-first. Not every registered reaction template can run pair generation; confirm with `cofkit build list-templates` and inspect `supports_pair_generation` plus `supports_topology_guided_generation`.
- `single-pair` auto-detects motif kinds by default. If the user already knows the motif kinds, pass them explicitly to reduce ambiguity.
- Unless the user asks to restrict topology selection, leaving `--topology` unset keeps `--all-topologies` enabled by default. In the current CLI that means both `single-pair` and `batch-binary-bridge` enumerate every applicable topology per monomer pair; add `--no-all-topologies` to keep only the best topology per pair.
- Use `--auto-detect-libraries` only for raw generic `.txt` SMILES libraries. Do not use it on `examples/default_monomers_library`; that directory is already grouped detector output.
- Prefer the grouped generic commands over legacy shortcuts. Deprecated flat aliases such as `cofkit batch-imine` still work, but default to `cofkit build batch-binary-bridge --template-id imine_bridge`.
- Read the generated artifacts before reporting success. Single-pair runs write `summary.json`. Batch runs write `manifest.jsonl` and `summary.md`. `batch-all-binary-bridges` also writes `combined_summary.json`. Classification writes `classification_manifest.jsonl`. Zeo++ writes `zeopp_report.json` plus raw outputs and logs. LAMMPS writes `lammps_report.json`, an optimized CIF, the generated data/input files, logs, and a dump trajectory. `graspa-widom` writes `graspa_widom_report.json`, staged `eqeq/` and `widom/` directories, backend logs, raw `Output/**/*.data`, and `widom/Output/results.csv`. `graspa-isotherm` writes `graspa_isotherm_report.json`, staged `eqeq/` and `isotherm/pressure_*/` directories, backend logs, raw `Output/**/*.data` per pressure point, and `isotherm/results.csv`. `graspa-mixture` writes `graspa_mixture_report.json`, staged `eqeq/` and `mixture/pressure_*/` directories, backend logs, raw `Output/**/*.data` per pressure point, plus `mixture/component_results.csv` and `mixture/selectivity_results.csv`. Library rebuilding writes `README.md`, `registry.jsonl`, and `failures.jsonl`.
- Keep current constraints visible: `stacking_mode` stays `"disabled"`, RDKit-backed monomer construction is part of the normal install and is the required path for SMILES-driven workflows, the shipped example library currently exposes `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge`, and the public build workflow is still binary-bridge-first.
- Do not route users to the internal sulfur-enabled benzothiazole conversion prototype through the public CLI. That path remains internal-only until its geometry is more reliable.
- `cofkit analyze zeopp` requires the Zeo++ `network` binary through `COFKIT_ZEOPP_PATH` or `--zeopp-path`. Without any `--probe-radius`, it still writes a point-probe baseline. Probe scans can fail partially; inspect `zeopp_report.json` instead of assuming one success count captures the whole run.
- `cofkit calculate lammps-optimize` requires an explicit-bond `P1` CIF with `_geom_bond_*` and `_ccdc_geom_bond_type`. Both `DREIDING` and `UFF` are available. Prefer `--forcefield dreiding` for real optimization work unless the user explicitly asks for UFF or is reproducing older results. Treat `UFF` as experimentally supported for now. `COFKIT_LMP_PATH` or `--lmp-path` selects the executable, and if `OMP_NUM_THREADS` is unset, cofkit defaults LAMMPS to half the machine core count.
- `cofkit calculate graspa-widom` requires one CIF, an EQeq executable through `COFKIT_EQEQ_PATH` or `--eqeq-path`, and either a gRASPA executable through `COFKIT_GRASPA_PATH` / `--graspa-path` or a RASPA2 executable through `COFKIT_RASPA2_PATH` / `--raspa2-path` with `--backend raspa2`. The current wrapper runs EQeq directly on the CIF, materializes `widom/framework.cif`, derives `UnitCells` from the cell lengths plus the larger cutoff, and expects backend output under `widom/Output/**/*.data`.
- `cofkit calculate graspa-isotherm` requires one CIF, an EQeq executable through `COFKIT_EQEQ_PATH` or `--eqeq-path`, a selected Monte Carlo backend executable, exactly one packaged `--component NAME`, and one or more `--pressure` points. The current wrapper runs EQeq directly on the CIF, materializes `isotherm/framework.cif`, derives shared `UnitCells` from the cell lengths plus the larger cutoff, stages one backend run per pressure under `isotherm/pressure_*/`, and expects output under each pressure-point `Output/**/*.data`. `--fugacity-coefficient` accepts either a positive float or `PR-EOS`; RASPA2 uses internal fugacity calculation when `PR-EOS` is requested.
- `cofkit calculate graspa-mixture` requires one CIF, an EQeq executable through `COFKIT_EQEQ_PATH` or `--eqeq-path`, a selected Monte Carlo backend executable, at least two packaged `--component NAME:FRACTION` values, and one or more `--pressure` points. The current wrapper stages one backend run per pressure under `mixture/pressure_*/`, expects output under each pressure-point `Output/**/*.data`, and computes pairwise selectivities from parsed component loadings.
- The current public gRASPA/RASPA2 surface is intentionally narrow: `graspa-widom` runs one bundled Widom-template family and parses Widom / Henry summaries, `graspa-isotherm` runs one packaged adsorbate at a time over explicit pressure points and parses absolute-loading / heat-of-adsorption summaries, and `graspa-mixture` runs packaged multi-component feeds and parses component loadings plus selectivity. Packaged probes currently cover `TIP4P`, `CO2`, `H2`, `N2`, `SO2`, `Xe`, and `Kr`. Prefer `--forcefield dreiding` for real adsorption work unless the user explicitly asks for UFF or needs a comparison run; treat `UFF` as experimentally supported for now. Non-finite fields can appear in raw backend output; check the `graspa_*_report.json` warnings and expect those values to be recorded as `null`. RASPA2 is parsed from its own output schema: Widom / Henry and mol/kg loading are supported, but RASPA2 `graspa-isotherm` reports no gRASPA-style `g/L` loading, and RASPA2 `graspa-mixture` currently records component `g/L` loading and heat fields as `null`. Ratios formed from separate pure-component `graspa-isotherm` runs are loading ratios, not mixture selectivities.

## Minimal Input Checklist

Ask only for the missing field that blocks execution:

- one-pair generation: template ID or desired chemistry, two SMILES strings, optional topology or dimensionality, output directory. If topology is omitted, the current default is to evaluate all applicable topologies for the pair.
- library screening: input directory, template ID or "all supported", whether the libraries are already grouped, optional topology or dimensionality, output directory. If topology is omitted, the current default is to evaluate all applicable topologies for each compatible pair.
- output triage: existing batch output directory and optional classification output directory
- Zeo++ pore analysis: CIF path, optional output directory, optional probe radii, and Zeo++ binary path only if `COFKIT_ZEOPP_PATH` is not already configured
- LAMMPS optimization: explicit-bond `P1` CIF path, optional output directory, optional optimization tuning flags, and LAMMPS binary path only if `COFKIT_LMP_PATH` is not already configured
- gRASPA/RASPA2 Widom insertion: CIF path, optional output directory, optional EQeq / Widom tuning flags, backend selection, and binary paths only if `COFKIT_EQEQ_PATH` plus the selected backend env var are not already configured
- gRASPA/RASPA2 adsorption isotherm: CIF path, one packaged component name, one or more pressure points in Pa, optional output directory, optional EQeq / adsorption tuning flags, backend selection, and binary paths only if `COFKIT_EQEQ_PATH` plus the selected backend env var are not already configured
- gRASPA/RASPA2 mixture adsorption/selectivity: CIF path, at least two packaged `NAME:FRACTION` components, one or more pressure points in Pa, optional output directory, optional EQeq / mixture tuning flags, backend selection, and binary paths only if `COFKIT_EQEQ_PATH` plus the selected backend env var are not already configured

Default routing preference for forcefields:

- for `lammps-optimize`, `graspa-widom`, `graspa-isotherm`, and `graspa-mixture`, prefer `--forcefield dreiding` unless the user explicitly requests `UFF` or asks to reproduce legacy `UFF` outputs
- describe `UFF` as experimentally supported rather than the recommended production path

## References

- For exact intent-to-command routing and output interpretation, read `references/intent-routing.md`.
- For broader project scope and examples, use `../../README.md` and `../../USER_MANUAL.md`.
- For code-editing or chemistry-extension tasks rather than end-user operation, use `../../agent-docs/AGENT_CODEBASE_MAP.md` and `../../agent-docs/ADDING_LINKAGES_AND_MONOMERS.md`.
