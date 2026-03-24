---
name: cofkit-navigator
description: Translate natural-language requests about using cofkit into the correct CLI or Python workflow. Use when a user wants to generate or screen COFs from SMILES or monomer-library text files, inspect supported reaction templates, classify finished batch outputs, rebuild the detector-scanned default library, or choose between cofkit's CLI, BatchStructureGenerator, and COFEngine surfaces.
license: MIT
---

# Cofkit Navigator

## Overview

Use this skill to operate `cofkit` for end-user tasks. Translate the request into the narrowest real `cofkit` workflow, run it, then summarize the generated artifacts instead of answering from memory.

Prefer the installed `cofkit` CLI for human-facing requests. If you are working inside the repository without an editable install, use `PYTHONPATH=src .venv/bin/python -m cofkit.cli ...` as the equivalent launcher.

## Choose The Surface

- `cofkit list-templates`: capability discovery. Use this first when the user asks what chemistries are available or when you need to confirm whether a template supports topology-guided pair generation.
- `cofkit single-pair`: one pair of monomers from SMILES. This is the default for requests like "build me a COF from these two monomers".
- `cofkit batch-binary-bridge`: one binary-bridge template over a monomer-library directory.
- `cofkit batch-all-binary-bridges`: one monomer-library directory across every currently supported binary-bridge template that the input can satisfy.
- `cofkit classify-output`: post-process an existing batch output tree into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid`.
- `cofkit build-default-library`: rebuild the detector-scanned grouped example library from `examples/batch_test_monomers`.
- `BatchStructureGenerator`: use this when code already has `MonomerSpec` objects or when you need library/template discovery from Python.
- `COFEngine`: use this when the user already knows the target topology and wants the direct project-style API.

## Operating Rules

- Treat the current practical generation path as binary-bridge-first. Not every registered reaction template can run pair generation; confirm with `list-templates` and inspect `supports_pair_generation` plus `supports_topology_guided_generation`.
- `single-pair` auto-detects motif kinds by default. If the user already knows the motif kinds, pass them explicitly to reduce ambiguity.
- Use `--auto-detect-libraries` only for raw generic `.txt` SMILES libraries. Do not use it on `examples/default_monomers_library`; that directory is already grouped detector output.
- Prefer the generic commands over legacy shortcuts. `batch-imine` is only an imine alias; default to `batch-binary-bridge --template-id imine_bridge` unless the user explicitly asks for the shortcut.
- Read the generated artifacts before reporting success. Single-pair runs write `summary.json`. Batch runs write `manifest.jsonl` and `summary.md`. `batch-all-binary-bridges` also writes `combined_summary.json`. Library rebuilding writes `README.md`, `registry.jsonl`, and `failures.jsonl`.
- Keep current constraints visible: `stacking_mode` stays `"disabled"`, RDKit is the practical path for SMILES-driven workflows, and the shipped example library currently exposes `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge`.

## Minimal Input Checklist

Ask only for the missing field that blocks execution:

- one-pair generation: template ID or desired chemistry, two SMILES strings, optional topology or dimensionality, output directory
- library screening: input directory, template ID or "all supported", whether the libraries are already grouped, output directory
- output triage: existing batch output directory and optional classification output directory

## References

- For exact intent-to-command routing and output interpretation, read `references/intent-routing.md`.
- For broader project scope and examples, use `../../README.md` and `../../USER_MANUAL.md`.
- For code-editing or chemistry-extension tasks rather than end-user operation, use `../../agent-docs/AGENT_CODEBASE_MAP.md` and `../../agent-docs/ADDING_LINKAGES_AND_MONOMERS.md`.
