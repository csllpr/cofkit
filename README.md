# cofkit

A reaction-aware periodic assembly toolkit for covalent organic frameworks (COFs).

## Overview

`cofkit` brings together chemistry-aware monomer definitions, reaction templates, topology selection, periodic assembly, scoring, validation, and CIF export in one package.

Start here based on what you need:

- use the installed `cofkit` CLI for practical single-pair and batch workflows
- use [USER_MANUAL.md](USER_MANUAL.md) for broader Python and workflow examples
- use [skills/cofkit-navigator/SKILL.md](skills/cofkit-navigator/SKILL.md) if you are operating the repo through Codex or another compatible agent
- use [docs/README.md](docs/README.md) for technical background extracted from this README

For current capabilities and limits, see [docs/CURRENT_SCOPE.md](docs/CURRENT_SCOPE.md).

## Installation and dependencies

Core install:

```bash
python3 -m pip install -e .
```

After editable install, the unified CLI is available as:

```bash
cofkit --help
```

Core `cofkit` has no mandatory third-party runtime dependencies beyond Python `3.10+`.

Optional tools you may want in your environment:

- `gemmi` for rebuilding topology metadata with generic space-group expansion, CIF-backed coarse validation, and broader `2D` / `3D` two-monomer compatibility scans
- `RDKit` for SMILES-based monomer construction and the practical batch-generation workflows
- `pytest` to run the test suite

The bundled topology repository under [`src/cofkit/data/topologies`](src/cofkit/data/topologies) is sufficient for normal use. External RCSR archives and topology environment variables are optional advanced inputs, not required setup steps.

## CLI

After editable install, the main user-facing interface is the `cofkit` command.

Inspect the available commands with:

```bash
cofkit --help
cofkit list-templates
```

The most useful subcommands are:

- `cofkit single-pair`
- `cofkit batch-binary-bridge`
- `cofkit batch-all-binary-bridges`
- `cofkit classify-output`
- `cofkit build-default-library`

### Single pair

Direct single-pair generation with motif autodetection and default topology enumeration:

```bash
cofkit single-pair \
  --template-id imine_bridge \
  --first-smiles 'C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N' \
  --second-smiles 'C1=C(C=C(C=C1C=O)C=O)C=O' \
  --first-id tapb \
  --second-id tfb \
  --output-dir out/cli_single_pair
```

That command autodetects the monomer roles (`amine` / `aldehyde` for this example), evaluates the current default topology pool for the pair, writes `summary.json`, and exports CIFs into validation buckets under `cifs/valid`, `cifs/warning`, or `cifs/invalid`.

If you want to force one topology instead of using the default topology pool:

```bash
cofkit single-pair \
  --template-id imine_bridge \
  --first-smiles 'C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N' \
  --second-smiles 'C1=C(C=C(C=C1C=O)C=O)C=O' \
  --first-id tapb \
  --second-id tfb \
  --topology hcb \
  --output-dir out/cli_single_pair_hcb
```

### One linkage over a library

Run one binary-bridge linkage over a monomer-library directory:

```bash
cofkit batch-binary-bridge \
  --template-id imine_bridge \
  --input-dir examples/default_monomers_library \
  --output-dir out/batch_imine_from_cli \
  --max-workers 8
```

### All available binary-bridge linkages

Run all currently available binary-bridge templates over one library directory:

```bash
cofkit batch-all-binary-bridges \
  --input-dir examples/default_monomers_library \
  --output-dir out/available_binary_bridge_batches \
  --max-workers 8
```

That CLI path reproduces the current per-template summary results for the supported `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge` batches in the shipped example library.

### Classify a finished output tree

```bash
cofkit classify-output \
  out/full_cif_generation_default_selector_20260320 \
  --output-dir out/full_cif_generation_default_selector_20260320_coarse_validation_triage
```

### Rebuild the detector-scanned example library

```bash
cofkit build-default-library
```

## Skill integration

If you are using Codex inside this repository, load [skills/cofkit-navigator/SKILL.md](skills/cofkit-navigator/SKILL.md). The skill routes natural-language requests into the narrowest supported `cofkit` workflow and is intended for:

- single-pair generation from two SMILES strings
- library-scale binary-bridge screening
- classification of finished output trees
- rebuilding the detector-scanned default monomer library
- choosing between the CLI, `BatchStructureGenerator`, and `COFEngine`

Typical prompts:

- "Build a single-pair imine COF from these two SMILES."
- "Run all supported binary-bridge batches over this library directory."
- "Classify this finished output tree into validation buckets."
- "Should I use the CLI, `BatchStructureGenerator`, or `COFEngine` for this task?"

The skill prefers the installed `cofkit` CLI. If the package is not installed in editable mode, it falls back to an equivalent `PYTHONPATH=src ...` launcher.

## Additional documentation

- [USER_MANUAL.md](USER_MANUAL.md) for end-user Python and workflow examples
- [docs/README.md](docs/README.md) for scope, wrapper scripts, topology notes, pipeline details, and validation thresholds
- [agent-docs/README.md](agent-docs/README.md) for development and agent-facing navigation guides
- [ARCHITECTURE.md](ARCHITECTURE.md) and [CHANGELOG.md](CHANGELOG.md) for broader project context
