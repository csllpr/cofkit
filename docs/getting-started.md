# Getting Started

## Install

Use the lockfile-backed `uv` environment:

```bash
uv sync --locked
uv run cofkit --help
```

For local development and verification:

```bash
uv sync --locked --extra dev
uv run pytest -q
```

If you explicitly need an editable install inside an existing Python environment, `python3 -m pip install -e .` still works. The repository-standard path is `uv run ...` because it keeps helper executables aligned with the repo-managed environment.

## CLI Shape

The installable CLI is grouped by task:

```bash
cofkit build --help
cofkit analyze --help
cofkit calculate --help
cofkit validate --help
cofkit build list-templates
```

Legacy flat aliases such as `cofkit single-pair`, `cofkit batch-imine`, and `cofkit classify-output` still work for now, but they emit deprecation warnings.

## First Single-Pair Build

```bash
cofkit build single-pair \
  --template-id imine_bridge \
  --first-smiles 'C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N' \
  --second-smiles 'C1=C(C=C(C=C1C=O)C=O)C=O' \
  --first-id tapb \
  --second-id tfb \
  --output-dir out/cli_single_pair
```

The command autodetects motif roles, evaluates all applicable default topologies, writes `summary.json`, and exports CIFs under `cifs/valid`, `cifs/warning`, `cifs/needs_optimization`, or `cifs/hard_invalid`.

## Environment Files

For CLI convenience, `cofkit` automatically loads the nearest `.env` file by searching upward from the current working directory. Values already present in the shell environment are left unchanged, so explicit exports still take precedence.

Common optional variables:

```bash
COFKIT_ZEOPP_PATH=/path/to/network
COFKIT_LMP_PATH=/path/to/lmp
COFKIT_EQEQ_PATH=/path/to/eqeq
COFKIT_GRASPA_PATH=/path/to/nvc_main.x
COFKIT_RASPA2_PATH=/path/to/simulate
```

See [external-tools.md](external-tools.md) for install notes.
