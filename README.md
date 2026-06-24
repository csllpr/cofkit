# cofkit

`cofkit` is a reaction-aware periodic assembly toolkit for covalent organic frameworks (COFs). It combines monomer detection, binary-bridge reaction templates, topology-aware assembly, CIF export, coarse validation, and optional external-tool wrappers.

## Quick Start

Install and verify the repository environment with `uv`:

```bash
uv sync --locked
uv run cofkit --help
```

For development and the local test suite:

```bash
uv sync --locked --extra dev
uv run pytest -q
```

The base package installs the required Python runtime dependencies: `gemmi`, `rdkit`, `openbabel`, `pandas`, and `pymatgen`.

## First Build

Generate one imine COF from two SMILES strings:

```bash
uv run cofkit build single-pair \
  --template-id imine_bridge \
  --first-smiles 'C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N' \
  --second-smiles 'C1=C(C=C(C=C1C=O)C=O)C=O' \
  --first-id tapb \
  --second-id tfb \
  --output-dir out/cli_single_pair
```

This autodetects the monomer roles, evaluates the default topology pool, writes `summary.json`, and exports CIFs under `cifs/valid`, `cifs/warning`, or `cifs/invalid`.

## Main CLI Groups

```bash
uv run cofkit build --help
uv run cofkit analyze --help
uv run cofkit calculate --help
uv run cofkit validate --help
uv run cofkit build list-templates
```

Common commands:

- `cofkit build single-pair`
- `cofkit build batch-binary-bridge`
- `cofkit build batch-all-binary-bridges`
- `cofkit build default-library`
- `cofkit analyze classify-output`
- `cofkit analyze decompose`
- `cofkit analyze zeopp`
- `cofkit calculate lammps-optimize`
- `cofkit calculate graspa-widom`
- `cofkit calculate graspa-isotherm`
- `cofkit calculate graspa-mixture`
- `cofkit calculate hybrid-mdmc`
- `cofkit validate simple`
- `cofkit validate optimize`

Legacy flat aliases such as `cofkit single-pair` and `cofkit classify-output` are still accepted for compatibility and emit deprecation warnings.

## Documentation

Start with the task page that matches what you are doing:

- [Getting started](docs/getting-started.md): install, verification, first commands
- [Building COFs](docs/building.md): single-pair, batch, topology, stacking, default libraries
- [Analysis and validation](docs/analysis-validation.md): classify, decompose, validate, Zeo++
- [Calculations](docs/calculations.md): LAMMPS, EQeq, gRASPA/RASPA2, hybrid MD/MC
- [Python API](docs/python-api.md): `COFEngine` and `BatchStructureGenerator`
- [External tools](docs/external-tools.md): Zeo++, LAMMPS, EQeq, gRASPA, RASPA2 setup
- [Current scope](docs/CURRENT_SCOPE.md): implemented capabilities and known limits
- [Documentation index](docs/README.md): complete docs map

## Optional External Tools

The Python package works without external simulation binaries for build and basic validation workflows. These optional tools enable wrapper commands:

- `Zeo++` for `cofkit analyze zeopp`, configured with `COFKIT_ZEOPP_PATH`
- `LAMMPS` for `cofkit calculate lammps-optimize`, `cofkit validate optimize`, and the MD segment of `cofkit calculate hybrid-mdmc`, configured with `COFKIT_LMP_PATH`
- `EQeq` for default charge staging in LAMMPS and gRASPA/RASPA2 workflows, configured with `COFKIT_EQEQ_PATH`
- `gRASPA` for Monte Carlo workflows, configured with `COFKIT_GRASPA_PATH`
- `RASPA2` as the CPU Monte Carlo backend, configured with `COFKIT_RASPA2_PATH`

The CLI automatically loads the nearest `.env` file by searching upward from the current working directory. Variables already present in the shell take precedence.

## Versioning

`cofkit` uses calendar versioning in `YYYY.M.D` form for normal releases, for example `2026.4.7`. If a second release is needed on the same day, use a PEP 440 post-release such as `2026.4.7.post1`.
