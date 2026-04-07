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

## Versioning

`cofkit` uses calendar versioning in `YYYY.M.D` form for normal releases, for example `2026.4.7`.

If a second release is needed on the same day, use a PEP 440 post-release such as `2026.4.7.post1`.

## Installation and dependencies

Canonical repository install:

```bash
uv sync --locked
uv run cofkit --help
```

`uv` is the canonical install method for this repository. The base environment installs the mandatory Python runtime dependencies `gemmi`, `rdkit`, `openbabel-wheel`, `pandas`, and `pymatgen`.

For local development and verification, add the `dev` extra:

```bash
uv sync --locked --extra dev
uv run pytest -q
```

`uv run ...` is the canonical local execution path for this repo because it uses the repo-managed `.venv` and keeps helper executables aligned with the same interpreter and dependencies.

If you explicitly want an editable install inside an existing Python environment, `python3 -m pip install -e .` still works, but it is no longer the canonical repo setup.

Optional external tools you may want in your environment:

- `Zeo++` for the initial `cofkit analyze zeopp` pore-property wrapper, with the binary path provided through `COFKIT_ZEOPP_PATH`
- `LAMMPS` for the initial `cofkit calculate lammps-optimize` local optimization wrapper, with the executable path provided through `COFKIT_LMP_PATH`
- `EQeq` for the default `cofkit calculate lammps-optimize` charge-assignment stage plus the `cofkit calculate graspa-widom` and `cofkit calculate graspa-isotherm` workflows, with the executable path provided through `COFKIT_EQEQ_PATH`
- `gRASPA` for the initial `cofkit calculate graspa-widom` Widom-insertion stage and the `cofkit calculate graspa-isotherm` adsorption stage, with the executable path provided through `COFKIT_GRASPA_PATH`
- `pytest` to run the local test suite when you install the `dev` extra

The bundled topology repository under [`src/cofkit/data/topologies`](src/cofkit/data/topologies) is sufficient for normal use. External RCSR archives and topology environment variables are optional advanced inputs, not required setup steps.

### Supported External Tool Installs

For reproducible wrapper behavior, prefer these exact upstreams and binary names.

#### Zeo++

Use Zeo++ `0.3` from the upstream tarball:

- download: `http://www.zeoplusplus.org/zeo++-0.3.tar.gz`
- expected binary: `network`

Typical build:

```bash
curl -LO http://www.zeoplusplus.org/zeo++-0.3.tar.gz
gunzip zeo++-0.3.tar.gz
tar xvf zeo++-0.3.tar
cd zeo++-0.3/voro++/src
make
cd ../..
make
export COFKIT_ZEOPP_PATH="$PWD/network"
```

#### LAMMPS

Install LAMMPS from `conda-forge`:

```bash
conda create -n cofkit-lammps -c conda-forge lammps
conda activate cofkit-lammps
export COFKIT_LMP_PATH="$(command -v lmp || command -v lmp_mpi)"
```

Conda is still useful here because LAMMPS is an external executable, not part of the canonical `uv`-managed Python environment for `cofkit`. `cofkit` accepts either `lmp` or `lmp_mpi` as long as `COFKIT_LMP_PATH` points to a working executable.

#### EQeq

Use the CIF-capable fork:

- source: `https://github.com/csllpr/EQeq`
- expected binary: `eqeq`

Typical build:

```bash
git clone https://github.com/csllpr/EQeq.git
cd EQeq
g++ main.cpp -O3 -o eqeq
export COFKIT_EQEQ_PATH="$PWD/eqeq"
```

That fork reads `data/ionizationdata.dat` and `data/chargecenters.dat` relative to the executable by default, so extra path arguments are usually unnecessary. `EQEQ_IONIZATION_DATA_PATH` and `EQEQ_CHARGE_CENTERS_PATH` remain available if you need to override them.

#### gRASPA

Use one of these repositories:

- preferred: `https://github.com/csllpr/gRASPA`
- fallback upstream: `https://github.com/snurr-group/gRASPA`
- expected binary: `nvc_main.x`

The preferred fork is the one we use for higher-throughput GPU scheduling because it carries local performance tweaks. gRASPA is source-first, so the exact build depends on your NVIDIA HPC SDK / CUDA installation. The checked-in `NVC_COMPILE` script in both repositories is the reference starting point and produces `nvc_main.x`.

Typical flow:

```bash
git clone https://github.com/csllpr/gRASPA.git
cd gRASPA
# Edit NVC_COMPILE if your NVIDIA HPC SDK / CUDA paths differ.
bash NVC_COMPILE
export COFKIT_GRASPA_PATH="$PWD/nvc_main.x"
```

If your local build places the binary somewhere else, such as `src_clean/nvc_main.x`, point `COFKIT_GRASPA_PATH` there instead.

## CLI

After editable install, the main user-facing interface is the `cofkit` command.

Inspect the available commands with:

```bash
cofkit --help
cofkit build --help
cofkit analyze --help
cofkit calculate --help
cofkit build list-templates
```

The most useful grouped commands are:

- `cofkit build single-pair`
- `cofkit build batch-binary-bridge`
- `cofkit build batch-all-binary-bridges`
- `cofkit analyze classify-output`
- `cofkit analyze zeopp`
- `cofkit calculate lammps-optimize`
- `cofkit calculate graspa-widom`
- `cofkit calculate graspa-isotherm`
- `cofkit build default-library`

Legacy flat aliases such as `cofkit single-pair` and `cofkit classify-output` are still accepted for compatibility and emit deprecation warnings.

### Single pair

Direct single-pair generation with motif autodetection and default topology enumeration:

```bash
cofkit build single-pair \
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
cofkit build single-pair \
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
cofkit build batch-binary-bridge \
  --template-id imine_bridge \
  --input-dir examples/default_monomers_library \
  --output-dir out/batch_imine_from_cli \
  --max-workers 8
```

### All available binary-bridge linkages

Run all currently available binary-bridge templates over one library directory:

```bash
cofkit build batch-all-binary-bridges \
  --input-dir examples/default_monomers_library \
  --output-dir out/available_binary_bridge_batches \
  --max-workers 8
```

That CLI path reproduces the current per-template summary results for the supported `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge` batches in the shipped example library.

### Classify a finished output tree

```bash
cofkit analyze classify-output \
  out/full_cif_generation_default_selector_20260320 \
  --output-dir out/full_cif_generation_default_selector_20260320_coarse_validation_triage
```

### Run initial Zeo++ pore analysis on one CIF

```bash
export COFKIT_ZEOPP_PATH=/path/to/zeo++/network

cofkit analyze zeopp \
  out/cli_single_pair_hcb/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_zeopp \
  --json
```

By default, the Zeo++ wrapper writes a point-probe baseline:

- basic pore metrics from `-res` and `-resex`
- point-probe channel summary from `-chan 0`
- point-probe surface area from `-sa 0 0 ...`
- point-probe pore volume from `-vol 0 0 ...`

If you also want accessibility-aware probe scans, repeat `--probe-radius`:

```bash
cofkit analyze zeopp \
  out/cli_single_pair_hcb/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_zeopp \
  --probe-radius 1.20 \
  --probe-radius 1.86 \
  --json
```

Each requested probe radius adds one scan with channel, surface-area, volume, and accessibility outputs for that probe. The wrapper keeps the raw Zeo++ outputs plus stdout/stderr logs in the output directory and records a `zeopp_report.json` summary there.

### Run an initial LAMMPS local cleanup on one CIF

```bash
export COFKIT_LMP_PATH=/path/to/lmp_mpi

cofkit calculate lammps-optimize \
  out/cli_single_pair_hcb/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_lammps_opt \
  --forcefield uff \
  --json
```

The `lammps-optimize` CLI already exposes the main runtime and minimization controls:

- `--forcefield {uff}` selects the force-field backend
- `--pair-cutoff` sets the global LJ cutoff
- `--position-restraint-force-constant` controls the stage-1 local `spring/self` restraint
- `--pre-minimization-steps` plus the `--pre-minimization-*` flags tune the default restrained `langevin + nve/limit` prerun; the default is `10000` steps and `--pre-minimization-steps 0` disables it
- `--two-stage` and `--no-two-stage` plus the `--stage2-*` flags control the second minimization stage; it is enabled by default and unrestrained unless you set a stage-2 restraint explicitly
- `--energy-tolerance`, `--force-tolerance`, `--max-iterations`, `--max-evaluations`, and `--min-style` control stage 1
- `--timestep` and the `--min-modify-*` flags expose the main LAMMPS minimizer tuning knobs
- `--relax-cell` and `--no-relax-cell` control the final `fix box/relax` stage; cell relaxation is enabled by default and the `--box-relax-*` flags tune it
- `--timeout-seconds` caps wall-clock subprocess time
- `--lmp-path` overrides `COFKIT_LMP_PATH` for one run

If `OMP_NUM_THREADS` is not already set in the environment, `cofkit` now launches LAMMPS with a default of half the machine core count.

For example, a longer staged run can be launched as:

```bash
cofkit calculate lammps-optimize \
  out/cli_single_pair_hcb/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_lammps_opt \
  --pair-cutoff 14.0 \
  --position-restraint-force-constant 0.10 \
  --pre-minimization-steps 100 \
  --pre-minimization-temperature 300 \
  --pre-minimization-damping 100 \
  --pre-minimization-displacement-limit 0.05 \
  --stage2-position-restraint-force-constant 0.02 \
  --energy-tolerance 1e-8 \
  --force-tolerance 1e-8 \
  --max-iterations 5000 \
  --max-evaluations 50000 \
  --min-style fire \
  --timestep 0.5 \
  --min-modify-dmax 0.15 \
  --min-modify-fire-integrator verlet \
  --min-modify-fire-tmax 4.0 \
  --box-relax-min-style cg \
  --box-relax-vmax 0.001
```

This first LAMMPS wrapper is still conservative, but it no longer invents generic bonded/nonbonded coefficients:

- it requires a `P1` CIF with an explicit `_geom_bond_*` loop
- explicit CIF bond types in `_ccdc_geom_bond_type` are required, and `cofkit` atomistic exports now write them by default
- it can run fixed-cell minimization only, or append an optional final `fix box/relax` stage
- it builds a bonded LAMMPS data file from the explicit CIF bond graph
- the public force-field backend is `UFF`, using Open Babel UFF atom typing plus formulas and parameter tables aligned with the bundled Open Babel `UFF.prm` and the reference `lammps_interface` logic
- `UFF` is the default and currently only implemented force-field backend in the public workflow
- the current UFF path writes bond, angle, dihedral, improper, and van der Waals terms, plus optional local position restraints whose energy is explicitly included in minimization through `fix_modify energy yes`
- the default LAMMPS path stages EQeq before export, writes charged `atom_style full` data, and enables Coulomb terms in the generated LAMMPS input
- `--charge-model none` is still available when you explicitly want an uncharged export
- it writes an updated CIF, the generated LAMMPS data/input files, logs, a trajectory dump, and `lammps_report.json`

This is still a topology-preserving pre-optimization step, not a full production force-field workflow. The current UFF path is now explicit-bond-order-driven and includes torsions and impropers. Charges are now assigned through EQeq by default, but the workflow should still be treated as a serious cleanup / pre-optimization protocol rather than a final force-field-quality optimization.

### Run EQeq + gRASPA Widom insertion on one CIF

```bash
export COFKIT_EQEQ_PATH=/path/to/eqeq
export COFKIT_GRASPA_PATH=/path/to/nvc_main.x

cofkit calculate graspa-widom \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_graspa \
  --component CO2 \
  --component N2 \
  --widom-moves-per-component 300000 \
  --json
```

The `graspa-widom` wrapper runs one staged workflow:

- copy the input CIF into an `eqeq/` run directory and assign framework charges with EQeq directly from the CIF
- copy the charged framework into `widom/framework.cif`
- materialize the packaged gRASPA Widom template files in `widom/`
- compute `UnitCells` from the charged CIF cell lengths and the larger of `CutOffVDW` / `CutOffCoulomb`
- run gRASPA and parse the Widom summary values from `widom/Output/*.data`

Packaged Widom probe definitions are available for `TIP4P`, `CO2`, `H2`, `N2`, `SO2`, `Xe`, and `Kr`. Activate only the probes you want with repeated `--component NAME` flags or `--all-components`. `--widom-moves-per-component` sets the target sampling per active component, and `cofkit` derives `NumberOfProductionCycles` from that selection. The bundled wrapper now defaults `NumberOfBlocks` to `5`.

The wrapper writes:

- `eqeq/` logs plus the charged CIF emitted by EQeq
- `widom/simulation.input`, `widom/framework.cif`, and the packaged `.def` files
- `widom/Output/results.csv` with the parsed Widom energy and Henry coefficient summaries
- `graspa_widom_report.json` with paths, settings, parsed component results, and warnings

`COFKIT_EQEQ_PATH` and `COFKIT_GRASPA_PATH` can both be overridden per run with `--eqeq-path` and `--graspa-path`. On many installations the gRASPA executable is named `nvc_main.x`. If gRASPA emits non-finite uncertainty fields, `cofkit` keeps the raw data file and records those specific values as `null` in the JSON report instead of failing the whole run.

### Run EQeq + gRASPA single-component adsorption isotherms on one CIF

```bash
export COFKIT_EQEQ_PATH=/path/to/eqeq
export COFKIT_GRASPA_PATH=/path/to/nvc_main.x

cofkit calculate graspa-isotherm \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_isotherm \
  --component CO2 \
  --pressure 10000 \
  --pressure 100000 \
  --pressure 1000000 \
  --json
```

The `graspa-isotherm` wrapper is the first real-adsorption gRASPA path in `cofkit`. It runs one staged workflow:

- copy the input CIF into an `eqeq/` run directory and assign framework charges with EQeq directly from the CIF
- copy the charged framework into `isotherm/framework.cif`
- materialize the packaged gRASPA component definitions into one per-pressure run directory under `isotherm/`
- compute `UnitCells` from the charged CIF cell lengths and the larger of `CutOffVDW` / `CutOffCoulomb`
- run one single-component GCMC adsorption simulation per requested pressure point
- parse absolute loading in `mol/kg` and `g/L`, plus heat of adsorption, from each pressure-point `Output/*.data`

The current public scope is intentionally narrow: one packaged component at a time, one or more user-supplied pressure points, and one parser focused on absolute loading plus heat-of-adsorption block averages. Packaged components currently match the Widom wrapper: `TIP4P`, `CO2`, `H2`, `N2`, `SO2`, `Xe`, and `Kr`.

The wrapper writes:

- `eqeq/` logs plus the charged CIF emitted by EQeq
- `isotherm/framework.cif`
- one staged gRASPA run directory per pressure under `isotherm/pressure_*/`
- `isotherm/results.csv` with the parsed pressure/loading summary
- `graspa_isotherm_report.json` with paths, settings, parsed point results, and warnings

`--pressure` repeats to define the isotherm grid. `--production-cycles` applies per pressure point. `--fugacity-coefficient` accepts either a positive float or `PR-EOS`. If gRASPA emits non-finite loading or heat values, `cofkit` keeps the raw data file and records those specific fields as `null` in the JSON report instead of failing the whole run.

This wrapper produces pure-component adsorption points. Ratios such as `Xe/Kr` computed from separate `graspa-isotherm` runs are loading ratios, not mixture selectivities.

### Rebuild the detector-scanned example library

```bash
cofkit build default-library
```

## Skill integration

If you are using Codex inside this repository, load [skills/cofkit-navigator/SKILL.md](skills/cofkit-navigator/SKILL.md). The skill routes natural-language requests into the narrowest supported `cofkit` workflow and is intended for:

- single-pair generation from two SMILES strings
- library-scale binary-bridge screening
- classification of finished output trees
- EQeq + gRASPA Widom insertion on an exported CIF
- EQeq + gRASPA single-component adsorption isotherms on an exported CIF
- rebuilding the detector-scanned default monomer library
- choosing between the CLI, `BatchStructureGenerator`, and `COFEngine`

Typical prompts:

- "Build a single-pair imine COF from these two SMILES."
- "Run all supported binary-bridge batches over this library directory."
- "Classify this finished output tree into validation buckets."
- "Run the gRASPA Widom wrapper on this optimized CIF."
- "Run the gRASPA isotherm wrapper on this optimized CIF."
- "Should I use the CLI, `BatchStructureGenerator`, or `COFEngine` for this task?"

The skill prefers the installed `cofkit` CLI. If the package is not installed in editable mode, it falls back to an equivalent `PYTHONPATH=src ...` launcher.

## Additional documentation

- [USER_MANUAL.md](USER_MANUAL.md) for end-user Python and workflow examples
- [docs/README.md](docs/README.md) for scope, wrapper scripts, topology notes, pipeline details, and validation thresholds
- [agent-docs/README.md](agent-docs/README.md) for development and agent-facing navigation guides
- [ARCHITECTURE.md](ARCHITECTURE.md) and [CHANGELOG.md](CHANGELOG.md) for broader project context
