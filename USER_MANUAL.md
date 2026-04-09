# User Manual

This manual covers the current practical workflows for generating COF structures in `cofkit`.

The examples below focus on the repository's current supported generation path:

- binary-bridge COFs generated from one role-specific monomer on each side
- current practical chemistries: imine, hydrazone, azine, beta-ketoenamine, boronate ester, and vinylene
- supported one-node topology families in `2D` and `3D`
- CIF export for inspection
- optional post-generation triage into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid`

`cofkit` exposes a generic binary-bridge interface for registered linkage profiles. Imine examples still make the simplest introduction, but the same batch and single-pair path also supports the currently implemented non-imine binary bridges listed above.

`COFEngine` / `COFProject.stacking_mode` must remain `"disabled"`. Current stacking support is opt-in post-build enumeration for eligible `2D` batch and `single-pair` outputs via `BatchGenerationConfig(stacking_ids=...)` or CLI `--stacking ...`.

## Choose the right interface

Use `COFEngine` when:

- you already know which topology you want
- you want a normal single-pair `CandidateEnsemble`
- you want direct `3D` single-pair generation for supported defaults such as `dia` and `pcu`

Use `BatchStructureGenerator.generate_monomer_pair_candidate(...)` or `generate_monomer_pair_candidates(...)` when:

- you already have `MonomerSpec` objects
- you want the same topology-family-aware generation used by the batch pipeline
- you want one best structure or all supported structures for one monomer pair
- you want opt-in `2D` stacking enumeration on exported candidates

Use `BatchStructureGenerator.run_imine_batch(...)` or [examples/run_batch_imine_generation.py](examples/run_batch_imine_generation.py) when:

- you have text libraries of monomers
- you want to generate all compatible monomer pairs
- you want `manifest.jsonl`, `summary.md`, and per-structure CIFs

Use `BatchStructureGenerator.run_binary_bridge_batch(...)` or [examples/run_binary_bridge_generation.py](examples/run_binary_bridge_generation.py) when:

- you want the same batch workflow through the new generic registry-based API
- you want to choose the linkage template explicitly with `template_id`
- you want the same imine workflow behavior through the generic template-driven interface
- you want to infer monomer roles from generic SMILES libraries with `--auto-detect-libraries`

Use [examples/classify_batch_output.py](examples/classify_batch_output.py) when:

- you already have a finished batch output directory with `manifest.jsonl` and `cifs/`
- you want a coarse quantitative screen before manual inspection
- you want a four-way split into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid`

## Prerequisites

Use any Python `3.10+` environment.

Canonical repository install:

```bash
uv sync --locked
uv run cofkit --help
```

Installed mandatory runtime dependencies:

- `gemmi` for CIF-backed coarse validation and the broader topology scan / symmetry-expansion utilities
- `RDKit` for SMILES-based monomer construction and the practical batch workflows
- `openbabel-wheel`, `pandas`, and `pymatgen` for the current UFF/DREIDING-backed LAMMPS workflow and related report/data handling
- these are installed automatically by `uv sync --locked`

For local development and verification in this repository, use the lockfile-backed `uv` environment with the `dev` extra:

```bash
uv sync --locked --extra dev
uv run pytest -q
```

That is the canonical local test path for this repo and avoids PATH-sensitive behavior from helper executables spawned during wrapper tests. If you explicitly need an editable install inside an existing Python environment, `python3 -m pip install -e .` still works, but it is no longer the canonical repo setup.

Optional external add-ons:

- `Zeo++` if you want to use `cofkit analyze zeopp`
- `LAMMPS` if you want to use `cofkit calculate lammps-optimize`
- `EQeq` if you want to use `cofkit calculate graspa-widom`, `cofkit calculate graspa-isotherm`, or `cofkit calculate graspa-mixture` for framework charge assignment
- `gRASPA` if you want to use `cofkit calculate graspa-widom`, `cofkit calculate graspa-isotherm`, or `cofkit calculate graspa-mixture`
- `pytest` if you want to run the tests locally via the `dev` extra

The bundled topology data shipped with the package is enough for normal structure generation. External RCSR archives are optional advanced inputs, not required setup.

### Recommended external installs

For the current wrappers, use these upstreams instead of ad hoc local builds:

- `Zeo++`: use `http://www.zeoplusplus.org/zeo++-0.3.tar.gz`, build the `network` binary, then point `COFKIT_ZEOPP_PATH` at that binary
- `LAMMPS`: install from `conda-forge`, then point `COFKIT_LMP_PATH` at `lmp` or `lmp_mpi`
- `EQeq`: use `https://github.com/csllpr/EQeq`, build the `eqeq` executable, then point `COFKIT_EQEQ_PATH` at it
- `gRASPA`: prefer `https://github.com/csllpr/gRASPA`, with `https://github.com/snurr-group/gRASPA` as the fallback upstream; build `nvc_main.x` and point `COFKIT_GRASPA_PATH` at that binary

Concrete shell commands for those installs are documented in [README.md](README.md) under `Supported External Tool Installs`.

When you use the `cofkit` CLI, it automatically loads the nearest `.env` file by searching upward from the current working directory. Explicit environment variables already present in the shell still take precedence.

Typical invocation style:

```bash
python3 your_script.py
```

## Workflow 1: Generate One Pair With `COFEngine`

This is the most direct route when you already know the target topology.

### 1. Build monomers

You can create `MonomerSpec` objects manually, or build them from SMILES with RDKit.

Example from SMILES:

```python
from cofkit import COFEngine, COFProject, build_rdkit_monomer, write_candidate_cif

tapb = build_rdkit_monomer(
    "tapb",
    "TAPB",
    "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N",
    "amine",
)
tfb = build_rdkit_monomer(
    "tfb",
    "TFB",
    "C1=C(C=C(C=C1C=O)C=O)C=O",
    "aldehyde",
)
```

### 2. Choose dimensionality and topology

Current practical rules:

- `hcb`, `hca`, `fes`, `fxt`, `sql`, `kgm`, `htb`, `hxl` are `2D`
- `dia` and `pcu` are `3D`
- `3+3` node-node pairs can use `hcb`, `fes`, `fxt`
- `3+2` node-linker pairs can use `hcb`, `hca`, `fes`, `fxt`
- `4+4` and `4+2` default to `dia` in the high-connectivity `3D` route
- `6+2` defaults to `pcu` in the high-connectivity `3D` route
- chemistry-compatible indexed topologies from the bundled repository can also be requested explicitly, and the default selector now includes a curated indexed subset when the current chemistry metadata and builder support agree

### 3. Run the engine

Example: explicit `2D` `fes`

```python
project = COFProject(
    monomers=(tapb, tfb),
    allowed_reactions=("imine_bridge",),
    target_dimensionality="2D",
    target_topologies=("fes",),
)

ensemble = COFEngine().run(project)
best = ensemble.top(1)[0]

print(best.metadata["net_plan"]["topology"])
print(best.score)
write_candidate_cif("out/tapb_tfb_fes.cif", best, project.monomers)
```

Example: explicit `3D` `dia`

```python
tetra_amine = build_rdkit_monomer(
    "tetra_amine",
    "tetra_amine",
    "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1",
    "amine",
)
tetra_aldehyde = build_rdkit_monomer(
    "tetra_aldehyde",
    "tetra_aldehyde",
    "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1",
    "aldehyde",
)

project = COFProject(
    monomers=(tetra_amine, tetra_aldehyde),
    allowed_reactions=("imine_bridge",),
    target_dimensionality="3D",
    target_topologies=("dia",),
)

best = COFEngine().run(project).top(1)[0]
write_candidate_cif("out/tetra_dia.cif", best, project.monomers)
```

### 4. Request multiple topologies explicitly

`COFEngine` can return multiple candidates if you request multiple supported topologies:

```python
project = COFProject(
    monomers=(tapb, tfb),
    allowed_reactions=("imine_bridge",),
    target_dimensionality="2D",
    target_topologies=("hcb", "fes", "fxt"),
)

ensemble = COFEngine().run(project)
for candidate in ensemble.candidates:
    print(candidate.metadata["net_plan"]["topology"], candidate.score)
```

For full topology-family screening of one pair, the batch-style pair generator is usually more convenient.

## Workflow 2: Generate One Pair Across All Supported Topologies

This uses the same topology-aware pair generation logic as the batch pipeline, but on a single monomer pair.

### 1. Build or prepare `MonomerSpec` objects

```python
from cofkit import BatchGenerationConfig, BatchStructureGenerator, build_rdkit_monomer

tapb = build_rdkit_monomer(
    "tapb",
    "TAPB",
    "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N",
    "amine",
)
tpal = build_rdkit_monomer("tpal", "TPAL", "O=Cc1ccc(C=O)cc1", "aldehyde")

generator = BatchStructureGenerator(
    BatchGenerationConfig(
        write_cif=True,
        enumerate_all_topologies=True,
    )
)
```

### 2. Get the single best structure for the pair

```python
summary, candidate = generator.generate_monomer_pair_candidate(
    tapb,
    tpal,
    out_dir="out/single_pair_best",
)

print(summary.structure_id)
print(summary.topology_id)
print(summary.cif_path)
print(candidate.score if candidate is not None else None)
```

This writes one CIF for the best topology if `write_cif` is enabled.

### 3. Get all supported structures for the pair

```python
summaries, candidates, attempted = generator.generate_monomer_pair_candidates(
    tapb,
    tpal,
    out_dir="out/single_pair_all",
)

print("attempted:", attempted)
for summary in summaries:
    print(summary.topology_id, summary.score, summary.cif_path)
```

This is the easiest way to inspect all generated topologies for one pair.

### 4. Restrict topology selection

You can limit generation to specific topology ids:

```python
generator = BatchStructureGenerator(
    BatchGenerationConfig(
        single_node_topology_ids=("hcb", "fes"),
        write_cif=True,
    )
)
```

That restriction applies to:

- `generate_pair_candidate(...)`
- `generate_pair_candidates(...)`
- `generate_monomer_pair_candidate(...)`
- `generate_monomer_pair_candidates(...)`
- `run_imine_batch(...)`

## Workflow 3: Generate From Text Libraries

The explicit batch pipeline auto-discovers library files named like:

- `amines_count_2.txt`
- `amines_count_3.txt`
- `amines_count_4.txt`
- `amines_count_6.txt`
- `aldehydes_count_2.txt`
- `aldehydes_count_3.txt`
- `aldehydes_count_4.txt`
- `aldehydes_count_6.txt`
- `hydrazides_count_2.txt`
- `boronic_acids_count_2.txt`
- `boronic_acids_count_3.txt`
- `catechols_count_2.txt`
- `catechols_count_3.txt`
- `keto_aldehydes_count_1.txt`
- `keto_aldehydes_count_2.txt`
- `keto_aldehydes_count_3.txt`
- `keto_aldehydes_count_4.txt`
- `activated_methylenes_count_N.txt` when those libraries are present

Each file should contain one SMILES string per line. A leading `smiles` header is allowed.

Example file:

```text
smiles
Nc1ccc(N)cc1
Nc1ccnc(N)c1
```

### 1. Legacy imine convenience entry point

```bash
python3 examples/run_batch_imine_generation.py \
  --input-dir examples/batch_test_monomers \
  --output-dir out/batch_imine_generation \
  --max-cif-exports 20000 \
  --num-conformers 2
```

### 2. Generic binary-bridge entry point

This is the same practical workflow through the generic template-driven API:

```bash
python3 examples/run_binary_bridge_generation.py \
  --template-id imine_bridge \
  --input-dir examples/batch_test_monomers \
  --output-dir out/binary_bridge_generation \
  --max-cif-exports 20000 \
  --num-conformers 2
```

With `--template-id imine_bridge`, this produces the same structures as the imine convenience entry point while exercising the registry-driven binary-bridge path.

The generic batch path now uses the default `8`-worker process pool for pair generation unless you override `--max-workers` or the Python config.

If you use a different registered binary-bridge template, make sure the input directory contains the expected role-specific libraries for that template. Atomistic reaction realization is currently implemented for `imine_bridge`, `hydrazone_bridge`, `azine_bridge`, `keto_enamine_bridge`, `boronate_ester_bridge`, and `vinylene_bridge`.

### Auto-detected libraries

You can also point the generic runner at raw plain `.txt` SMILES libraries and let the detector regroup them by monomer role/connectivity:

```bash
python3 examples/run_binary_bridge_generation.py \
  --template-id keto_enamine_bridge \
  --input-dir /path/to/raw_generic_smiles_libraries \
  --output-dir out/binary_bridge_generation_keto \
  --auto-detect-libraries \
  --num-conformers 2
```

This is useful when the filenames do not already encode role information. The autodetect path is strict: every non-empty row must resolve to one of the selected template's supported motif kinds, or the batch run fails. The detector currently recognizes the practical binary-bridge motif kinds supported by the RDKit path:

- `amine`
- `aldehyde`
- `hydrazine`
- `hydrazide`
- `boronic_acid`
- `catechol`
- `keto_aldehyde`
- `activated_methylene`

When a monomer matches more than one currently supported role, the autodetect path records it as an ambiguity instead of silently choosing.

### Generated example library

The repository also ships a detector-scanned explicit library at [examples/default_monomers_library](examples/default_monomers_library). It contains:

- regrouped `*_count_N.txt` libraries derived from the current detector
- `registry.jsonl` with detector metadata and source provenance
- `failures.jsonl` with failed or ambiguous autodetections

Use that directory directly in batch commands without `--auto-detect-libraries`; it is already grouped output from the detector.

To regenerate it from the raw fixture source set:

```bash
python3 examples/build_default_monomers_library.py
```

### CLI usage

The installable CLI is grouped under `cofkit build`, `cofkit analyze`, and `cofkit calculate`. Legacy flat aliases such as `cofkit single-pair` and `cofkit classify-output` still work for now, but they emit deprecation warnings.

```bash
cofkit build batch-binary-bridge \
  --template-id imine_bridge \
  --input-dir examples/batch_test_monomers \
  --output-dir out/batch_imine_generation \
  --num-conformers 4
```

Useful flags:

- `--max-pairs 100` to cap the number of attempted monomer pairs
- `--no-write-cif` to skip CIF export
- `--max-cif-exports 500` to limit export count
- `--no-all-topologies` to keep only the best topology per pair
- `--auto-detect-libraries` to infer monomer roles/connectivities from SMILES instead of relying only on filename prefixes

For one batch per currently available binary-bridge linkage:

```bash
cofkit build batch-all-binary-bridges \
  --input-dir examples/default_monomers_library \
  --output-dir out/available_binary_bridge_batches
```

With the shipped default library snapshot, that command currently discovers `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge`.

For one direct pair:

```bash
cofkit build single-pair \
  --template-id imine_bridge \
  --first-smiles 'C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N' \
  --second-smiles 'C1=C(C=C(C=C1C=O)C=O)C=O' \
  --first-id tapb \
  --second-id tfb \
  --output-dir out/cli_single_pair
```

The single-pair CLI autodetects motif kinds by default. If you want to force them explicitly, pass `--first-motif-kind` and `--second-motif-kind`. If you want to force topology selection instead of using the default topology pool, pass one or more `--topology` values. The direct single-pair output writes:

- `summary.json`
- one CIF per written topology under `cifs/valid`, `cifs/warning`, or `cifs/invalid`

For eligible `2D` outputs, you can also request named bilayer registries during export:

```bash
cofkit build single-pair \
  --template-id imine_bridge \
  --first-smiles 'C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N' \
  --second-smiles 'C1=C(C=C(C=C1C=O)C=O)C=O' \
  --first-id tapb \
  --second-id tfb \
  --topology hcb \
  --stacking AA \
  --stacking AB \
  --output-dir out/cli_single_pair_stacked
```

That path keeps the normal build pipeline unchanged, then emits one exported structure per requested registry. Written CIFs append the registry tag to the `COFid` comment line, for example `# COFid: ... stacking=AA`.

The example scripts in `examples/` are still available, but they now act as thin wrappers over the shared grouped `cofkit` CLI.

### Python usage

```python
from cofkit import BatchGenerationConfig, BatchStructureGenerator

generator = BatchStructureGenerator(
    BatchGenerationConfig(
        rdkit_num_conformers=4,
        enumerate_all_topologies=True,
        write_cif=True,
    )
)

summary = generator.run_imine_batch(
    "examples/batch_test_monomers",
    "out/batch_imine_generation",
)

print(summary.attempted_pairs)
print(summary.successful_structures)
print(summary.manifest_path)
```

For the generic autodetect path:

```python
summary = generator.run_binary_bridge_batch(
    "examples/default_monomers_library",
    "out/binary_bridge_generation_keto",
    template_id="keto_enamine_bridge",
    auto_detect_libraries=True,
)
```

## Output Files

### Single-pair generation

Depending on the function you call:

- one best-topology CIF in the chosen output directory
- or one CIF per generated topology

## Workflow 5: Literature-Style Smoke Examples

For a quick sanity check across the currently implemented binary-bridge chemistries, run:

```bash
python3 examples/run_experimental_examples.py
```

This writes per-example output plus a small report covering:

- `COF-300` as an imine `dia` case
- `TpPa-1` as a beta-ketoenamine `hcb` case
- `COF-42` as a hydrazone `hcb` case
- `COF-5` as a boronate-ester `hcb` case
- `V-COF-1` as a vinylene `hcb` case

### Batch generation

The batch workflow writes:

- `manifest.jsonl`: one JSON object per generated structure
- `summary.md`: human-readable run summary
- `cifs/valid/`, `cifs/warning/`, and `cifs/invalid/`: one CIF per exported structure if CIF export is enabled

Structures classified as `hard_hard_invalid` are still recorded in the manifest, but their CIF export is blocked and `cif_export_blocked = true` is set in the per-structure metadata.

## Workflow 4: Classify A Finished Batch Output

Once you have a large CIF set, the classifier can split it into clearly acceptable, borderline, and clearly broken categories for manual review.

### CLI usage

```bash
cofkit analyze classify-output \
  out/full_cif_generation_default_selector_20260320 \
  --output-dir out/full_cif_generation_default_selector_20260320_coarse_validation_triage
```

By default this writes symlink-based views of the source CIFs instead of copying them.

### Output structure

The classifier writes:

- `classification_manifest.jsonl`: one JSON object per classified structure
- `valid/manifest.jsonl` plus `valid/cifs/`
- `warning/manifest.jsonl`, `warning/cifs/`, and `warning/reasons/<reason>/`
- `hard_hard_invalid/manifest.jsonl`, `hard_hard_invalid/cifs/`, and `hard_hard_invalid/reasons/<reason>/`
- `hard_invalid/manifest.jsonl`, `hard_invalid/cifs/`, and `hard_invalid/reasons/<reason>/`

To regenerate the detector-scanned default library from the raw fixture inputs:

```bash
cofkit build default-library
```

## Workflow 6: Initial Zeo++ Pore Analysis

The first `analyze`-namespace external-tool wrapper is Zeo++. It now writes a point-probe pore baseline for one CIF at a time and can optionally add repeated accessibility-aware probe scans.

### Environment setup

Point `cofkit` at the Zeo++ `network` binary through an environment variable:

```bash
export COFKIT_ZEOPP_PATH=/path/to/zeo++/network
```

### CLI usage

```bash
cofkit analyze zeopp \
  out/cli_single_pair_hcb/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_zeopp \
  --json
```

Default output includes:

- largest included sphere
- largest free sphere
- largest included sphere along the maximum free-sphere path
- axis-resolved free-sphere and included-sphere-along-path values from `-resex`
- point-probe channel summary from `-chan 0`
- point-probe surface area from `-sa 0 0 ...`
- point-probe pore volume from `-vol 0 0 ...`

To add accessibility-aware probe scans, repeat `--probe-radius`:

```bash
cofkit analyze zeopp \
  out/cli_single_pair_hcb/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_zeopp \
  --probe-radius 1.20 \
  --probe-radius 1.86 \
  --json
```

Each requested probe radius adds:

- probe-specific channel summary
- probe-specific accessible/non-accessible surface area
- probe-specific accessible/non-accessible pore volume
- Voronoi-node accessibility counts from `-axs`

The output directory stores:

- raw Zeo++ output files for the baseline and each probe scan
- stdout/stderr logs for every Zeo++ subprocess call
- `zeopp_report.json`

Current scope note: the public wrapper currently focuses on basic pore metrics plus surface-area / pore-volume summaries. Richer Zeo++ modes such as PSD histograms, grids, ray analyses, ZeoVis exports, and other hidden commands remain out of the public `cofkit` CLI for now.

## Workflow 7: Initial LAMMPS CIF Cleanup

The first `calculate`-namespace external-tool wrapper is a conservative LAMMPS local cleanup for explicit-bond CIFs.

### Environment setup

Point `cofkit` at the LAMMPS executable through an environment variable:

```bash
export COFKIT_LMP_PATH=/path/to/lmp_mpi
```

Alternatively, pass `--lmp-path` per run if you do not want to export `COFKIT_LMP_PATH`.

### CLI usage

```bash
cofkit calculate lammps-optimize \
  out/cli_single_pair_hcb/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_lammps_opt \
  --forcefield uff \
  --json
```

### Exposed CLI controls

The public `lammps-optimize` command already exposes the main optimization settings directly:

- `--forcefield {uff,dreiding}` selects the force-field backend
- `--pair-cutoff` sets the LJ cutoff used by the current backend
- `--position-restraint-force-constant` controls the stage-1 local `spring/self` restraint
- `--pre-minimization-steps` plus the `--pre-minimization-*` flags tune the default restrained `langevin + nve/limit` prerun; the default is `10000` steps and `--pre-minimization-steps 0` disables it
- `--two-stage` and `--no-two-stage` plus the `--stage2-*` flags control the second minimization stage; it is enabled by default and unrestrained unless you set a stage-2 restraint explicitly
- `--energy-tolerance`, `--force-tolerance`, `--max-iterations`, `--max-evaluations`, and `--min-style` control stage 1
- `--timestep` and the `--min-modify-*` flags expose the main LAMMPS minimizer tuning controls
- `--relax-cell` and `--no-relax-cell` control the final `fix box/relax` stage; cell relaxation is enabled by default and the `--box-relax-*` flags tune it
- `--timeout-seconds` sets the subprocess timeout
- `--lmp-path` overrides `COFKIT_LMP_PATH` for one run

If `OMP_NUM_THREADS` is not already set, `cofkit` launches LAMMPS with a default of half the machine core count.

A more explicitly tuned run looks like:

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

Current behavior:

- input must be a `P1` CIF with an explicit `_geom_bond_*` loop
- `cofkit` atomistic CIF exports now write `_ccdc_geom_bond_type`, and the LAMMPS wrapper requires that explicit bond-type field for every bond
- the run can stay fixed-cell or append a final `fix box/relax` stage
- the explicit CIF bond graph drives the bonded topology written into the LAMMPS data file
- `UFF` is the default backend and uses Open Babel UFF atom typing plus formulas and parameter tables aligned with the bundled pinned Open Babel `UFF.prm` reference
- `DREIDING` is also available and uses Open Babel UFF atom typing mapped onto a pinned `lammps-interface` DREIDING table plus DREIDING-style coefficient formulas
- both public backends currently write bond, angle, dihedral, improper, and van der Waals terms
- for real optimization work, prefer `DREIDING`; `UFF` remains available for compatibility and comparison runs, but should currently be treated as experimental support
- the default LAMMPS path runs EQeq before export, writes charged `atom_style full` data, and enables Coulomb terms in the generated LAMMPS input
- `--charge-model none` remains available when you explicitly want an uncharged export
- optional `spring/self` restraints keep the optimization local by default, and their energy is included in the minimization objective via `fix_modify energy yes`

Temporary parameter review note:

- the current `DREIDING` implementation is pinned to `lammps-interface` commit `255f027cb76142d39c050a6810404debc6a06562`, and that upstream is unmaintained
- a few heavier-atom DREIDING entries there are explicitly heuristic (`Cu`, `Ni`, `Mg`)
- hydrogen-bond-specific DREIDING parameters from that upstream are not part of the current `cofkit` export
- some historical rounded DREIDING/gRASPA tables differ slightly from the current generated values
- current `UFF` gRASPA rows are generated from the bundled Open Babel `UFF.prm`, so they may differ slightly from rounded example files

The output directory stores:

- the generated `lammps_input.data`
- the generated `lammps_minimize.in`
- `lammps.log`, `lammps.stdout.log`, and `lammps.stderr.log`
- `lammps_trajectory.lammpstrj`
- `lammps_report.json`
- one updated `*_lammps_optimized.cif`

Current scope note: this is a topology-preserving local cleanup step for generated explicit-bond COF CIFs. The current `UFF` and `DREIDING` exports are explicit-bond-order-driven and include torsion and improper terms. EQeq charges are now part of the default LAMMPS export, but even with staged minimization and optional box relaxation it should still be treated as a pre-optimization candidate generator rather than a final optimized structure.

## Workflow 8: EQeq + gRASPA Widom Insertion

The second `calculate`-namespace external-tool wrapper is a staged `EQeq -> gRASPA` Widom insertion workflow for one CIF.

### Environment setup

Point `cofkit` at the EQeq and gRASPA executables through environment variables:

```bash
export COFKIT_EQEQ_PATH=/path/to/eqeq
export COFKIT_GRASPA_PATH=/path/to/nvc_main.x
```

The gRASPA executable is commonly named `nvc_main.x`. The CLI also auto-loads the nearest `.env` file, so these variables can live there instead of being exported in every shell. Both executables can also be overridden per run with `--eqeq-path` and `--graspa-path`.

### CLI usage

```bash
cofkit calculate graspa-widom \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_graspa \
  --forcefield dreiding \
  --component CO2 \
  --component N2 \
  --widom-moves-per-component 300000 \
  --json
```

The wrapper runs these stages:

- copy the input CIF into `eqeq/` and run EQeq directly on that CIF
- take the charged EQeq CIF output and copy it to `widom/framework.cif`
- copy the packaged gRASPA adsorbate/pseudo-atom `.def` files into `widom/` and generate `force_field_mixing_rules.def` for the selected framework forcefield
- render `widom/simulation.input`
- derive `UnitCells` from the charged CIF cell lengths and the larger of `CutOffVDW` / `CutOffCoulomb`
- run gRASPA inside `widom/`
- parse `widom/Output/*.data` into structured component results

No Open Babel conversion is part of this workflow. The current public wrapper assumes EQeq accepts ordinary CIF inputs directly.

### Exposed CLI controls

The public `graspa-widom` command exposes the main staging and runtime controls directly:

- `--eqeq-path` and `--graspa-path` override `COFKIT_EQEQ_PATH` and `COFKIT_GRASPA_PATH`
- `--eqeq-lambda`, `--eqeq-h-i0`, `--eqeq-charge-precision`, `--eqeq-method`, `--eqeq-real-space-cells`, `--eqeq-reciprocal-space-cells`, and `--eqeq-eta` control the EQeq stage
- `--component NAME` repeats to activate packaged Widom probe molecules on demand, and `--all-components` activates the full packaged set
- `--forcefield {dreiding,uff}` selects the generated framework mixing rules
- `--widom-moves-per-component` controls the target Widom sampling per active component; `cofkit` derives `NumberOfProductionCycles` from that target unless `--production-cycles` is explicitly provided
- `--temperature`, `--pressure`, `--initialization-cycles`, `--equilibration-cycles`, `--production-cycles`, `--trial-positions`, `--trial-orientations`, `--cutoff-vdw`, `--cutoff-coulomb`, and `--ewald-precision` control the generated gRASPA `simulation.input`
- `--eqeq-timeout-seconds` and `--graspa-timeout-seconds` cap the subprocess wall-clock runtime

### Available packaged Widom probes

The packaged Widom template currently ships definitions for:

- `TIP4P`
- `CO2`
- `H2`
- `N2`
- `SO2`
- `Xe`
- `Kr`

Activate probes explicitly with repeated `--component NAME` flags or `--all-components`. The generated `simulation.input` always uses CIF-backed framework input plus charges read from that charged CIF. `NumberOfBlocks` now defaults to `5`, and when `--production-cycles` is omitted, `cofkit` sets `NumberOfProductionCycles = (--widom-moves-per-component) * (number of active components)`.

For real adsorption calculations, prefer `--forcefield dreiding`. `UFF` is available for comparison and early support, but should currently be treated as experimental.

### Output structure

The output directory stores:

- `eqeq/` with the copied input CIF, EQeq stdout/stderr logs, the charged CIF, and the optional EQeq JSON output if the executable writes it
- `widom/framework.cif`
- `widom/simulation.input`
- the packaged Widom adsorbate `.def` files plus one generated `force_field_mixing_rules.def` copied into `widom/`
- `widom/graspa.stdout.log` and `widom/graspa.stderr.log`
- one or more `widom/Output/*.data` files from gRASPA
- `widom/Output/results.csv`
- `graspa_widom_report.json`

`results.csv` records the parsed Widom energy and Henry coefficient summaries for each component. `graspa_widom_report.json` also records the resolved executable paths, the computed `UnitCells`, the effective EQeq and Widom settings, and any warnings collected during parsing.

### Current scope note

This wrapper is intentionally narrow. It currently exposes one packaged Widom-template family, one selectable packaged probe set, one selectable framework forcefield family (`DREIDING` or `UFF`), and one parser focused on Widom energy plus Henry coefficient summaries. If gRASPA emits non-finite uncertainty values, `cofkit` preserves the raw `.data` file, records those specific fields as `null` in `graspa_widom_report.json`, and leaves the rest of the parsed result intact. The temporary parameter review note from the LAMMPS section applies here as well: `DREIDING` is generated from the pinned `lammps-interface` reference and `UFF` is generated from the bundled Open Babel `UFF.prm`.

## Workflow 9: EQeq + gRASPA Single-Component Adsorption Isotherms

The third `calculate`-namespace external-tool wrapper is a staged `EQeq -> gRASPA` adsorption workflow for one CIF and one packaged adsorbate component over one or more pressure points.

### Environment setup

Point `cofkit` at the EQeq and gRASPA executables through environment variables:

```bash
export COFKIT_EQEQ_PATH=/path/to/eqeq
export COFKIT_GRASPA_PATH=/path/to/nvc_main.x
```

The gRASPA executable is commonly named `nvc_main.x`. The CLI also auto-loads the nearest `.env` file, so these variables can live there instead of being exported in every shell. Both executables can also be overridden per run with `--eqeq-path` and `--graspa-path`.

### CLI usage

```bash
cofkit calculate graspa-isotherm \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_isotherm \
  --forcefield dreiding \
  --component CO2 \
  --pressure 10000 \
  --pressure 100000 \
  --pressure 1000000 \
  --json
```

The wrapper runs these stages:

- copy the input CIF into `eqeq/` and run EQeq directly on that CIF
- take the charged EQeq CIF output and copy it to `isotherm/framework.cif`
- for each requested pressure point, copy the packaged gRASPA component `.def` files plus one generated framework `force_field_mixing_rules.def` and `framework.cif` into one `isotherm/pressure_*/` run directory
- render one pressure-specific `simulation.input` per run directory
- derive one shared `UnitCells` setting from the charged CIF cell lengths and the larger of `CutOffVDW` / `CutOffCoulomb`
- run one single-component GCMC adsorption simulation per pressure point
- parse `Output/*.data` into structured pressure/loading/heat results

No Open Babel conversion is part of this workflow. The current public wrapper assumes EQeq accepts ordinary CIF inputs directly.

### Exposed CLI controls

The public `graspa-isotherm` command exposes the main staging and runtime controls directly:

- `--eqeq-path` and `--graspa-path` override `COFKIT_EQEQ_PATH` and `COFKIT_GRASPA_PATH`
- `--eqeq-lambda`, `--eqeq-h-i0`, `--eqeq-charge-precision`, `--eqeq-method`, `--eqeq-real-space-cells`, `--eqeq-reciprocal-space-cells`, and `--eqeq-eta` control the EQeq stage
- `--component NAME` selects one packaged adsorbate definition
- `--forcefield {dreiding,uff}` selects the generated framework mixing rules
- `--pressure PA` repeats to define one or more pressure points in Pa for the isotherm grid
- `--fugacity-coefficient VALUE` accepts either a positive float or `PR-EOS`
- `--temperature`, `--initialization-cycles`, `--equilibration-cycles`, `--production-cycles`, `--trial-positions`, `--trial-orientations`, `--cutoff-vdw`, `--cutoff-coulomb`, and `--ewald-precision` control the generated gRASPA `simulation.input`
- `--eqeq-timeout-seconds` and `--graspa-timeout-seconds` cap the subprocess wall-clock runtime

### Available packaged adsorption components

The packaged adsorption template currently ships definitions for:

- `TIP4P`
- `CO2`
- `H2`
- `N2`
- `SO2`
- `Xe`
- `Kr`

Select exactly one component with `--component NAME`. The generated `simulation.input` always uses CIF-backed framework input plus charges read from that charged CIF. `NumberOfBlocks` defaults to `5`, and `--production-cycles` applies per pressure point rather than across the full pressure sweep.

For real adsorption calculations, prefer `--forcefield dreiding`. `UFF` is available for comparison and early support, but should currently be treated as experimental.

### Output structure

The output directory stores:

- `eqeq/` with the copied input CIF, EQeq stdout/stderr logs, the charged CIF, and the optional EQeq JSON output if the executable writes it
- `isotherm/framework.cif`
- one per-pressure run directory under `isotherm/pressure_*/`
- one `simulation.input` plus gRASPA stdout/stderr logs per pressure-point directory
- one or more `Output/*.data` files under each pressure-point directory
- `isotherm/results.csv`
- `graspa_isotherm_report.json`

`results.csv` records one parsed adsorption point per requested pressure, including absolute loading in `mol/kg` and `g/L` plus heat of adsorption in `kJ/mol`. `graspa_isotherm_report.json` also records the resolved executable paths, the computed `UnitCells`, the effective EQeq and adsorption settings, and any warnings collected during parsing.

### Current scope note

This wrapper is intentionally narrow. It currently exposes one packaged adsorbate at a time, one or more explicit pressure points, one selectable framework forcefield family (`DREIDING` or `UFF`), and one parser focused on absolute loading plus heat-of-adsorption block averages. It does not yet implement restart-chained pressure stepping, excess-loading post-processing, or more advanced gRASPA sampling modes. If gRASPA emits non-finite loading or heat uncertainty values, `cofkit` preserves the raw `.data` file, records those specific fields as `null` in `graspa_isotherm_report.json`, and leaves the rest of the parsed result intact. Ratios computed from separate pure-component `graspa-isotherm` runs should be treated as loading ratios only; use `graspa-mixture` for true mixed-feed selectivity.

## Workflow 10: EQeq + gRASPA Mixture Adsorption and Selectivity

The fourth `calculate`-namespace gRASPA wrapper is a staged `EQeq -> gRASPA` mixture adsorption workflow for one CIF, two or more packaged adsorbates, and one or more pressure points.

### CLI usage

```bash
cofkit calculate graspa-mixture \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_mixture \
  --forcefield dreiding \
  --component Kr:0.1 \
  --component Xe:0.9 \
  --pressure 10000 \
  --pressure 100000 \
  --fugacity-coefficient PR-EOS \
  --json
```

The wrapper runs these stages:

- copy the input CIF into `eqeq/` and run EQeq directly on that CIF
- take the charged EQeq CIF output and copy it to `mixture/framework.cif`
- for each requested pressure point, copy the packaged gRASPA adsorbate `.def` files plus one generated framework `force_field_mixing_rules.def` and `framework.cif` into one `mixture/pressure_*/` run directory
- render one pressure-specific multi-component `simulation.input` per run directory, including `MolFraction`, `IdentityChangeProbability`, and `SwapProbability` entries for every adsorbate
- derive one shared `UnitCells` setting from the charged CIF cell lengths and the larger of `CutOffVDW` / `CutOffCoulomb`
- run one multi-component GCMC adsorption simulation per pressure point
- parse per-component loading/heat summaries, compute adsorbed mole fractions, and compute pairwise selectivities `(x_i / x_j) / (y_i / y_j)`

### Exposed CLI controls

The public `graspa-mixture` command exposes the same EQeq staging controls as the Widom/isotherm wrappers plus:

- repeated `--component NAME:FRACTION` flags to define the mixed feed
- `--forcefield {dreiding,uff}` to select the generated framework mixing rules
- repeated `--pressure PA` flags to define the pressure grid
- `--fugacity-coefficient VALUE` to apply either a positive float or `PR-EOS` to every component
- `--translation-probability`, `--rotation-probability`, `--reinsertion-probability`, `--identity-change-probability`, `--swap-probability`, and `--create-number-of-molecules` to control the generated per-component GCMC move set

### Output structure

The output directory stores:

- `eqeq/` with the copied input CIF, EQeq stdout/stderr logs, the charged CIF, and the optional EQeq JSON output if the executable writes it
- `mixture/framework.cif`
- one per-pressure run directory under `mixture/pressure_*/`
- one `simulation.input` plus gRASPA stdout/stderr logs per pressure-point directory
- one or more `Output/*.data` files under each pressure-point directory
- `mixture/component_results.csv`
- `mixture/selectivity_results.csv`
- `graspa_mixture_report.json`

`component_results.csv` records one parsed adsorption row per pressure/component, including feed mole fraction, adsorbed mole fraction, loading in `mol/kg` and `g/L`, and heat of adsorption. `selectivity_results.csv` records ordered pairwise selectivities plus propagated selectivity error bars. `graspa_mixture_report.json` records the staged paths, computed `UnitCells`, effective EQeq and mixture settings, parsed point results, and warnings.

Real gRASPA runs may still print a missing-restartfile line to stderr even when `RestartFile no` is set. `cofkit` currently preserves that stderr content and surfaces it as a warning if the simulation otherwise completes and parses successfully. The same temporary parameter review note from the Widom/isotherm sections applies here.

For real adsorption calculations, prefer `--forcefield dreiding`. `UFF` is available for comparison and early support, but should currently be treated as experimental.

### Quantitative rules

Current `warning` criteria:

- any bridge distance residual `> 0.75 A`
- mean bridge distance residual `> 0.35 A`
- more than `25%` of bridge events with residual `> 0.50 A`

Current `hard_invalid` criteria:

- any unreacted motifs
- missing or unparsable CIF
- any bridge distance residual `> 1.00 A`
- mean bridge distance residual `> 0.60 A`
- any bridge `actual_distance / target_distance < 0.70`
- any bridge `actual_distance / target_distance > 1.60`
- disconnected monomer-instance graph reconstructed from CIF bonding
- any nonbonded heavy-atom contact `< 1.05 A`
- `2D` cell area `< 10.0 A^2`
- `3D` cell volume `< 20.0 A^3`

Current `hard_hard_invalid` criterion:

- any actual bridge distance `>= 2.5 A`

Hydrogen atoms are ignored in the clash check. Warning-level cases still undergo CIF-backed checks, so a structure can be promoted from `warning` to `hard_invalid` if the exported network is disconnected or physically impossible. During generation, `hard_hard_invalid` structures are not written as CIFs at all.

## How to inspect results

For a generated `Candidate`, the most useful fields are:

```python
candidate.metadata["net_plan"]["topology"]
candidate.metadata["embedding"]["placement_mode"]
candidate.metadata["score_breakdown"]
candidate.metadata["score_metadata"]["bridge_geometry_residual"]
candidate.metadata["graph_summary"]
```

For pair-level summaries from the batch generator:

```python
summary.topology_id
summary.score
summary.cif_path
summary.metadata["embedding"]
summary.metadata["score_breakdown"]
```

## Common patterns

### Screen one pair across all supported topologies

Use:

```python
generator.generate_monomer_pair_candidates(...)
```

### Build one known structure for export

Use:

```python
COFEngine().run(project)
```

with an explicit `target_topologies=(...)`.

### Screen a whole library

Use:

```python
generator.run_binary_bridge_batch(...)
```

## Common pitfalls

- The current practical route is binary-bridge generation with one resolved motif role on each side. The shipped batch-library examples currently cover imine, hydrazone, and keto-enamine workflows; boronate ester and vinylene are demonstrated through the single-pair smoke examples.
- `COFEngine` / `COFProject.stacking_mode` must still stay `"disabled"`; the current stacking support is only the post-build `BatchStructureGenerator` / CLI registry enumeration path.
- If RDKit is unavailable, SMILES-based monomer construction will fail.
- High-connectivity `4`- and `6`-connected cases are usually better handled in the `3D` route (`dia`, `pcu`) than by forcing `2D` topologies.
- `COFEngine` without explicit `target_topologies` still follows its own normal planner behavior for `2D`; if you want a full supported-topology sweep for one pair, use `BatchStructureGenerator.generate_monomer_pair_candidates(...)`.
- A `warning` classification means “inspect manually”; a `hard_invalid` classification means the current coarse screen found a broken network, impossible bridge length, or severe heavy-atom clash.
- A `hard_hard_invalid` classification means the bridge geometry is so extreme that the batch writer refuses to emit a CIF.

## Repository examples

See also:

- [examples/mock_imine.py](examples/mock_imine.py)
- [examples/export_mock_imine_cif.py](examples/export_mock_imine_cif.py)
- [examples/run_tapb_tfb_hcb_case.py](examples/run_tapb_tfb_hcb_case.py)
- [examples/run_batch_imine_generation.py](examples/run_batch_imine_generation.py)
- [examples/run_binary_bridge_generation.py](examples/run_binary_bridge_generation.py)
- [examples/run_available_binary_bridge_batches.py](examples/run_available_binary_bridge_batches.py)
- [examples/classify_batch_output.py](examples/classify_batch_output.py)
