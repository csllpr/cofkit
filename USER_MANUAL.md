# User Manual

This manual covers the current practical workflows for generating COF structures in `cofkit`.

The examples below focus on the repository's current supported generation path:

- binary-bridge COFs generated from one role-specific monomer on each side
- current practical chemistries: imine, beta-ketoenamine, hydrazone, boronate ester, and vinylene
- supported one-node topology families in `2D` and `3D`
- CIF export for inspection
- optional post-generation triage into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid`

`cofkit` exposes a generic binary-bridge interface for registered linkage profiles. Imine examples still make the simplest introduction, but the same batch and single-pair path also supports the currently implemented non-imine binary bridges listed above.

`stacking_mode` must remain `"disabled"`.

## Choose the right interface

Use `COFEngine` when:

- you already know which topology you want
- you want a normal single-pair `CandidateEnsemble`
- you want direct `3D` single-pair generation for supported defaults such as `dia` and `pcu`

Use `BatchStructureGenerator.generate_monomer_pair_candidate(...)` or `generate_monomer_pair_candidates(...)` when:

- you already have `MonomerSpec` objects
- you want the same topology-family-aware generation used by the batch pipeline
- you want one best structure or all supported structures for one monomer pair

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

Core install:

```bash
python3 -m pip install -e .
```

Optional tools:

- `gemmi` for CIF-backed coarse validation and the broader topology scan / symmetry-expansion utilities
- `RDKit` for SMILES-based monomer construction and the practical batch workflows
- `pytest` if you want to run the tests locally

The bundled topology data shipped with the package is enough for normal structure generation. External RCSR archives are optional advanced inputs, not required setup.

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

If you use a different registered binary-bridge template, make sure the input directory contains the expected role-specific libraries for that template. Atomistic reaction realization is currently implemented for `imine_bridge`, `hydrazone_bridge`, `keto_enamine_bridge`, `boronate_ester_bridge`, and `vinylene_bridge`.

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

```bash
cofkit batch-binary-bridge \
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
cofkit batch-all-binary-bridges \
  --input-dir examples/default_monomers_library \
  --output-dir out/available_binary_bridge_batches
```

With the shipped default library snapshot, that command currently discovers `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge`.

For one direct pair:

```bash
cofkit single-pair \
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

The example scripts in `examples/` are still available, but they now act as thin wrappers over the shared `cofkit` CLI.

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
cofkit classify-output \
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
cofkit build-default-library
```

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
- `stacking_mode` must stay `"disabled"`.
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
