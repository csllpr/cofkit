# Intent Routing

This reference maps natural-language requests onto the current `cofkit` interfaces.

## Command Launchers

Use one of these equivalent launch styles:

- installed CLI: `cofkit ...`
- repo-local fallback: `PYTHONPATH=src python -m cofkit.cli ...`

The repo-local fallback still requires the project dependencies to be installed in that interpreter.

## Request Patterns

### Show what the toolkit can do today

Use the grouped help plus build-template discovery:

```bash
cofkit --help
cofkit analyze --help
cofkit calculate --help
cofkit build list-templates --json
```

Interpretation rules:

- `cofkit build list-templates --json` only describes the build chemistry surface, not every public CLI workflow.
- Focus on templates where `supports_pair_generation` is `true` for practical structure generation.
- Registered templates with `supports_pair_generation=false` are still part of the reaction library, but they are not on the current topology-guided pair-generation path.
- `cofkit analyze --help` is where current pore-analysis and output-triage workflows show up.
- `cofkit calculate --help` is where current external optimization workflows show up.

### Build one COF from two monomers

Use `build single-pair`.

```bash
cofkit build single-pair \
  --template-id imine_bridge \
  --first-smiles '<FIRST_SMILES>' \
  --second-smiles '<SECOND_SMILES>' \
  --first-id monomer_a \
  --second-id monomer_b \
  --output-dir out/single_pair_run
```

When `--topology` is omitted, the current CLI keeps `--all-topologies` enabled and enumerates every applicable topology for the pair.

Useful variants:

- add `--first-motif-kind ...` and `--second-motif-kind ...` when the roles are already known
- add `--topology hcb` to force one topology
- add `--target-dimensionality 3D` for supported higher-connectivity cases
- add `--no-all-topologies` when the user wants only the best topology instead of the full supported set

Primary artifact:

- `summary.json`

### Screen one linkage over a library directory

Use `build batch-binary-bridge`.

```bash
cofkit build batch-binary-bridge \
  --template-id hydrazone_bridge \
  --input-dir <LIBRARY_DIR> \
  --output-dir out/hydrazone_batch \
  --max-workers 8
```

When `--topology` is omitted, the current CLI keeps `--all-topologies` enabled and attempts every applicable topology for each compatible monomer pair.

Add `--auto-detect-libraries` only when `<LIBRARY_DIR>` contains raw generic `.txt` SMILES files rather than grouped files like `amines_count_3.txt`.

Useful variants:

- add one or more `--topology ...` values to restrict the run to explicit topologies
- add `--no-all-topologies` when the user wants only the best topology per pair instead of the full supported set

Primary artifacts:

- `manifest.jsonl`
- `summary.md`
- `cifs/`

### Run every supported binary-bridge workflow over a library directory

Use `build batch-all-binary-bridges`.

```bash
cofkit build batch-all-binary-bridges \
  --input-dir <LIBRARY_DIR> \
  --output-dir out/all_binary_bridge_batches \
  --max-workers 8 \
  --template-workers 1
```

Primary artifacts:

- `combined_summary.json`
- `<output-dir>/<template_id>/manifest.jsonl`
- `<output-dir>/<template_id>/summary.md`

### Find which binary-bridge templates a library can support

There is no dedicated CLI command for discovery-only by library. Use the Python surface:

```python
from cofkit import BatchGenerationConfig, BatchStructureGenerator

generator = BatchStructureGenerator(BatchGenerationConfig())
template_ids = generator.available_binary_bridge_template_ids(
    "examples/default_monomers_library",
    auto_detect_libraries=False,
)
print(template_ids)
```

Switch `auto_detect_libraries=True` only for raw generic libraries.

### Rebuild the grouped example library

Use `build default-library`.

```bash
cofkit build default-library
```

Primary artifacts:

- `examples/default_monomers_library/README.md`
- `examples/default_monomers_library/registry.jsonl`
- `examples/default_monomers_library/failures.jsonl`

### Classify an existing batch output tree

Use `analyze classify-output`.

```bash
cofkit analyze classify-output \
  <BATCH_OUTPUT_DIR> \
  --output-dir out/classified_batch
```

Primary artifacts:

- `classification_manifest.jsonl`
- `valid/manifest.jsonl`
- `warning/manifest.jsonl`
- `hard_invalid/manifest.jsonl`
- `hard_hard_invalid/manifest.jsonl`

### Run Zeo++ pore analysis on one CIF

Use `analyze zeopp`.

```bash
export COFKIT_ZEOPP_PATH=/path/to/zeo++/network

cofkit analyze zeopp \
  <STRUCTURE_CIF> \
  --output-dir out/zeopp_run \
  --json
```

Useful variants:

- repeat `--probe-radius ...` to request accessibility-aware scans in addition to the default point-probe baseline
- add `--channel-radius ...` when channel scans should use a shared radius instead of matching each probe radius
- tune `--surface-samples-per-atom` and `--volume-samples-total` when Zeo++ Monte Carlo settings need to change
- add `--zeopp-path ...` if `COFKIT_ZEOPP_PATH` is not configured

Primary artifacts:

- `zeopp_report.json`
- raw Zeo++ output files and stdout/stderr logs under the output directory

### Run LAMMPS local optimization on one CIF

Use `calculate lammps-optimize`.

```bash
export COFKIT_LMP_PATH=/path/to/lmp_mpi

cofkit calculate lammps-optimize \
  <STRUCTURE_CIF> \
  --output-dir out/lammps_opt \
  --json
```

Useful variants:

- add `--pre-minimization-steps ...` plus the other `--pre-minimization-*` flags for a short restrained prerun before minimization
- add `--two-stage` plus `--stage2-*` for a weaker or unrestrained second minimization stage
- add `--relax-cell` plus `--box-relax-*` for a final box/relax stage
- tune `--energy-tolerance`, `--force-tolerance`, `--max-iterations`, `--max-evaluations`, `--min-style`, `--timestep`, and `--min-modify-*` as needed
- add `--lmp-path ...` if `COFKIT_LMP_PATH` is not configured

Requirements:

- input must be an explicit-bond `P1` CIF
- `_ccdc_geom_bond_type` must be present on every bond row
- the current public force-field backend is `UFF` only

Primary artifacts:

- `lammps_report.json`
- `*_lammps_optimized.cif`
- `lammps_input.data`
- `lammps_minimize.in`
- `lammps.log`, `lammps.stdout.log`, `lammps.stderr.log`
- `lammps_trajectory.lammpstrj`

### Work from Python instead of the CLI

Choose the API by how much the user already knows:

- `COFEngine`: explicit topology, project-style generation
- `BatchStructureGenerator.generate_monomer_pair_candidate(...)`: best structure for one pair
- `BatchStructureGenerator.generate_monomer_pair_candidates(...)`: all supported structures for one pair
- `BatchStructureGenerator.run_binary_bridge_batch(...)`: one template over libraries

## Output Reading Rules

- For `summary.json`, report `template_id`, input monomer IDs, attempted and successful structures, and each result's `topology_id`, `status`, `score`, and `cif_path`.
- For `manifest.jsonl`, treat each line as one generated structure summary. Use it when the printed terminal summary is not enough.
- For `summary.md`, treat it as the human-readable batch overview, not the detailed machine-readable source of truth.
- For `combined_summary.json`, summarize each template separately and do not collapse per-template counts together without saying so.
- For classification outputs, keep the four-way split explicit: `valid`, `warning`, `hard_invalid`, `hard_hard_invalid`.
- For `zeopp_report.json`, keep the point-probe baseline separate from any requested probe scans. Do not assume `probe_scans_successful` means no useful probe data exists; inspect individual scan status and parsed fields.
- For `lammps_report.json`, report the optimized CIF path, atom/bond/angle/dihedral/improper counts, the forcefield backend, the key settings used, and any warnings.

## Common Traps

- Do not treat `cofkit build list-templates` as a complete picture of the whole toolkit; it only covers the build chemistry surface.
- Do not treat every template returned by `cofkit build list-templates` as immediately runnable for structure generation.
- Do not assume an unspecified topology means one implicit default topology; the current `single-pair` and `batch-binary-bridge` CLIs enumerate all applicable topologies unless the user passes `--no-all-topologies` or explicit `--topology` values.
- Do not use `--auto-detect-libraries` on `examples/default_monomers_library`.
- Do not claim stacking workflows are supported; `stacking_mode` remains `"disabled"`.
- Do not route users to the current internal benzothiazole/sulfur-enabled conversion prototype through the public CLI.
- Do not send legacy atomistic CIFs without `_ccdc_geom_bond_type` into `cofkit calculate lammps-optimize`.
- Do not answer from docs alone when the command can be run and the artifacts can be inspected directly.
- Do not default to deprecated flat aliases such as `cofkit batch-imine` when the grouped binary-bridge interface expresses the same task more clearly.
