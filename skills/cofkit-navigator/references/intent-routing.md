# Intent Routing

This reference maps natural-language requests onto the current `cofkit` interfaces.

## Command Launchers

Use one of these equivalent launch styles:

- installed CLI: `cofkit ...`
- repo-local fallback: `PYTHONPATH=src .venv/bin/python -m cofkit.cli ...`

## Request Patterns

### Show what the toolkit can do today

Use:

```bash
cofkit build list-templates --json
```

Interpretation rules:

- Focus on templates where `supports_pair_generation` is `true` for practical structure generation.
- Registered templates with `supports_pair_generation=false` are still part of the reaction library, but they are not on the current topology-guided pair-generation path.

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

Add `--auto-detect-libraries` only when `<LIBRARY_DIR>` contains raw generic `.txt` SMILES files rather than grouped files like `amines_count_3.txt`.

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

## Common Traps

- Do not treat every template returned by `cofkit build list-templates` as immediately runnable for structure generation.
- Do not use `--auto-detect-libraries` on `examples/default_monomers_library`.
- Do not claim stacking workflows are supported; `stacking_mode` remains `"disabled"`.
- Do not answer from docs alone when the command can be run and the artifacts can be inspected directly.
- Do not default to deprecated flat aliases such as `cofkit batch-imine` when the grouped binary-bridge interface expresses the same task more clearly.
