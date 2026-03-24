# Wrapper Scripts

The example scripts in `examples/` are still available as thin wrappers over the shared CLI.

## Imine batch wrapper

```bash
python3 examples/run_batch_imine_generation.py \
  --input-dir examples/batch_test_monomers \
  --output-dir out/batch_imine_generation_full_cif \
  --max-cif-exports 20000 \
  --num-conformers 2
```

That workflow writes a `manifest.jsonl`, a `summary.md`, and by default one CIF per successful topology-specific structure. Use `--no-write-cif` to disable export or `--no-all-topologies` to keep only the best topology per pair. The example script is now a thin wrapper over the shared `cofkit` CLI.

## Generic binary-bridge wrapper

```bash
python3 examples/run_binary_bridge_generation.py \
  --template-id imine_bridge \
  --input-dir examples/batch_test_monomers \
  --output-dir out/binary_bridge_generation_full_cif \
  --max-cif-exports 20000 \
  --num-conformers 2
```

With `--template-id imine_bridge`, that script reproduces the same topology-specific outputs as the imine wrapper while using the template-driven registry path. The installable `cofkit` CLI remains the primary user-facing interface.

The generic script can also infer monomer roles directly from raw generic SMILES libraries when the filenames do not already encode role information:

```bash
python3 examples/run_binary_bridge_generation.py \
  --template-id keto_enamine_bridge \
  --input-dir /path/to/raw_generic_smiles_libraries \
  --output-dir out/binary_bridge_generation_keto \
  --auto-detect-libraries \
  --num-conformers 2
```

Use `--auto-detect-libraries` only on raw generic inputs whose non-empty rows all resolve to one of the selected template's supported motif kinds. Do not combine it with [`../examples/default_monomers_library`](../examples/default_monomers_library), which is already detector-grouped explicit library output.

## Rebuild the bundled example library

```bash
python3 examples/build_default_monomers_library.py
```

That script scans [`../examples/batch_test_monomers`](../examples/batch_test_monomers), writes regrouped role/count libraries under [`../examples/default_monomers_library`](../examples/default_monomers_library), and records detector provenance in `registry.jsonl` plus failures/ambiguities in `failures.jsonl`.

## Experimental smoke examples

For literature-style single-pair smoke cases across the currently implemented binary-bridge chemistries, see:

```bash
python3 examples/run_experimental_examples.py
```

To run all currently available binary-bridge batches from the detector-scanned default library with the default `8`-worker pair pool through the wrapper:

```bash
python3 examples/run_available_binary_bridge_batches.py \
  --input-dir examples/default_monomers_library \
  --output-dir out/available_binary_bridge_batches
```

The shipped [`../examples/default_monomers_library`](../examples/default_monomers_library) snapshot currently exposes `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge` batches. `boronate_ester_bridge` and `vinylene_bridge` remain covered by the direct single-pair smoke examples above.

## Output classification wrapper

To classify a finished batch output into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid` buckets through the wrapper:

```bash
python3 examples/classify_batch_output.py \
  out/full_cif_generation_default_selector_20260320 \
  --output-dir out/full_cif_generation_default_selector_20260320_coarse_validation_triage
```

That workflow writes one classification manifest plus four CIF trees:

- `valid/`
- `warning/`
- `hard_hard_invalid/`
- `hard_invalid/`

The classifier keeps the source CIFs untouched and materializes the categorized views through symlinks by default.
