# Building COFs

`cofkit` currently focuses on binary-bridge COFs generated from one role-specific monomer on each side. The practical binary-bridge templates are:

- `imine_bridge`
- `hydrazone_bridge`
- `azine_bridge`
- `keto_enamine_bridge`
- `boronate_ester_bridge`
- `vinylene_bridge`

Discover the current registry with:

```bash
cofkit build list-templates
```

## Single Pair

Use `single-pair` when you have two SMILES strings:

```bash
cofkit build single-pair \
  --template-id imine_bridge \
  --first-smiles '<SMILES_A>' \
  --second-smiles '<SMILES_B>' \
  --first-id monomer_a \
  --second-id monomer_b \
  --output-dir out/single_pair
```

By default, the CLI autodetects motif kinds with `--auto-detect-motifs`, enumerates all applicable topologies with `--all-topologies`, and writes CIFs. Useful options:

- `--first-motif-kind` and `--second-motif-kind` force role assignment
- `--topology hcb` restricts topology selection; repeat for multiple topologies
- `--target-dimensionality 2D` or `--target-dimensionality 3D` controls topology selection
- `--no-all-topologies` keeps only the best generated topology
- `--no-write-cif` skips CIF export
- `--max-cif-exports N` caps exported CIFs

You can also build from a COFid:

```bash
cofkit build single-pair \
  --cofid '<COFID>' \
  --output-dir out/from_cofid
```

`--cofid` defines the monomers, topology, and linkage, so it cannot be combined with direct SMILES or topology overrides.

## Stacking Variants

For eligible `2D` outputs, request named bilayer registries during export:

```bash
cofkit build single-pair \
  --template-id imine_bridge \
  --first-smiles '<SMILES_A>' \
  --second-smiles '<SMILES_B>' \
  --topology hcb \
  --stacking AA \
  --stacking AB \
  --output-dir out/single_pair_stacked
```

This keeps the build pipeline unchanged and emits one exported structure per requested registry. CIF comment lines append the registry tag, for example `# COFid: ... stacking=AA`.

## One Linkage Over A Library

Run one binary-bridge template over a grouped monomer-library directory:

```bash
cofkit build batch-binary-bridge \
  --template-id imine_bridge \
  --input-dir examples/default_monomers_library \
  --output-dir out/batch_imine \
  --max-workers 8
```

Useful batch flags:

- `--max-pairs N` caps attempted monomer pairs
- `--no-write-cif` skips CIF export
- `--max-cif-exports N` caps exported CIFs
- `--no-all-topologies` keeps only the best topology per pair
- `--topology ID` restricts topology selection; repeat for multiple ids
- `--auto-detect-libraries` infers roles/connectivities from raw SMILES libraries

Grouped library files are named by role and connectivity, for example:

```text
amines_count_3.txt
aldehydes_count_2.txt
hydrazides_count_2.txt
keto_aldehydes_count_4.txt
```

Each file contains one SMILES string per line. A leading `smiles` header is allowed.

## All Available Binary Bridges

Run every binary-bridge template that the input library can satisfy:

```bash
cofkit build batch-all-binary-bridges \
  --input-dir examples/default_monomers_library \
  --output-dir out/available_binary_bridge_batches \
  --max-workers 8
```

With the shipped default library snapshot, this currently discovers `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge`. Other registered templates need their expected role libraries, such as `hydrazines`, `boronic_acids`, `catechols`, or `activated_methylenes`.

## Default Library

Regenerate the detector-scanned example library:

```bash
cofkit build default-library
```

The generated library under `examples/default_monomers_library` contains grouped `*_count_N.txt` files, `registry.jsonl`, and `failures.jsonl`.

## Outputs

Single-pair runs write:

- `summary.json`
- exported CIFs under `cifs/valid`, `cifs/warning`, `cifs/needs_optimization`, or `cifs/hard_invalid`

Batch runs write:

- `manifest.jsonl`
- `summary.md`
- exported CIFs under `cifs/valid`, `cifs/warning`, `cifs/needs_optimization`, or `cifs/hard_invalid`

Structures classified as `hard_hard_invalid` are recorded in the manifest, but CIF export is blocked and `cif_export_blocked = true` is set in the per-structure metadata.
