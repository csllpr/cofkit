# Building COFs

`cofkit` supports binary-bridge COFs generated from one role-specific monomer on each side and topology-creating cyclotrimerizations from a single precursor. The practical binary-bridge templates are:

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

The `first.geometry` and `second.geometry` objects in `summary.json` record the RDKit embedding method, whether a fallback was used, and force-field optimization status. Standard ETKDGv3 remains the primary path. If it fails, cofkit tries bounded random-coordinate ETKDGv3 and ETKDGv2 embeddings. A final 2D coordinate fallback is restricted to planar metal-containing or formally charged precursors, is not force-field minimized, and is labeled explicitly in the summary.

## Ring-Forming Linkages

Use `ring-forming` for boroxine or triazine cyclotrimerization. Three-connected product rings are represented as virtual topology nodes; real ditopic precursor molecules occupy the hcb edges. A primitive hcb build therefore contains three precursor instances and two three-participant ring events.

```bash
cofkit build ring-forming \
  --template-id boroxine_trimerization \
  --smiles 'OB(O)c1ccc(B(O)O)cc1' \
  --topology hcb \
  --output-dir out/cof1
```

For CTF-1-like chemistry:

```bash
cofkit build ring-forming \
  --template-id triazine_trimerization \
  --smiles 'N#Cc1ccc(C#N)cc1' \
  --topology hcb \
  --output-dir out/ctf1
```

The command writes `summary.json` and an atomistic CIF. Ring result rows use the same success contract as single-pair builds (`status: "ok"`), and the report includes attempted, successful, and CIF-written counters. Ring validation reports radial, angular, and planarity residuals. Boroxine realization removes three waters per ring event; triazine realization preserves atoms and rewrites the three nitrile bonds into alternating C–N ring bonds. One-precursor COFids using linkage codes `boroxine` and `triazine` are accepted through `--cofid`.

Request explicit bilayer registries with repeatable `--stacking` options:

```bash
cofkit build ring-forming \
  --template-id triazine_trimerization \
  --smiles 'N#Cc1ccc(C#N)cc1' \
  --stacking AA \
  --stacking AB \
  --output-dir out/ctf1_stacked
```

For hexagonal cells, `AA` uses no lateral shift, `AB` uses `(1/3, 1/3)`, and `slipped` uses `(1/2, 0)` in the fractional in-plane basis. Each output contains two explicit layers, duplicates the periodic reaction graph and atomistic ring realization, sets a product-aware bilayer `c` length, and embeds `stacking=<registry>` in the CIF COFid comment. Ring validation remains layer-aware.

The Python builder also supports a trigonal precursor on one hcb sublattice and indexed mixed-connectivity nets whose other nodes are 3-connected virtual rings; `kgd` is the supported ideal `3+6` example. Less compatible net/precursor geometries are emitted with `ring_geometry_rejected` rather than accepted silently.

## Stacking Variants

For eligible binary-bridge `2D` outputs, request named bilayer registries during export:

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

The command resets its output directory before writing. As a safety guard, it refuses to delete a non-empty directory that was not created by a previous run of this command (previous runs leave a `.cofkit-default-library` marker file); pass `--force` to replace such a directory anyway.

## Outputs

Single-pair runs write:

- `summary.json`
- exported CIFs under `cifs/valid`, `cifs/warning`, `cifs/needs_optimization`, or `cifs/hard_invalid`

Ring-forming runs write `summary.json` and, unless disabled, `ring-candidate-1.cif` directly under the selected output directory.

Batch runs write:

- `manifest.jsonl`
- `summary.md`
- exported CIFs under `cifs/valid`, `cifs/warning`, `cifs/needs_optimization`, or `cifs/hard_invalid`

Structures classified as `hard_hard_invalid` are recorded in the manifest, but CIF export is blocked and `cif_export_blocked = true` is set in the per-structure metadata.
