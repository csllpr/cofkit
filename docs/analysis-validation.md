# Analysis And Validation

## Classify A Finished Output Tree

Use the coarse classifier after a batch run:

```bash
cofkit analyze classify-output \
  out/batch_imine \
  --output-dir out/batch_imine_coarse_validation
```

By default, this writes symlink-based views of the source CIFs. Use `--link-mode copy` or `--link-mode hardlink` when needed.

The classifier writes:

- `classification_manifest.jsonl`
- `valid/manifest.jsonl` plus `valid/cifs/`
- `warning/manifest.jsonl`, `warning/cifs/`, and `warning/reasons/<reason>/`
- `needs_optimization/manifest.jsonl`, `needs_optimization/cifs/`, and `needs_optimization/reasons/<reason>/`
- `hard_hard_invalid/manifest.jsonl`, `hard_hard_invalid/cifs/`, and `hard_hard_invalid/reasons/<reason>/`
- `hard_invalid/manifest.jsonl`, `hard_invalid/cifs/`, and `hard_invalid/reasons/<reason>/`

See [COARSE_VALIDATION.md](COARSE_VALIDATION.md) for bucket definitions and thresholds.

## Decompose One CIF To COFid

```bash
cofkit analyze decompose \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif
```

When `--topology` is omitted, `cofkit` first reads an embedded `# COFid:` comment when present, then attempts conservative topology detection from the recovered periodic linkage graph. Pass `--topology` to force the topology token used in the output COFid:

```bash
cofkit analyze decompose \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif \
  --topology hcb
```

By default, this prints only the recovered COFid. Add `--json` for the full decomposition payload, including topology-detection diagnostics when auto mode is used:

```bash
cofkit analyze decompose \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif \
  --json
```

Current scope:

- input must be an atomistic CIF for one of the supported binary-bridge or ring-forming structures
- explicit `_geom_bond_*` connectivity is preferred
- if the bond loop is absent, `cofkit` falls back to periodic distance-based bond detection
- use `--bond-mode distance` to force periodic distance-based bond detection even when explicit CIF bond rows are present
- if bond labels are present but `_ccdc_geom_bond_type` / `_geom_bond_type` is absent, `cofkit` infers bond orders from local geometry
- supported canonical linkage codes are `imine`, `hydrazone`, `azine`, `boest`, `bken`, `vinylene`, `boroxine`, and `triazine`
- overlapping nitrogen-containing C=N bonds are classified per bond through two independent branches: `azine > hydrazone > imine` for N-N environments and `bken > imine` for the keto-enamine environment; azine requires `C=N-N=C`, hydrazone requires `C=N-N-C(=O)`, and beta-ketoenamine requires the keto-aldehyde-derived ring-carbonyl environment
- a C=N bond claimed by a more specific branch is excluded from generic imine decomposition; a bond matching both independent branches is reported as cross-branch ambiguous instead of being assigned by a global priority
- vinylene decomposition rejects C=C bonds in five- or six-membered rings, requires carbon anchors at both endpoints, and requires exactly one endpoint to have the stronger activated-methylene environment (a carbon anchor activated by multiple nitrogens or a multiple bond to N/O/S); the raw and rejected C=C counts are reported under `metadata["vinylene_linkage_detection"]`
- recognized boronate-ester or boroxine chemistry takes priority over a vinylene match in the same CIF; `recognized_boron_linkages`, `boron_linkage_override_applied`, and `boron_linkage_rejected_bond_count` record that resolution under `metadata["vinylene_linkage_detection"]`
- independently recoverable imine or post-boron-resolution vinylene chemistry takes priority over triazine decomposition in the same CIF; the evaluations and winning linkage families are reported under `metadata["triazine_linkage_resolution"]`, and raw C=N candidates that do not recover imine monomers do not suppress triazine
- template-id aliases such as `hydrazone_bridge`, `boronate_ester_bridge`, `boroxine_trimerization`, and `triazine_trimerization` are accepted through `--linkage`
- boroxine and triazine decomposition recognizes complete six-membered product rings, restores the boronic-acid or nitrile precursor groups, and reconstructs the net through virtual three-connected ring nodes
- equivalent disconnected layer graphs in named stacked ring-forming outputs are normalized before topology ranking
- topology can be supplied through `--topology`; otherwise auto-detection is attempted against cofkit's available topology repository
- auto-detection is conservative and may report ambiguity; special decorated routes such as `bex` may still need explicit `--topology bex`

The decomposition workflow was adapted from deCOFpose: <https://github.com/r-fedorov/deCOFpose>.

## Validate A CIF Against COFid

Simple mode forces distance-inferred decomposition and compares recovered monomer blocks and linkage to the supplied COFid. Topology is intentionally not compared in validation mode; use `analyze decompose --json` when topology-detection diagnostics are needed.

```bash
cofkit validate simple \
  '<COFID>' \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif
```

Optimize mode first runs the default LAMMPS optimization pipeline, then decomposes the optimized CIF with the same distance-inferred comparison:

```bash
cofkit validate optimize \
  '<COFID>' \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_validate_lammps
```

Use `--json` on either mode for structured diagnostics.

## Zeo++ Pore Analysis

Configure the Zeo++ `network` binary:

```bash
export COFKIT_ZEOPP_PATH=/path/to/zeo++/network
```

Run the point-probe baseline:

```bash
cofkit analyze zeopp \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_zeopp \
  --json
```

Add accessibility-aware probe scans with repeated `--probe-radius`:

```bash
cofkit analyze zeopp \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_zeopp \
  --probe-radius 1.20 \
  --probe-radius 1.86 \
  --json
```

The wrapper keeps raw Zeo++ outputs, stdout/stderr logs, and `zeopp_report.json` in the output directory.
