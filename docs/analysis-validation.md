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

When `--topology` is omitted, `cofkit` validates an embedded `# COFid:` topology against the recovered periodic graph when present, then falls back to conservative graph-based topology detection if the comment is incompatible. Pass `--topology` to require that topology in the output COFid:

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

The established per-family engine remains the default. Use the experimental event/hypothesis engine explicitly for side-by-side evaluation:

```bash
cofkit analyze decompose \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif \
  --decomposition-mode event \
  --json
```

Event mode emits structured local events and globally validated reconstruction hypotheses under `metadata["event_detection"]` and `metadata["hypotheses"]`. It has not replaced legacy mode because the independent labelled benchmark set is still pending. See [event-decomposition.md](event-decomposition.md).

Current scope:

- input must be an atomistic CIF for one of the supported binary-bridge or ring-forming structures
- input atom sites must be expanded in `P1`; non-identity space-group operations and self-image bond rows are rejected explicitly rather than silently flattened
- explicit `_geom_bond_*` connectivity is preferred
- if the bond loop is absent, `cofkit` falls back to periodic distance-based bond detection
- use `--bond-mode distance` to force periodic distance-based bond detection even when explicit CIF bond rows are present
- if bond labels are present but `_ccdc_geom_bond_type` / `_geom_bond_type` is absent, `cofkit` infers bond orders from local geometry
- periodic distances are evaluated with a triclinic-safe nearest-lattice-image search, and distinct explicit bonds between the same atom labels are retained when their periodic image gains differ
- supported canonical linkage codes are `imine`, `hydrazone`, `azine`, `boest`, `bken`, `vinylene`, `boroxine`, and `triazine`
- `--linkage auto` is the default; all supported families are evaluated and success requires exactly one complete periodic precursor decomposition
- `--decomposition-mode legacy` is the default; `--decomposition-mode event` opts into the experimental event/hypothesis engine without changing legacy behavior
- overlapping nitrogen-containing C-N/C=N bonds are classified per bond through two independent branches: `azine > hydrazone > imine` for N-N environments and `bken > imine` for the keto-enamine environment; azine requires `C=N-N=C`, hydrazone requires `C=N-N-C(=O)`, and beta-ketoenamine accepts the conventional C-N/C=C keto-enamine tautomer next to the keto-aldehyde-derived ring-carbonyl environment
- a C=N bond claimed by a more specific branch is excluded from generic imine decomposition; a bond matching both independent branches is reported as cross-branch ambiguous instead of being assigned by a global priority
- vinylene decomposition rejects C=C bonds in five- or six-membered rings, requires carbon anchors at both endpoints, and requires exactly one endpoint to have the stronger activated-methylene environment (a carbon anchor activated by multiple nitrogens or a multiple bond to N/O/S); the raw and rejected C=C counts are reported under `metadata["vinylene_linkage_detection"]`
- recognized boronate-ester or boroxine chemistry is reported under `metadata["vinylene_linkage_detection"]` but does not suppress a separate structurally valid vinylene candidate
- triazine assignment is reported as ambiguous only when an imine or vinylene strategy independently recovers valid precursor motifs and a rank-2 or rank-3 periodic backbone; the evaluations and coexisting linkage families are reported under `metadata["triazine_linkage_resolution"]`
- template-id aliases such as `hydrazone_bridge`, `boronate_ester_bridge`, `boroxine_trimerization`, and `triazine_trimerization` are accepted through `--linkage`
- boroxine and triazine decomposition recognizes complete six-membered product rings, restores the boronic-acid or nitrile precursor groups, and reconstructs the net through virtual three-connected ring nodes
- equivalent disconnected layer graphs in named stacked ring-forming outputs are normalized before topology ranking
- topology can be supplied through `--topology`; otherwise auto-detection is attempted against cofkit's available topology repository
- supplied and comment-derived topologies are validated against the recovered connectivity mode and exact periodic cycle-gain rank; rank-0 molecules and rank-1 chains cannot be accepted as 2D or 3D frameworks
- auto-detection is conservative and may report ambiguity; special decorated routes such as `bex` may still need explicit `--topology bex`
- success requires the complete linkage-specific precursor set (two correctly typed blocks for binary bridges, one for ring formation), exact agreement between declared connectivity and rebuilt RDKit motif counts, and a COFid accepted by the build-input validator; partial or chemically inconsistent recovery is returned as `skipped`, while malformed CIF/topology input is returned as `error` in JSON mode

The decomposition workflow was adapted from deCOFpose: <https://github.com/r-fedorov/deCOFpose>.

## Validate A CIF Against COFid

Simple mode forces distance-inferred decomposition and compares recovered monomer blocks and linkage to the supplied COFid. The expected topology is enforced as a periodic-rank and connectivity compatibility constraint, although the validation result does not separately score topology-token equality. Use `analyze decompose --json` when detailed topology diagnostics are needed.

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
