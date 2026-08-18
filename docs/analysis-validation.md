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

The event/hypothesis engine is the default. Use the retained per-family engine explicitly when compatibility or side-by-side evaluation is needed:

```bash
cofkit analyze decompose \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif \
  --decomposition-mode legacy \
  --json
```

Default event mode emits structured local events and globally validated reconstruction hypotheses under `metadata["event_detection"]` and `metadata["hypotheses"]`. See [event-decomposition.md](event-decomposition.md).

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
- `--decomposition-mode event` is the default; `--decomposition-mode legacy` selects the retained per-family compatibility engine
- topology reconstruction follows chemical event semantics: the two B-O cuts in one boronate-ester ring contribute one quotient-graph edge, while multi-cut events such as azine retain distinct precursor connections; inconsistent periodic gains within an equivalent-bond event are rejected
- default event mode resolves nitrogen overlaps locally, including conventional and keto-tautomerized azine/acylhydrazone N-N products before an overlapping beta-ketoenamine interpretation; the retained legacy engine exposes its independent `azine > hydrazone > imine` and `bken > imine` per-bond branches under `metadata["nitrogen_linkage_detection"]`
- vinylene decomposition rejects C=C bonds in five- or six-membered rings, requires carbon anchors at both endpoints, and supports activated methyl groups through direct withdrawing groups or conjugated aza/cyano aromatic systems; event mode branches tied orientations for global validation
- recognized boronate-ester or boroxine chemistry is reported under `metadata["vinylene_linkage_detection"]` but does not suppress a separate structurally valid vinylene candidate
- default event mode treats a detected triazine ring as a monomer motif when another complete supported-family reconstruction retains it intact; explicit legacy mode retains its narrower imine/vinylene competition rule
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
