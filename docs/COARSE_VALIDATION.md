# Coarse Validation

The optional coarse validator is meant to separate clearly broken outputs from inspectable but strained ones. It uses six classes:

- `valid`: no warning or hard-invalid criteria triggered and every required check was measured
- `warning`: no hard-invalid criteria triggered, but at least one soft bridge-geometry criterion triggered
- `needs_optimization`: the graph appears intact, but bridge-distance geometry is outside hard thresholds and is expected to need local optimization before use
- `hard_invalid`: at least one clearly broken-network, impossible-geometry, or clash criterion triggered
- `hard_hard_invalid`: bridge geometry is so extreme that CIF export is blocked during generation
- `unvalidated`: no hard failure was found, but at least one required check could not be measured (`is_valid` is null, never counted as valid); the missing checks are listed under `unmeasured_required_checks`

## Final-geometry measurement

Bridge-linkage verdicts are recomputed from the final exported CIF: the inter-monomer linkage bond distances are measured from the atom-site coordinates and the explicit periodic image codes in the `_geom_bond` loop (the `_geom_bond_distance` column is not trusted), and residuals / ratios are checked against the linkage-profile target distance (or the per-template realized-bond distance window). The assembly-time ("seed") bridge metrics from `score_metadata.bridge_event_metrics` are recorded only as informational `seed_*` metrics (`bridge_metrics_source: "final_cif"` marks the authoritative channel). A repair pass that moves atoms — the soft-relax pass or the LAMMPS geometry-repair route — is therefore judged on the geometry it actually produced, and a distorted repaired linkage cannot inherit a passing seed verdict.

Every validation record carries a per-check `coverage` map with one of these statuses per check:

- `measured`: the check ran and produced values
- `no_contacts`: the contact search ran to completion and found no assessable neighbor inside its search radius (distinct from missing data)
- `missing_data`: a required check could not be evaluated (missing/unparsable CIF, no bond loop to measure, or no linkage target/window for the template)
- `not_applicable`: the check does not apply to this record (e.g. no bridge events claimed and no inter-monomer bonds present)
- `skipped`: the metadata-dependent check was deliberately skipped because metadata-level hard reasons already determine the verdict

## Warning thresholds

- any bridge distance residual `> 0.75 A`
- mean bridge distance residual `> 0.35 A`
- more than `25%` of measured bridge bonds with residual `> 0.50 A`

## Hard-invalid thresholds

- any unreacted motifs
- missing or unparsable CIF
- disconnected monomer-instance graph reconstructed from CIF bonding
- any nonbonded heavy-atom contact `< 1.05 A`
- `2D` cell area `< 10.0 A^2`
- `3D` cell volume `< 20.0 A^3`

## Needs-optimization thresholds

These graph-intact bridge-geometry failures are classified as `needs_optimization` rather than `hard_invalid`:

- any bridge distance residual `> 1.00 A`
- mean bridge distance residual `> 0.60 A`
- any bridge `actual_distance / target_distance < 0.70`
- any bridge `actual_distance / target_distance > 1.60`
- any realized inter-monomer linkage bond outside its per-template distance window (currently `boronate_ester_bridge`: B-O bonds outside `1.25-1.65 A`; the ring-closure realization targets ~1.44 A)

## Hard-hard-invalid threshold

- any measured inter-monomer bridge bond distance `>= 2.5 A`

Hydrogen atoms are ignored in the clash check. Warning-level bridge drift still goes through CIF-backed checks, so a candidate can be promoted from `warning` to `hard_invalid` if the exported structure shows a broken network or impossible heavy-atom contact. During generation, `valid` / `warning` / `needs_optimization` / `hard_invalid` CIFs are written into separate subdirectories under `cifs/`, while `hard_hard_invalid` structures stay manifest-only with `cif_export_blocked = true`.
