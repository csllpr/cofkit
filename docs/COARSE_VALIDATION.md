# Coarse Validation

The optional coarse validator is meant to separate clearly broken outputs from inspectable but strained ones. It uses five classes:

- `valid`: no warning or hard-invalid criteria triggered
- `warning`: no hard-invalid criteria triggered, but at least one soft bridge-geometry criterion triggered
- `needs_optimization`: the graph appears intact, but bridge-distance geometry is outside hard thresholds and is expected to need local optimization before use
- `hard_invalid`: at least one clearly broken-network, impossible-geometry, or clash criterion triggered
- `hard_hard_invalid`: bridge geometry is so extreme that CIF export is blocked during generation

## Warning thresholds

- any bridge distance residual `> 0.75 A`
- mean bridge distance residual `> 0.35 A`
- more than `25%` of bridge events with residual `> 0.50 A`

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

## Hard-hard-invalid threshold

- any actual bridge distance `>= 2.5 A`

Hydrogen atoms are ignored in the clash check. Warning-level bridge drift still goes through CIF-backed checks, so a candidate can be promoted from `warning` to `hard_invalid` if the exported structure shows a broken network or impossible heavy-atom contact. During generation, `valid` / `warning` / `needs_optimization` / `hard_invalid` CIFs are written into separate subdirectories under `cifs/`, while `hard_hard_invalid` structures stay manifest-only with `cif_export_blocked = true`.
