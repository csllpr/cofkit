# Calculations

The `calculate` namespace wraps optional external tools. Build and basic validation workflows do not require these binaries.

## LAMMPS Optimization

Configure LAMMPS:

```bash
export COFKIT_LMP_PATH=/path/to/lmp
```

Run a local cleanup on one explicit-bond `P1` CIF:

```bash
cofkit calculate lammps-optimize \
  out/cli_single_pair/cifs/valid/tapb__tfb__hcb.cif \
  --output-dir out/tapb_tfb_lammps_opt \
  --forcefield dreiding \
  --json
```

Important controls:

- `--forcefield` selects a registered force-field family or stable parameter-set ID
- `--charge-model {none,eqeq}` controls charge staging; default is `eqeq`
- `--pair-cutoff`, `--coulomb-cutoff`, and `--ewald-precision` tune nonbonded settings
- `--pre-minimization-steps` and `--pre-minimization-*` tune the restrained prerun
- `--two-stage` / `--no-two-stage` and `--stage2-*` tune staged minimization
- `--timestep` and `--min-modify-*` expose LAMMPS minimizer controls
- `--relax-cell` / `--no-relax-cell` and `--box-relax-*` control the final `fix box/relax` stage
- `--timeout-seconds` caps the LAMMPS subprocess

For real optimization work, prefer `--forcefield dreiding`. `UFF` remains available for compatibility and comparison runs, but should be treated as experimental support.

## Force-Field Metadata

Framework force fields are registered as versioned parameter-set implementations rather than bare family names. Each registry entry records a stable ID and aliases, primary citations, checksummed parameter artifacts, atom-typing implementation, charge assumptions, nonbonded mixing rules, intramolecular scaling, element coverage, validation status, and parameter-data license status. The packaged registry is available through `cofkit.load_packaged_forcefield_metadata()` and `cofkit.resolve_forcefield_metadata()`.

Current default IDs are `uff-openbabel-3.1.0-cofkit-1.0` and `dreiding-standard-1990-cofkit-1.0`. Family selectors such as `uff` and `dreiding` resolve to their registered defaults; stable IDs can be used anywhere a framework force field is selected. Calculation reports embed the complete resolved metadata snapshot. Monte Carlo run directories also contain `framework_forcefield_metadata.json`, separate from guest provenance in `guest_forcefield_metadata.json`.

## EQeq

LAMMPS and Monte Carlo wrappers use EQeq for framework charge assignment when charges are enabled:

```bash
export COFKIT_EQEQ_PATH=/path/to/eqeq
```

Common EQeq controls are shared across wrappers:

- `--eqeq-path`
- `--eqeq-lambda`
- `--eqeq-h-i0`
- `--eqeq-charge-precision`
- `--eqeq-method`
- `--eqeq-real-space-cells`
- `--eqeq-reciprocal-space-cells`
- `--eqeq-eta`
- `--eqeq-timeout-seconds`

## Widom Insertion

Configure a Monte Carlo backend:

```bash
export COFKIT_GRASPA_PATH=/path/to/nvc_main.x
# or, for RASPA2:
export COFKIT_RASPA2_PATH=/path/to/simulate
```

Run Widom insertion:

```bash
cofkit calculate graspa-widom \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_widom \
  --forcefield dreiding \
  --component CO2 \
  --component N2 \
  --widom-moves-per-component 300000 \
  --json
```

Use `--backend raspa2` for RASPA2. The selected executable can be overridden with `--graspa-path`, `--raspa2-path`, or backend-neutral `--raspa-path`.

Packaged guest parameters carry explicit family, source, and framework-compatibility metadata. `TIP4P`, `CO2`, `H2`, `N2`, and `SO2` are registered in the DREIDING family and are accepted only with `--forcefield dreiding`. `Xe` and `Kr` are registered in the UFF family and are accepted with both UFF and DREIDING because DREIDING does not provide parameters for them. Their current values come from the RASPA `GenericMOFs` force field and are distinct from the bundled Open Babel `UFF.prm`; this distinction is recorded in `guest_forcefield_metadata.json`. `--all-components` therefore requires DREIDING; a UFF run can select `Xe` and/or `Kr`.

External parameterized guests use repeated `--guest-bundle path/to/guest.json` flags and are selected by bundle `name` or alias. Guest-bundle schema version 2 requires top-level `parameter_family`, `parameter_source`, and `compatible_framework_forcefields` fields; incompatible selections fail before a backend is executed. Legacy version 1 bundles must be upgraded rather than silently treated as compatible.

Outputs include staged `eqeq/` and `widom/` directories, backend logs, raw `Output/**/*.data`, `widom/Output/results.csv`, `widom/framework_forcefield_metadata.json`, `widom/guest_forcefield_metadata.json`, and `graspa_widom_report.json`.

## Single-Component Isotherms

```bash
cofkit calculate graspa-isotherm \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_isotherm \
  --forcefield dreiding \
  --component CO2 \
  --pressure 10000 \
  --pressure 100000 \
  --pressure 1000000 \
  --json
```

Use exactly one `--component NAME` and one or more repeated `--pressure PA` values. `--fugacity-coefficient` accepts a positive float or `PR-EOS`; for RASPA2, `PR-EOS` omits an explicit fugacity coefficient so RASPA2 can use its internal calculation.

Outputs include staged `eqeq/`, `isotherm/framework.cif`, one `isotherm/pressure_*/` directory per pressure point with framework and guest force-field metadata JSON files, `isotherm/results.csv`, and `graspa_isotherm_report.json`.

## Mixture Adsorption

```bash
cofkit calculate graspa-mixture \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_mixture \
  --forcefield dreiding \
  --component Kr:0.1 \
  --component Xe:0.9 \
  --pressure 10000 \
  --pressure 100000 \
  --fugacity-coefficient PR-EOS \
  --json
```

Mixtures require at least two repeated `--component NAME:FRACTION` values and one or more pressure points. The wrapper computes adsorbed mole fractions and pairwise selectivities using `(x_i / x_j) / (y_i / y_j)`.

Outputs include `mixture/component_results.csv`, `mixture/selectivity_results.csv`, staged per-pressure backend directories with framework and guest force-field metadata JSON files, and `graspa_mixture_report.json`.

## Hybrid MD/MC

```bash
cofkit calculate hybrid-mdmc \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_hybrid_mdmc \
  --cycles 5 \
  --component CO2 \
  --pressure 100000 \
  --lammps-forcefield dreiding \
  --raspa-forcefield dreiding \
  --md-steps 1000 \
  --gcmc-production-cycles 200000 \
  --json
```

The default exchange mode is `framework`: each cycle carries the MD-updated framework CIF into the next GCMC segment. Add `--exchange-mode guest-restart` to carry final gRASPA guest snapshots into the following LAMMPS MD segment and write post-MD guest coordinates back to the next gRASPA MC segment.

Guest restart currently requires `--backend graspa`; RASPA2 restartfile staging is not supported. This is an alternating MD/GCMC workflow, not dynamic GCMC inside LAMMPS.
