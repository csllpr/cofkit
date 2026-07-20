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
  --component CO2_DREIDING \
  --component N2_DREIDING \
  --widom-moves-per-component 300000 \
  --json
```

Use `--backend raspa2` for RASPA2. The selected executable can be overridden with `--graspa-path`, `--raspa2-path`, or backend-neutral `--raspa-path`.

Every packaged guest selector is force-field tagged, and the same tag is stored as validated `model_tag` metadata. The default truncated-LJ probe set is `TIP4P_DREIDING`, `CO2_DREIDING`, `H2_DREIDING`, `N2_DREIDING`, `SO2_DREIDING`, `Xe_GENERICMOFS`, and `Kr_GENERICMOFS`; `--all-components` selects this default set. Untagged packaged names are intentionally not aliases, so a selector cannot silently choose a model. The guest tag identifies its parameter model or source; it does not select or restrict the framework `--forcefield`. The user explicitly chooses that pairing, and cofkit records both selections without applying a guest/framework compatibility allowlist. The RASPA GenericMOFs Xe/Kr values remain distinct from the bundled Open Babel `UFF.prm` values.

The optional RASPA2 ExampleMoleculeForceField set provides `He_RASPA`, `Ar_RASPA`, `CH4_RASPA`, `O2_RASPA`, `CO2_RASPA`, and `N2_RASPA`. These models use shifted LJ interactions, no tail corrections, and Lorentz-Berthelot mixing, and are not selected by `--all-components`. The framework force field remains an independent, explicit user choice. Because RASPA has one global shifted/truncated setting, every component in a mixture must use the same convention. For example, use `CH4_RASPA` with `CO2_RASPA`, not with the truncated-LJ `CO2_DREIDING`:

```bash
cofkit calculate graspa-mixture framework.cif \
  --forcefield dreiding \
  --component CH4_RASPA:0.5 \
  --component CO2_RASPA:0.5 \
  --pressure 100000
```

The RASPA example models are supported by the Widom, isotherm, and mixture workflows with either gRASPA or RASPA2. They are intentionally marked unsupported for hybrid `guest_restart` exchange because cofkit's current LAMMPS handoff does not reproduce shifted LJ interactions or massless charge-only sites. Hybrid `framework` exchange remains available because guest coordinates are not injected into LAMMPS. RASPA describes these source files as examples that require validation for the target system; cofkit preserves that warning in the staged `guest_forcefield_metadata.json` provenance.

External parameterized guests use repeated `--guest-bundle path/to/guest.json` flags and are selected by bundle `name` or alias. Guest-bundle schema version 2 requires top-level `parameter_family` and `parameter_source` provenance fields. Optional `vdw_treatment` (`truncated` or `shifted`), `tail_corrections` (boolean), and `mixing_rule` (`lorentz_berthelot`) fields declare the bundle's global RASPA convention and default to cofkit's legacy truncated/no-tail/Lorentz-Berthelot behavior. The framework force field is selected independently. Legacy version 1 bundles must still be upgraded because they lack the required provenance contract.

Outputs include staged `eqeq/` and `widom/` directories, backend logs, raw `Output/**/*.data`, `widom/Output/results.csv`, `widom/framework_forcefield_metadata.json`, `widom/guest_forcefield_metadata.json`, and `graspa_widom_report.json`.

## Single-Component Isotherms

```bash
cofkit calculate graspa-isotherm \
  out/tapb_tfb_lammps_opt/tapb__tfb__hcb_lammps_optimized.cif \
  --output-dir out/tapb_tfb_isotherm \
  --forcefield dreiding \
  --component CO2_DREIDING \
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
  --component Kr_GENERICMOFS:0.1 \
  --component Xe_GENERICMOFS:0.9 \
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
  --component CO2_DREIDING \
  --pressure 100000 \
  --lammps-forcefield dreiding \
  --raspa-forcefield dreiding \
  --md-steps 1000 \
  --gcmc-production-cycles 200000 \
  --json
```

The default exchange mode is `framework`: each cycle carries the MD-updated framework CIF into the next GCMC segment. Add `--exchange-mode guest-restart` to carry final gRASPA guest snapshots into the following LAMMPS MD segment and write post-MD guest coordinates back to the next gRASPA MC segment.

Guest restart currently requires `--backend graspa`; RASPA2 restartfile staging is not supported. This is an alternating MD/GCMC workflow, not dynamic GCMC inside LAMMPS.
