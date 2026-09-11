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
- `--dreiding-hbond` / `--no-dreiding-hbond` toggles DREIDING hydrogen-bond export; default is enabled
- `--pair-cutoff`, `--coulomb-cutoff`, and `--ewald-precision` tune nonbonded settings
- `--pre-minimization-steps` and `--pre-minimization-*` tune the restrained prerun
- `--two-stage` / `--no-two-stage` and `--stage2-*` tune staged minimization
- `--timestep` and `--min-modify-*` expose LAMMPS minimizer controls
- `--relax-cell` / `--no-relax-cell` and `--box-relax-*` control the final `fix box/relax` stage
- `--timeout-seconds` caps the LAMMPS subprocess

`dreiding` is the shared CLI/API optimization default. `UFF` remains available for compatibility and comparison runs, but should be treated as experimental support.

Optimization exports check final-stage convergence evidence from `lammps.log`, finite global force norm within the requested force tolerance, and (for cell relaxation) virial stress within `--pressure-tolerance` (default 1 atm) for the relaxed components. The default energy tolerance is zero, enabling force-only stopping; stage overrides preserve zero. These are numerical acceptance criteria, not evidence of force-field predictive accuracy. A LAMMPS run that completes but fails these checks — most commonly stopping on the energy criterion or the iteration cap with forces still above tolerance — no longer hard-stops: the optimized CIF and `lammps_report.json` are still written, but the result is explicitly marked unconverged (`converged: false` in `lammps_report.json` and `convergence.json`, a `warnings` entry, a `# cofkit warning:` comment in the optimized CIF, and a stderr `warning:` line from the CLI). Runs whose convergence diagnostics are unreadable (missing stage markers or non-finite metrics) still raise an error without publishing an optimized CIF, since the final geometry cannot be trusted. Trajectories default to every 100 steps and always include a final snapshot with image flags. This interval limits intermediate I/O; it is not an MD sampling recommendation. Parsing retains one frame in memory.

`validate optimize` uses the same defaults and accepts `--settings-json path.json`, a JSON object of `LammpsOptimizationSettings` overrides. The Python wrapper accepts `settings=LammpsOptimizationSettings(...)`. Periodic self-bonds or multiple bonds between the same atom IDs require an explicit larger supercell; unsupported primitive representations are rejected.

With the DREIDING backend, hydrogen-bond terms follow the original paper's Table V convention (Mayo et al., *J. Phys. Chem.* 1990, 94, 8897-8909) and are exported for the LAMMPS backend only: hydrogens bonded to N/O/F-typed atoms are retyped `H__HB` (near-zero Lennard-Jones well depth), and donor-acceptor interactions are computed with a 12-10 `hbond/dreiding/lj` term (`pair_style hybrid/overlay`), with `Rhb = 2.75` angstrom, a `cos^4(theta_DHA)` angular factor, cosine periodicity 4, inner cutoff 9.0 angstrom, outer cutoff 11.0 angstrom, and an angle cutoff of 90 degrees (the LAMMPS documentation example values for this style). `Dhb` is 7.0 kcal/mol when the run carries charges (input-CIF or EQeq, following the paper's Gasteiger-like charge convention) and 9.0 kcal/mol for charge-free runs. The export is enabled by default; `--no-dreiding-hbond` opts out on `cofkit calculate lammps-optimize`, `cofkit calculate hybrid-mdmc`, and `cofkit validate optimize` (`--dreiding-hbond` / `--no-dreiding-hbond`). The flag has no effect with the UFF backend, and the RASPA/gRASPA framework assets have no hydrogen-bond term. This changes optimization numerics for structures containing N-H, O-H, or F-H bonds relative to earlier cofkit releases.

The DREIDING Cu/Ni/Mg rows are UFF parameter substitutions (Rappe et al., *J. Am. Chem. Soc.* 1992, 114, 10024-10035, DOI 10.1021/ja00051a040, transcribed from the bundled pinned `UFF.prm`) because those elements are outside the DREIDING paper coverage; DREIDING-organic plus UFF-metal mixing is standard practice but remains a model choice. A numerically converged calculation does not by itself validate the underlying force-field model.

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

EQeq now defaults to six charge decimals and verifies complete finite charges, an atom bijection by identity/elements/periodic coordinates, and intended total charge. `--eqeq-target-charge` defaults to 0 electrons; `--eqeq-net-charge-tolerance` defaults to 0.001 electrons. The latter is a configurable numerical charge-conservation budget, not a claim about EQeq charge accuracy. For larger frameworks, increase precision if rounding exceeds the budget. The wrapper never silently neutralizes charges or maps ambiguous atoms by row order. Fractional occupancies require an explicitly resolved input configuration.

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

Every packaged guest selector is force-field tagged, and the same tag is stored as validated `model_tag` metadata. The default truncated-LJ probe set is `TIP4P_DREIDING`, `CO2_DREIDING`, `N2_DREIDING`, `SO2_DREIDING`, `Xe_GENERICMOFS`, and `Kr_GENERICMOFS`; `--all-components` selects this default set. H2 is excluded because upstream gRASPA does not implement its requested Feynman-Hibbs potential; explicitly selecting `H2_DREIDING` with gRASPA fails before execution. Select a supporting backend explicitly for that model. Untagged packaged names are intentionally not aliases, so a selector cannot silently choose a model. The guest tag identifies its parameter model or source; it does not select or restrict the framework `--forcefield`. The user explicitly chooses that pairing, and cofkit records both selections without applying a guest/framework compatibility allowlist. The RASPA GenericMOFs Xe/Kr values remain distinct from the bundled Open Babel `UFF.prm` values.

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

Use exactly one `--component NAME` and one or more repeated `--pressure PA` values. `--fugacity-coefficient` defaults to `PR-EOS` and accepts a positive float or `PR-EOS`; for RASPA2, `PR-EOS` omits an explicit fugacity coefficient so RASPA2 can use its internal calculation.

Outputs include staged `eqeq/`, `isotherm/framework.cif`, one `isotherm/pressure_*/` directory per pressure point with framework and guest force-field metadata JSON files, `isotherm/results.csv`, and `graspa_isotherm_report.json`.

Adsorption calculations use uncapped native gRASPA cycles by default; Widom retains a deliberate one-insertion-move budget. Nominal cycle counts are not interchangeable between backends and do not certify equilibration. PR-EOS requires valid critical constants in the selected model; missing data is an error, with explicit coefficients available as an alternative. At zero pressure the ideal limit is rendered explicitly. Independent pressure points use reproducible child seeds; the rendered inputs record the effective seeds. RASPA2 now receives its requested seed.

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

Mixtures require at least two repeated `--component NAME:FRACTION` values and one or more pressure points. The wrapper computes adsorbed mole fractions and pairwise selectivities using `(x_i / x_j) / (y_i / y_j)`. Selectivity uncertainty is reported as unavailable (`null` in JSON) because the current backend summaries lack paired loading statistics and a common uncertainty convention. Marginal loading errors are not treated as independent, and zero observed uptake is not assigned zero uncertainty.

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

## Calculation attempts and compatibility

LAMMPS, EQeq, adsorption, and hybrid calculations claim fresh attempt directories. If an output directory already exists, the new result is written under `attempt-<unique-id>/`; use the returned/report paths to locate it. Existing files are preserved. Reusing a directory does not resume a calculation. `attempt.json` identifies the input, settings, and executable (for direct engine calls); `execution.json` hashes staged model/input files and records the command. Reports and optimized structures are published atomically. An incomplete manifest must not be treated as a completed result. Hybrid cycles derive distinct MD/MC seeds from the configured MD velocity seed and record actual per-stage settings.

Zeo++ calls now always request high-accuracy decomposition (`-ha`), recorded in the report. Probe definitions and sample counts remain explicit; production precision still requires a convergence study.

RDKit geometry preparation accepts `optimization_max_iterations` (default 500) and `optimization_attempts` (default 1) in its Python API. Only a converged minimizer return code can yield `optimized` status. Minimization that does not fully converge is downgraded to a warning rather than a hard stop: the build proceeds with the lowest-energy unconverged conformer (or the unminimized embedded conformer when no finite energy is available), labeled `unconverged` with diagnostics, mirroring the unconverged-warning behavior of LAMMPS optimizations. Deliberately unoptimized fallbacks remain labeled `skipped`.
