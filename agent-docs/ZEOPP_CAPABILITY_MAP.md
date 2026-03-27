# Zeo++ Capability Map

This note records the local Zeo++ `network` capability surface discovered during `cofkit` integration work on March 27, 2026. It is intentionally broader than the current public `cofkit analyze zeopp` wrapper.

## Binary and environment

- `cofkit` resolves Zeo++ from `COFKIT_ZEOPP_PATH`
- development machine path used during integration:
  `/Users/csllpr/software/zeo++-0.3/network`

## Important behavior notes

- The public `network` help text is incomplete relative to the real command surface in `main.cc`.
- This Zeo++ build appears to default to radii and high-accuracy handling already.
  Evidence:
  - `material.h` initializes `highAccuracy = true`
  - `arguments.cc` prints `High accuracy requested (DEF (default) settings).`
- The current public `cofkit` wrapper should not assume every listed Zeo++ mode is equally stable on generated CIFs.

## Publicly documented commands from Zeo++

### Format and structure exports

- `-cssr`
- `-cif`
- `-v1`
- `-xyz`
- `-superxyz`
- `-vtk`
- `-vis`
- `-nt2`
- `-mopac`
- `-supermopac`

### Basic pore metrics

- `-res`
- `-resex`
- `-chan`
- `-axs`

### Surface, volume, and distributions

- `-sa`
- `-vol`
- `-volpo`
- `-psd`
- `-block`

### Visualization and auxiliary outputs

- `-zvis`
- `-visVoro`
- `-sphericalSubstructures`
- `-findTetrahedra`
- `-cellmulti`
- `-gridBOV`
- `-gridG`
- `-gridGBohr`

### Control flags

- `-r`
- `-nor`
- `-mass`

## Additional commands present in source

These are implemented in the local Zeo++ source tree even though some are not shown in the printed help.

### Structure and metal-site diagnostics

- `-strinfo`
- `-strinfoex`
- `-oms`
- `-omsex`

### Rich pore / accessibility workflows

- `-poreinfo`
- `-poreinfoSummary`
- `-saex`
- `-zsa`
- `-vsa`
- `-lsa`
- `-zvol`
- `-vvol`
- `-lvol`
- `-zpsd`
- `-vpsd`
- `-zchan`
- `-zvor`

### Ray-based analysis

- `-ray_atom`
- `-zray_atom`
- `-ray_sphere`
- `-zray_sphere`
- `-ray_node`
- `-zray_node`
- `-ray_andrew_sphere`
- `-zray_andrew_sphere`
- `-ray_andrew_atom`
- `-zray_andrew_atom`

### Grids and projected data

- `-gridGAI`
- `-gridGAIBohr`
- `-gridGprojdata`
- `-gridGprojdataBohr`
- `-gridGprojdataperframe`
- `-gridGprojdataperframeBohr`

### Hologram / segment / feature paths

- `-holo`
- `-featholo`
- `-realfeatholo`
- `-zseg`
- `-zfeat`

### Miscellaneous internal or specialized flags

- `-ha`
- `-noha`
- `-nomass`
- `-allowAdjustCoordsAndCell`
- `-stripatomnames`
- `-cage`
- `-zcage`
- `-defining`
- `-simpl`
- `-sub`
- `-fsub`
- `-fsubM`

## Output shapes worth knowing

### `-res`

One line with:

- largest included sphere
- largest free sphere
- largest included sphere along the max free-sphere path

### `-resex`

Extends `-res` with axis-resolved values:

- free sphere along `a`, `b`, `c`
- included sphere along the best free path in `a`, `b`, `c`

### `-chan`

Parse-friendly text containing:

- number of channels
- channel dimensionality
- per-channel `Di / Df / Dif`
- summary max values
- probe radius / diameter

### `-axs`

One boolean line per Voronoi node:

- `true` for accessible
- `false` for inaccessible

### `-sa`

Compact summary lines containing:

- unit-cell volume
- density
- accessible surface area
- non-accessible surface area
- normalized surface metrics
- per-channel and per-pocket surface-area totals

### `-vol`

Compact summary lines containing:

- unit-cell volume
- density
- accessible pore volume
- non-accessible pore volume
- accessible volume fraction
- `cm^3/g` normalized volume
- per-channel and per-pocket volume totals

### `-psd`

Histogram output with:

- bin size
- bin count
- accessible sample count
- cumulative and derivative distributions

## Practical stability observations from local probing

On current `cofkit` CIF outputs:

- `-res`, `-resex`, `-chan`, and `-axs` behaved cleanly
- `-sa 0 0 ...` and `-vol 0 0 ...` also behaved cleanly and are useful as point-probe baselines
- positive-probe `-sa`, `-vol`, `-psd`, `-ray_atom`, and `-block` were much less stable on some generated CIFs and could abort with Voronoi volume-check failures

This is why the current public wrapper should remain selective and preserve raw logs.

## Recommended `cofkit` exposure levels

- `basic`
  `-res`, `-resex`, `-chan 0`, `-sa 0 0`, `-vol 0 0`
- `probe-scan`
  repeated `-chan`, `-sa`, `-vol`, `-axs` for user-specified probe radii
- `experimental`
  positive-probe Monte Carlo distributions and richer pore analysis
- `visual`
  ZeoVis / Visit / grid / Voronoi exports
- `expert`
  hidden diagnostic and zeolite-specific commands
