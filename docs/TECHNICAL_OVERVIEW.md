# Technical Overview

## Design notes

The long-term direction is:

1. reactant-aware input models
2. reaction grammar / template library
3. discrete assignment over motifs and reaction events
4. periodic product graph construction
5. initial embedding followed by a lightweight continuous optimization pass over cell scale, rigid poses, and bridge geometry
6. ranked candidate ensembles instead of a single structure

For the current phase, stacking exploration is intentionally out of scope. `stacking_mode` must stay `"disabled"`; the engine now rejects other values rather than implying partial support.

## Topology repository

The topology layer now has three levels:

- `RCSRArchiveImporter` parses CGD-style topology records from RCSR zip bundles.
- `TopologyRepository` stores an index of available topologies and loads full definitions by id.
- `NetPlanner` filters repository index entries by dimensionality and monomer connectivity before it builds plans.

Each index entry exposes cheap metadata without reparsing the full topology:

- topology id / name
- dimensionality
- number of node definitions
- number of edge definitions
- expected connectivity for each node definition
- space group, cell parameters, and source archive/member when present
- precomputed `two_monomer_*` chemistry metadata in the bundled workspace index
- lower-level `zero_linker_*` periodic-graph diagnostics in the bundled workspace index

Builtins (`hcb`, `sql`, `kgm`) remain available as fallback hints when neither bundled data nor archive-backed entries provide a requested topology.

## Current pipeline

The engine currently runs:

1. reaction selection
2. net planning
3. assignment
4. product graph construction
5. initial embedding
6. dependency-free continuous optimization
7. scoring
8. optional post-generation coarse validation / triage over exported CIFs

The batch binary-bridge pipeline builds on top of that:

1. monomer-library loading
2. RDKit monomer construction / caching
3. pair enumeration from registered binary-bridge role prefixes and monomer connectivities
4. topology-family-aware candidate generation
5. manifest / summary writing
6. process-level pair execution by default in the practical batch CLIs (`8` workers unless overridden)
7. CIF export into `cifs/valid`, `cifs/warning`, or `cifs/invalid`
8. block CIF export for `hard_hard_invalid` structures while still recording them in the manifest
9. optional reclassification of a finished output tree with `examples/classify_batch_output.py`

For convenience, the documented imine workflow still uses `run_imine_batch(...)` and [`../examples/run_batch_imine_generation.py`](../examples/run_batch_imine_generation.py), but the shared implementation now lives under the generic binary-bridge path exposed by [`../examples/run_binary_bridge_generation.py`](../examples/run_binary_bridge_generation.py).

## Supported topology-family notes

For the current single-node topology families, `cofkit` expands supported CGD nets into explicit `P1` node/edge orbits before construction. In practice that means:

- `3+2` runs can currently target `hcb`, `hca`, `fes`, and `fxt`
- `3+3` runs can currently target `hcb`, `fes`, and `fxt`
- `4+2` and `4+4` runs default to `dia` in the batch workflow and are available directly through `COFEngine` with `target_dimensionality="3D"` or `target_topologies=("dia",)`
- `6+2` runs default to `pcu` in the batch workflow and are available directly through `COFEngine` with `target_dimensionality="3D"` or `target_topologies=("pcu",)`
- explicit `sql`, `kgm`, `htb`, and `hxl` requests still work in batch mode, but they are no longer the default route for the uploaded `4`- and `6`-connected fixture libraries
- `BatchStructureGenerator.generate_monomer_pair_candidate(s)` now exposes the full supported one-node family directly for single-pair imine generation, and `COFEngine.run(...)` reaches the same builders for explicit one-node topology requests plus supported `3D` single-pair defaults such as `dia` / `pcu`
- chemistry-compatible indexed topologies from the bundled repository can now also be requested explicitly in batch and direct single-pair flows, and the default selector includes a curated subset of those indexed nets when the current chemistry metadata and builder support agree

The optimizer is intentionally modest. It only refines bridge-forming candidates through lateral cell scaling, monomer translation updates, and lightweight orientation cleanup when those steps reduce bridge-geometry residuals. Ring-forming reactions stay coarse, and stacking remains completely out of scope.

For the current imine realization path, two geometry details are now important enough to treat as part of the documented behavior:

- template-specific imine motif-origin correction is applied in the supported `3D` builder paths as well as the earlier `2D` paths, so high-connectivity `dia`-style builds do not silently bypass the bent-linkage span correction
- periodic-image bridge events store realized atom overrides back in the base monomer-local frame before CIF export, which avoids pathological retained-hydrogen directions on image-crossing imine events

The CIF exporter is deliberately honest as well: if a `MonomerSpec` carries atom coordinates, it writes atomistic sites; if not, it falls back to a legal coarse CIF built from monomer centers and motif origins so the current assembly can still be inspected.
