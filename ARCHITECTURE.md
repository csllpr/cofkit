# cofkit architecture

## Current modules

- `cofkit.geometry`
  - vector helpers
  - local motif frames
- `cofkit.model`
  - domain dataclasses for motifs, monomers, reaction templates, reaction events, assembly state, and candidates
- `cofkit.reactions`
  - builtin reaction library plus linkage-profile registry for pair-role ordering, bridge targets, batch-library prefixes, and realization hooks
- `cofkit.product_graph`
  - periodic product graph with motif consumption tracking and validation
- `cofkit.planner`
  - topology hints, repository-backed net compatibility checks, and net planning
- `cofkit.topology_importers`
  - dependency-free RCSR CGD-style archive parsing and archive discovery/configuration
- `cofkit.topology_index`
  - topology definition and index dataclasses
- `cofkit.topology_repository`
  - index-backed topology listing, filtering, loading, and builtin fallback handling
- `cofkit.search`
  - discrete assignment plans and backtracking reaction-event assignment
- `cofkit.chem`
  - registry-backed motif metadata, lightweight molecule parsing, geometric motif detection, optional RDKit/SMARTS-backed binary-bridge monomer construction, and bridge-geometry helpers
- `cofkit.reaction_realization`
  - registry-backed atomistic bond/deletion realization for the currently implemented binary-bridge products
- `cofkit.batch_models`
  - neutral batch-facing dataclasses and compatibility aliases for pair and run summaries
- `cofkit.monomer_library`
  - extracted monomer-role resolution, detector-backed library loading, and explicit-library regrouping
- `cofkit.single_node_topologies`
  - symmetry-expanded `2D` one-node topology layouts and explicit `P1` node/edge-orbit reconstruction for shared batch and single-pair generation
- `cofkit.single_node_topologies_3d`
  - `3D` one-node topology layouts and explicit primitive-cell expansions for `dia` and `pcu`
- `cofkit.indexed_topology_layouts`
  - generic layout reconstruction from chemistry-compatible indexed topologies for more complex `node-node` and `node-linker` cases beyond the handcrafted one-node families
- `cofkit.topology_symmetry`
  - optional `gemmi`-backed space-group expansion of indexed `2D` and `3D` topologies into explicit `P1` node/edge graphs for metadata scans and future builder reuse
- `cofkit.topology_builders`
  - shared topology-family builder registry for supported one-node and indexed-topology pair generation across batch and direct single-pair flows
- `cofkit.embedding`
  - simple periodic cell and monomer-pose initialization from topology hints and motif/reaction geometry, including lower-symmetry `hcb` cells when asymmetric trigonal monomers cannot satisfy a symmetric hexagonal metric
- `cofkit.optimizer`
  - dependency-free continuous refinement of initial cell scale and rigid monomer poses
- `cofkit.scoring`
  - first-pass candidate ranking from event coverage, topology priors, bridge geometry, and unreacted penalties
- `cofkit.cif`
  - legal P1 CIF export with atomistic output when monomer coordinates exist and coarse fallback pseudo-sites otherwise
- `cofkit.engine`
  - orchestration layer that combines planning, assignment, product-graph construction, initial embedding, optimization, and candidate scoring
- `cofkit.batch`
  - batch binary-bridge generation over monomer libraries, including RDKit-backed monomer caching, detector-backed library autodiscovery, `3+3` / `3+2` / `4+4` / `4+2` / `6+2` pair handling, symmetry-expanded one-node topology builders, indexed-topology fallback, manifest output, process-level pair execution, and CIF export enabled by default
- `cofkit.validation`
  - coarse post-generation validation / triage into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid` buckets using bridge metrics plus CIF-backed network and clash checks
- `cofkit.cli`
  - installable top-level command-line router for `build`, `analyze`, and future `calculate` namespaces
- `cofkit.cli_build`
  - build-facing command registration and handlers for single-pair generation, batch workflows, template discovery, and default-library generation
- `cofkit.cli_analyze`
  - analyze-facing command registration and handlers for output classification
- `cofkit.cli_calculate`
  - reserved namespace for future external calculation-tool orchestration

## Current pipeline

1. reaction selection
2. net planning
3. assignment
4. product graph construction
5. initial embedding
6. continuous optimization
7. scoring

## Batch flow

1. load monomer libraries from `examples/batch_test_monomers/`-style explicit text inputs or detector-scanned generic `.txt` inputs
2. build/cache `MonomerSpec` objects, optionally through RDKit/SMARTS motif detection
3. enumerate supported pairings from the selected binary-bridge linkage profile and the discovered or detected role-specific libraries
4. enumerate all applicable default topologies per pair, spanning the supported one-node families plus the current curated indexed-topology subset
5. route supported one-node and indexed-topology pair cases through the shared topology-builder registry used by both batch generation and direct single-pair input
6. write `manifest.jsonl` and `summary.md`
7. run pair generation through a process pool by default in the practical batch CLIs (`8` workers unless overridden)
8. export one CIF per successful topology-specific structure by default, routing it into `cifs/valid`, `cifs/warning`, or `cifs/invalid`
9. block CIF export for `hard_hard_invalid` structures while still recording them in manifests
10. optionally reclassify a finished output tree into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid`

## CLI flow

1. `cofkit.cli` parses a grouped user-facing command surface instead of relying on one flat namespace or ad hoc example-specific argument handling.
2. `cofkit.cli_build` owns the grouped build surface under `cofkit build ...`, while `cofkit.cli_analyze` owns the grouped analysis surface under `cofkit analyze ...`.
3. Single-pair build commands route through the same `BatchStructureGenerator.generate_monomer_pair_candidate(s)` path used by the batch engine, so topology dispatch and validation/export behavior stay aligned.
4. Batch build commands delegate monomer-role resolution and generic-library autodetection to `cofkit.monomer_library` rather than keeping that logic inside the CLI or the example wrappers.
5. The legacy scripts in `examples/` are now thin wrappers over the shared CLI entry points, which reduces drift between documented workflows and the installed interface.

## Immediate next steps

1. **Motif detection**
   - broaden SMARTS/rule-based motif discovery beyond the current binary-bridge set and add better handling for mixed-role monomers
2. **Topology integration**
   - extend the current indexed-topology path from the present curated subset toward broader multi-node and higher-symmetry topology categories
3. **Reaction geometry models**
   - richer product-specific coordinate priors and multi-event consistency handling, especially for ring-forming and cyclization chemistries
4. **Continuous optimization**
   - extend the current dependency-free bridge optimizer toward richer torsion and ring-closure variables
5. **Validation / ranking**
   - tune the current `warning` / `hard_invalid` / `hard_hard_invalid` thresholds, add better planarity metrics, ring strain metrics, and stronger reaction-geometry residuals
6. **I/O**
   - move from the current legal P1 CIF export toward richer atomistic export, bond/provenance annotations, and downstream-friendly metadata

## Explicitly out of scope for the current phase

- layered stacking exploration
- interlayer registry search
- pi-stacking scoring
- any non-`disabled` stacking mode in the engine API

## Topology flow

1. `RCSRArchiveImporter` parses zip members or extracted `.cgd` files containing CGD-like `CRYSTAL` records.
2. Bundled topology data lives under `src/cofkit/data/topologies/` with extracted `2d/` and `3d/` directories plus an `index.json` metadata cache.
3. The bundled index now precomputes chemistry-facing `two_monomer_*` compatibility metadata plus lower-level `zero_linker_*` periodic-graph diagnostics, and `cofkit.topology_data` can rebuild that cache directly from the extracted `.cgd` files when the scan logic changes.
4. `TopologyRepository` prefers that bundled index by default and lazily loads full topology definitions from the extracted `.cgd` files.
5. `NetPlanner` prefilters index entries by target dimensionality and sorted monomer connectivities.
6. `cofkit.single_node_topologies` currently expands the supported one-node `2D` nets (`hcb`, `hca`, `fes`, `fxt`, plus explicit high-connectivity opt-ins) into explicit `P1` cells for batch and direct single-pair generation.
7. `cofkit.single_node_topologies_3d` currently expands the first supported one-node `3D` family (`dia`, `pcu`) into explicit primitive cells for high-connectivity generation.
8. `cofkit.topology_symmetry` now handles a much broader set of indexed `2D` and `3D` space groups for topology-chemistry scans by reconstructing full `P1` node/edge graphs through crystallographic symmetry operations.
9. `cofkit.indexed_topology_layouts` converts chemistry-compatible indexed topologies into builder-ready layout objects for the current explicit indexed-generation path.
10. `cofkit.topology_builders` chooses the current supported one-node or indexed-layout builder for each pair/topology request so batch and direct single-pair paths share the same construction logic.
11. Explicit external archive paths still work through `TopologyRepository.from_rcsr_archives(...)` and the archive discovery helpers, but they are optional advanced inputs rather than required setup.
12. If no real topology matches, builtin hints remain available as fallback.

The current bundled repository contains 194 indexed `2D` topologies and 2403 indexed `3D` topologies, with one extracted `3D` CGD file currently recorded as unsupported by the parser.

## Practical chemistry surface

The current end-to-end binary-bridge path now covers:

- `imine_bridge`
- `hydrazone_bridge`
- `keto_enamine_bridge`
- `boronate_ester_bridge`
- `vinylene_bridge`

The current batch autodetect surface is aligned to the practical RDKit motif registry:

- `amine`
- `aldehyde`
- `hydrazide`
- `boronic_acid`
- `catechol`
- `keto_aldehyde`
- `activated_methylene`

The generated example library under [`examples/default_monomers_library`](examples/default_monomers_library) is produced from the detector output and serves as a stable explicit-library snapshot of the current role-detection behavior.

## Design principles

- reaction-aware rather than connector-aware
- hypergraph-capable for ring-forming events
- ensemble output rather than single deterministic structure
- scientifically honest separation between seed assembly and physical relaxation
- keep scoring and optimization separate even when they share bridge-geometry metrics
