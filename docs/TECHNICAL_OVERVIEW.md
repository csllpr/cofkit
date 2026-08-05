# Technical Overview

## Design notes

The long-term direction is:

1. reactant-aware input models
2. reaction grammar / template library
3. discrete assignment over motifs and reaction events
4. periodic product graph construction
5. initial embedding followed by a lightweight continuous optimization pass over cell scale, rigid poses, and bridge geometry
6. ranked candidate ensembles instead of a single structure

`COFProject.stacking_mode` remains `"disabled"`; open-ended stacking search is still out of scope. Named ring-forming bilayers can be requested through `COFProject.stacking_ids`, while binary-bridge stacking remains in the batch/single-pair layer. Both routes reuse the same opt-in post-build registry enumerator for eligible `2D` outputs.

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
7. CIF export into `cifs/valid`, `cifs/warning`, `cifs/needs_optimization`, or `cifs/hard_invalid`
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

The optimizer is intentionally modest. Bridge-forming candidates use lateral cell scaling, monomer translation updates, and lightweight orientation cleanup. Ring-forming candidates use a separate virtual-node geometry profile and shared-precursor translation refinement, scored and validated against ring radius, 120-degree incidence, and planarity. Stacking scoring remains out of scope; named bilayer registries duplicate the validated layer graph and preserve ring-center offsets during atomistic realization.

For the current imine realization path, two geometry details are now important enough to treat as part of the documented behavior:

- template-specific imine motif-origin correction is applied in the supported `3D` builder paths as well as the earlier `2D` paths, so high-connectivity `dia`-style builds do not silently bypass the bent-linkage span correction
- periodic-image bridge events store realized atom overrides back in the base monomer-local frame before CIF export, which avoids pathological retained-hydrogen directions on image-crossing imine events

The CIF exporter is deliberately honest as well: if a `MonomerSpec` carries atom coordinates, it writes atomistic sites; if not, it falls back to a legal coarse CIF built from monomer centers and motif origins so the current assembly can still be inspected.

## CIF Decomposition

`cofkit analyze decompose` is the reverse path from an atomistic CIF back to COFid. The current implementation prefers explicit CIF bond loops, infers missing bond orders from local geometry when bond labels are present without type fields, and falls back to periodic distance-based bond detection when bond loops are absent. Binary-bridge strategies cut marked linkage bonds and repair the resulting monomer classes. Ring strategies recognize complete boroxine or triazine rings, remove their cyclic bonds, regenerate the boronic-acid or nitrile endpoints on the recovered precursor fragments, and serialize the result with either a caller-supplied topology or a conservatively auto-detected topology.

The nitrogen binary-linkage markers share one bond classifier so their assignments are mutually exclusive. For each eligible C=N bond, an N-N environment resolves through `azine > hydrazone > imine`, while the independent keto-aldehyde-derived environment resolves through `bken > imine`. Classification is local to the bond rather than promoted across the CIF; a bond that satisfies both independent branches is withheld from all four strategies and reported as cross-branch ambiguous in decomposition metadata.

The vinylene marker first collects formal C=C bonds, preserves the existing same-instance rejection when COFKit labels are available, and then rejects bonds belonging to local five- or six-membered rings. Remaining bonds require carbon anchors on both endpoints and a unique activated side whose anchor is electron-withdrawing through multiple nitrogen neighbors or a multiple bond to N/O/S. If the same molecular graph contains recognized boronate-ester B-O bonds or a complete boroxine ring, that boron chemistry takes priority and the structurally eligible vinylene bonds are suppressed. This prevents Kekule-encoded aromatic bonds, ordinary cyclic or unactivated alkenes, and C=C motifs embedded in boron-linked structures from being cut as vinylene linkages while retaining the rejected candidates in diagnostic metadata.

Triazine ring decomposition is evaluated after the binary imine and vinylene strategies. If either strategy independently recovers monomers, it takes priority and the triazine candidate is suppressed. The vinylene check includes the boron-linkage override above, so a C=C candidate already rejected in favor of boron chemistry cannot suppress triazine. Full monomer recovery is required for the imine branch because the alternating formal C=N bonds within a Kekule triazine ring are not, by themselves, evidence of an imine linkage.

COFid build input and reverse CIF decomposition support the binary-bridge linkage codes `imine`, `hydrazone`, `azine`, `boest`, `bken`, and `vinylene`, plus the one-precursor ring-forming codes `boroxine` and `triazine`. Template-id aliases resolve to their canonical COFid linkage codes.

Topology auto-detection first trusts embedded `# COFid:` comments, then ranks cofkit repository topologies against the recovered periodic linkage graph. Ring events become virtual three-connected topology nodes; two-connected precursor fragments are suppressed into net edges, while higher-connected precursors remain explicit nodes. Equivalent disconnected layer graphs from named stacked bilayers are normalized to one layer before ranking. Detection can operate without explicit CIF bond loops and without cofkit-style atom labels, but it is not a general net recognizer: it remains tied to the supported decomposition chemistry and to topologies available in cofkit's local topology repository. Ambiguous candidates are reported instead of guessed, and decorated special cases such as `bex` may still need an explicit `--topology`.

The CIF extraction layer is local to `cofkit` and uses `gemmi`; it does not add ASE as a runtime dependency. The binary-linkage decomposition approach was adapted from the deCOFpose project at <https://github.com/r-fedorov/deCOFpose>; cofkit adds the ring-event and virtual-node reverse strategy used for boroxine and triazine.

## COFid Validation

`cofkit validate simple` wraps decomposition as a direct CIF-vs-COFid check. It forces distance-inferred bond detection, even when explicit CIF bond rows exist, then compares recovered monomer blocks and linkage against the supplied COFid. The supplied topology is only used as an explicit decomposition input and is intentionally excluded from the validation comparison.

`cofkit validate optimize` runs the default LAMMPS optimization workflow first, then applies the same forced distance-inferred decomposition to the optimized CIF. This keeps the validation path focused on whether the post-optimization geometry still decomposes to the expected monomer/linkage identity.
