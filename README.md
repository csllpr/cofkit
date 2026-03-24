# cofkit

A reaction-aware periodic assembly toolkit for covalent organic frameworks (COFs).

## Overview

`cofkit` is a foundation for a reaction-aware COF toolkit. It brings together chemistry-aware monomer definitions, reaction templates, topology selection, periodic assembly, scoring, validation, and export in one package.

The current practical workflows focus on topology-guided binary-bridge generation, but the architecture is intended to grow toward broader chemistry support, richer topology coverage, and a wider structure-generation toolkit over time.

Core abstractions:

- **MonomerSpec** — a reactant monomer with annotated reactive motifs
- **ReactiveMotif** — a chemically typed reaction site with a local frame
- **ReactionTemplate** — a reaction grammar entry (bridge-forming or ring-forming)
- **ReactionEvent** — an instantiated reaction in a candidate periodic assembly
- **PeriodicProductGraph** — the framework connectivity after reactions are applied
- **COFEngine** — the orchestration layer for reaction selection, discrete assignment, and seed assembly

Current user-facing functionality centers on:

- direct project-style generation through `COFEngine`
- single-pair screening and library-scale runs through `BatchStructureGenerator`
- installed CLI workflows for single-pair generation, batch generation, classification, and example-library building

Usage guidance for practical structure generation is collected in [USER_MANUAL.md](USER_MANUAL.md).

Developer / agent navigation guides:

- [docs/AGENT_CODEBASE_MAP.md](docs/AGENT_CODEBASE_MAP.md)
- [docs/ADDING_LINKAGES_AND_MONOMERS.md](docs/ADDING_LINKAGES_AND_MONOMERS.md)
- [skills/cofkit-navigator/SKILL.md](skills/cofkit-navigator/SKILL.md)
  - Codex skill for routing natural-language end-user requests into the current `cofkit` CLI and Python workflows.

## Current scope

The current package still keeps mandatory runtime dependencies at zero so the scientific core can be exercised in a clean Python runtime. Optional extras add RDKit-backed monomer construction and gemmi-backed CIF/topology utilities for the practical workflows.

Implemented so far:

- geometry primitives and local frames
- core domain dataclasses
- builtin COF reaction library plus linkage-profile metadata for pair-role ordering, bridge targets, and realization hooks
- periodic product graph with motif usage validation
- topology repository/index support for RCSR CGD-style bundles, with bundled topology data preferred by default and builtin fallback hints
- a discrete assignment layer for motif-to-reaction event matching
- registry-backed motif metadata plus lightweight geometric motif detection for the currently supported fallback kinds
- optional RDKit/SMARTS-backed monomer construction for amine, aldehyde, hydrazide, boronic acid, catechol, keto aldehyde, and activated-methylene motifs, including conformer generation and bond retention when RDKit is installed
- initial linkage geometry helpers for bridge-forming reactions
- initial periodic embedding for monomer instances from topology hints and motif/reaction heuristics, including oblique `hcb` cells for asymmetric `3+3` and `3+2` cases when a symmetric hexagonal metric is too restrictive
- a dependency-free continuous optimization pass for seed cell/pose refinement after embedding
- first-pass candidate scoring with event coverage, bridge geometry, topology bonuses, and unreacted penalties
- candidate metadata carrying embedding provenance, optimizer metrics, and score breakdowns
- legal P1 CIF export, with atomistic output when monomer coordinates are available and a coarse fallback otherwise
- registry-backed atomistic reaction realization for the currently implemented binary-bridge products: imine, hydrazone, beta-ketoenamine, boronate ester, and vinylene
- reaction-aware batch binary-bridge generation over monomer libraries, with imine workflows as the primary documented path, including `3+3`, `3+2`, `4+4`, `4+2`, and `6+2` enumeration, manifest/summary writing, and CIF export enabled by default
- automatic monomer-role detection for batch library loading, so generic `.txt` SMILES libraries can be regrouped by detected role/connectivity instead of relying only on `*_count_N.txt` filenames
- an auto-generated example library under [`examples/default_monomers_library`](examples/default_monomers_library) with detector-scanned role/count registration metadata
- symmetry-expanded single-node generation across supported `2D` and `3D` one-node families, available through both batch runs and single-pair input, with `3+2` / `3+3` handling on `hcb` / `hca` / `fes` / `fxt`, plus `dia` and `pcu` builders for `4+2`, `4+4`, and `6+2` cases
- default topology-family routing that now sends higher-connectivity nonplanar inputs to `3D` one-node nets (`dia` for `4`-connected cases, `pcu` for `6+2`) while keeping older `2D` high-connectivity families opt-in
- shared topology-builder dispatch for the supported one-node families, reused by batch runs and direct single-pair generation
- indexed-topology layout generation for chemistry-compatible bundled topologies beyond the handcrafted one-node families, available in both batch and explicit single-pair generation
- compatibility-aware default topology selection that now includes curated indexed topologies such as `sql`, `kgm`, `hxl`, `pts`, `ctn`, `bor`, `kgd`, `tbo`, `dia`, `pcu`, `acs`, `lon`, and `qtz` when the current chemistry metadata and builder support permit them
- coarse post-generation validation with `valid` / `warning` / `hard_invalid` / `hard_hard_invalid` triage, categorized CIF output trees, and export blocking for obviously broken structures
- process-level batch pair generation with an `8`-worker default budget for the practical CLI workflows
- an installable `cofkit` CLI, including direct `single-pair` generation plus unified batch/classification/library-building entry points
- extracted batch-facing support layers (`cofkit.monomer_library`, `cofkit.batch_models`) so CLI, wrappers, and future linkage extensions share the same monomer-role resolution and summary schema
- tests covering core invariants

Not implemented yet:

- broader SMARTS/rule-based motif detection beyond the current binary-bridge set, especially for ring-forming and cyclization chemistries
- symmetry reduction beyond raw indexed topology filtering
- full general topology coverage beyond the current supported one-node families plus the present indexed-layout subset
- torsion-aware or force-field-backed optimization beyond the current lightweight pass
- ring-closure geometry models
- chemically faithful atomistic CIF generation for arbitrary monomers without fallback/pseudo-sites
- semiempirical / force-field cleanup
- any stacking exploration, registry search, or stacking score terms

## Installation and dependencies

Core install:

```bash
python3 -m pip install -e .
```

After editable install, the unified CLI is available as:

```bash
cofkit --help
```

Core `cofkit` has no mandatory third-party runtime dependencies beyond Python `3.10+`.

Optional extras you may want in your environment:

- `gemmi` for rebuilding topology metadata with generic space-group expansion, CIF-backed coarse validation, and broader `2D` / `3D` two-monomer compatibility scans
- `RDKit` for SMILES-based monomer construction and the practical batch-generation workflows
- `pytest` to run the test suite

The bundled topology repository under [`src/cofkit/data/topologies`](src/cofkit/data/topologies) is sufficient for normal use. External RCSR archives and topology environment variables are optional advanced inputs, not required setup steps.

## Quick example

```python
from cofkit import (
    COFEngine,
    COFProject,
    Frame,
    MonomerSpec,
    ReactiveMotif,
)

tri_amine = MonomerSpec(
    id="tapb",
    name="TAPB-like triamine",
    motifs=(
        ReactiveMotif(id="n1", kind="amine", atom_ids=(1,), frame=Frame.xy()),
        ReactiveMotif(id="n2", kind="amine", atom_ids=(2,), frame=Frame.yz()),
        ReactiveMotif(id="n3", kind="amine", atom_ids=(3,), frame=Frame.zx()),
    ),
)

tri_aldehyde = MonomerSpec(
    id="tfp",
    name="TFP-like trialdehyde",
    motifs=(
        ReactiveMotif(id="c1", kind="aldehyde", atom_ids=(4,), frame=Frame.xy()),
        ReactiveMotif(id="c2", kind="aldehyde", atom_ids=(5,), frame=Frame.yz()),
        ReactiveMotif(id="c3", kind="aldehyde", atom_ids=(6,), frame=Frame.zx()),
    ),
)

engine = COFEngine()
project = COFProject(
    monomers=(tri_amine, tri_aldehyde),
    allowed_reactions=("imine_bridge",),
    target_dimensionality="2D",
    target_topologies=("hcb",),
)

results = engine.run(project)
best = results.top(1)[0]
print(best.score)
```

`cofkit` ships with bundled RCSR-derived topology data under [`src/cofkit/data/topologies`](src/cofkit/data/topologies). The current bundle contains:

- `2d/` with 194 extracted `.cgd` files
- `3d/` with 2404 extracted `.cgd` files, with 2403 currently indexed and 1 recorded as unsupported by the parser
- `index.json` with precomputed connectivity/count metadata for fast repository filtering
- `import-metadata.json` recording the bundled import summary and any skipped files

Normal engine usage prefers that bundled data automatically, without any topology environment variables.

If you want to point `cofkit` at external RCSR archives explicitly instead of the bundled copy, you can still do that:

```bash
export COFKIT_TOPOLOGY_2D_ARCHIVE=/path/to/topo_2d.zip
export COFKIT_TOPOLOGY_3D_ARCHIVE=/path/to/topo_3d.zip
```

If the 3D bundle is not named `topo_3d.zip`, either set `COFKIT_TOPOLOGY_3D_ARCHIVE` directly or call `discover_rcsr_archives()` / `TopologyRepository.from_rcsr_archives(...)` with the path you want to use.

To write a legal CIF:

```python
from cofkit import write_candidate_cif

write_candidate_cif("out/mock_imine.cif", best, project.monomers)
```

## CLI

After editable install, the main user-facing interface is the `cofkit` command.

Inspect the available commands with:

```bash
cofkit --help
cofkit list-templates
```

The most useful subcommands are:

- `cofkit single-pair`
- `cofkit batch-binary-bridge`
- `cofkit batch-all-binary-bridges`
- `cofkit classify-output`
- `cofkit build-default-library`

### Single pair

Direct single-pair generation with motif autodetection and default topology enumeration:

```bash
cofkit single-pair \
  --template-id imine_bridge \
  --first-smiles 'C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N' \
  --second-smiles 'C1=C(C=C(C=C1C=O)C=O)C=O' \
  --first-id tapb \
  --second-id tfb \
  --output-dir out/cli_single_pair
```

That command autodetects the monomer roles (`amine` / `aldehyde` for this example), evaluates the current default topology pool for the pair, writes `summary.json`, and exports CIFs into validation buckets under `cifs/valid`, `cifs/warning`, or `cifs/invalid`.

If you want to force one topology instead of using the default topology pool:

```bash
cofkit single-pair \
  --template-id imine_bridge \
  --first-smiles 'C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N' \
  --second-smiles 'C1=C(C=C(C=C1C=O)C=O)C=O' \
  --first-id tapb \
  --second-id tfb \
  --topology hcb \
  --output-dir out/cli_single_pair_hcb
```

### One linkage over a library

Run one binary-bridge linkage over a monomer-library directory:

```bash
cofkit batch-binary-bridge \
  --template-id imine_bridge \
  --input-dir examples/default_monomers_library \
  --output-dir out/batch_imine_from_cli \
  --max-workers 8
```

### All available binary-bridge linkages

Run all currently available binary-bridge templates over one library directory:

```bash
cofkit batch-all-binary-bridges \
  --input-dir examples/default_monomers_library \
  --output-dir out/available_binary_bridge_batches \
  --max-workers 8
```

That installable CLI path has been verified against the older wrapper script and reproduces the same per-template summary results for the current `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge` batches.

### Classify a finished output tree

```bash
cofkit classify-output \
  out/full_cif_generation_default_selector_20260320 \
  --output-dir out/full_cif_generation_default_selector_20260320_coarse_validation_triage
```

### Rebuild the detector-scanned example library

```bash
cofkit build-default-library
```

## Wrapper scripts

The example scripts in `examples/` are still available as thin wrappers over the shared CLI.

For the imine batch wrapper:

```bash
python3 examples/run_batch_imine_generation.py \
  --input-dir examples/batch_test_monomers \
  --output-dir out/batch_imine_generation_full_cif \
  --max-cif-exports 20000 \
  --num-conformers 2
```

That workflow writes a `manifest.jsonl`, a `summary.md`, and by default one CIF per successful topology-specific structure. Use `--no-write-cif` to disable export or `--no-all-topologies` to keep only the best topology per pair. The example script is now a thin wrapper over the shared `cofkit` CLI.

The same workflow is also available through the generic binary-bridge interface:

```bash
python3 examples/run_binary_bridge_generation.py \
  --template-id imine_bridge \
  --input-dir examples/batch_test_monomers \
  --output-dir out/binary_bridge_generation_full_cif \
  --max-cif-exports 20000 \
  --num-conformers 2
```

With `--template-id imine_bridge`, that script reproduces the same topology-specific outputs as the imine wrapper while using the template-driven registry path. The installable `cofkit` CLI remains the primary user-facing interface.

The generic script can also infer monomer roles directly from raw generic SMILES libraries when the filenames do not already encode role information:

```bash
python3 examples/run_binary_bridge_generation.py \
  --template-id keto_enamine_bridge \
  --input-dir /path/to/raw_generic_smiles_libraries \
  --output-dir out/binary_bridge_generation_keto \
  --auto-detect-libraries \
  --num-conformers 2
```

Use `--auto-detect-libraries` only on raw generic inputs whose non-empty rows all resolve to one of the selected template's supported motif kinds. Do not combine it with [`examples/default_monomers_library`](examples/default_monomers_library), which is already detector-grouped explicit library output.

If you want a stable explicit library layout instead of detecting at batch runtime, regenerate the bundled example library with:

```bash
python3 examples/build_default_monomers_library.py
```

That script scans [`examples/batch_test_monomers`](examples/batch_test_monomers), writes regrouped role/count libraries under [`examples/default_monomers_library`](examples/default_monomers_library), and records detector provenance in `registry.jsonl` plus failures/ambiguities in `failures.jsonl`.

For literature-style single-pair smoke cases across the currently implemented binary-bridge chemistries, see:

```bash
python3 examples/run_experimental_examples.py
```

To run all currently available binary-bridge batches from the detector-scanned default library with the default `8`-worker pair pool through the wrapper:

```bash
python3 examples/run_available_binary_bridge_batches.py \
  --input-dir examples/default_monomers_library \
  --output-dir out/available_binary_bridge_batches
```

The shipped [`examples/default_monomers_library`](examples/default_monomers_library) snapshot currently exposes `hydrazone_bridge`, `imine_bridge`, and `keto_enamine_bridge` batches. `boronate_ester_bridge` and `vinylene_bridge` remain covered by the direct single-pair smoke examples above.

To classify a finished batch output into `valid`, `warning`, `hard_invalid`, and `hard_hard_invalid` buckets through the wrapper:

```bash
python3 examples/classify_batch_output.py \
  out/full_cif_generation_default_selector_20260320 \
  --output-dir out/full_cif_generation_default_selector_20260320_coarse_validation_triage
```

That workflow writes one classification manifest plus four CIF trees:

- `valid/`
- `warning/`
- `hard_hard_invalid/`
- `hard_invalid/`

The classifier keeps the source CIFs untouched and materializes the categorized views through symlinks by default.

## Design notes

The long-term direction is:

1. reactant-aware input models
2. reaction grammar / template library
3. discrete assignment over motifs and reaction events
4. periodic product graph construction
5. initial embedding followed by a lightweight continuous optimization pass over cell scale, rigid poses, and bridge geometry
6. ranked candidate ensembles instead of a single structure

For the current phase, stacking exploration is intentionally out of scope.
`stacking_mode` must stay `"disabled"`; the engine now rejects other values rather than implying partial support.

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
2. optional RDKit monomer construction / caching
3. pair enumeration from registered binary-bridge role prefixes and monomer connectivities
4. topology-family-aware candidate generation
5. manifest / summary writing
6. process-level pair execution by default in the practical batch CLIs (`8` workers unless overridden)
7. CIF export into `cifs/valid`, `cifs/warning`, or `cifs/invalid`
8. block CIF export for `hard_hard_invalid` structures while still recording them in the manifest
9. optional reclassification of a finished output tree with `examples/classify_batch_output.py`

For convenience, the documented imine workflow still uses `run_imine_batch(...)` and [examples/run_batch_imine_generation.py](examples/run_batch_imine_generation.py), but the shared implementation now lives under the generic binary-bridge path exposed by [examples/run_binary_bridge_generation.py](examples/run_binary_bridge_generation.py).

For the current single-node topology families, `cofkit` expands supported CGD nets into explicit `P1` node/edge orbits before construction. In practice that means:

- `3+2` runs can currently target `hcb`, `hca`, `fes`, and `fxt`
- `3+3` runs can currently target `hcb`, `fes`, and `fxt`
- `4+2` and `4+4` runs default to `dia` in the batch workflow and are available directly through `COFEngine` with `target_dimensionality="3D"` or `target_topologies=("dia",)`
- `6+2` runs default to `pcu` in the batch workflow and are available directly through `COFEngine` with `target_dimensionality="3D"` or `target_topologies=("pcu",)`
- explicit `sql`, `kgm`, `htb`, and `hxl` requests still work in batch mode, but they are no longer the default route for the uploaded `4`- and `6`-connected fixture libraries
- `BatchStructureGenerator.generate_monomer_pair_candidate(s)` now exposes the full supported one-node family directly for single-pair imine generation, and `COFEngine.run(...)` reaches the same builders for explicit one-node topology requests plus supported `3D` single-pair defaults such as `dia` / `pcu`
- chemistry-compatible indexed topologies from the bundled repository can now also be requested explicitly in batch and direct single-pair flows, and the default selector includes a curated subset of those indexed nets when the current chemistry metadata and builder support agree

The optimizer is intentionally modest. It only refines bridge-forming candidates through lateral cell scaling, monomer translation updates, and lightweight orientation cleanup when those steps reduce bridge-geometry residuals. Ring-forming reactions stay coarse, and stacking remains completely out of scope.

The CIF exporter is deliberately honest as well: if a `MonomerSpec` carries atom coordinates, it writes atomistic sites; if not, it falls back to a legal coarse CIF built from monomer centers and motif origins so the current assembly can still be inspected.

## Coarse validation

The optional coarse validator is meant to separate clearly broken outputs from inspectable but strained ones. It uses four classes:

- `valid`: no warning or hard-invalid criteria triggered
- `warning`: no hard-invalid criteria triggered, but at least one soft bridge-geometry criterion triggered
- `hard_invalid`: at least one clearly broken-network, impossible-geometry, or clash criterion triggered
- `hard_hard_invalid`: bridge geometry is so extreme that CIF export is blocked during generation

Current warning thresholds:

- any bridge distance residual `> 0.75 A`
- mean bridge distance residual `> 0.35 A`
- more than `25%` of bridge events with residual `> 0.50 A`

Current hard-invalid thresholds:

- any unreacted motifs
- missing or unparsable CIF
- any bridge distance residual `> 1.00 A`
- mean bridge distance residual `> 0.60 A`
- any bridge `actual_distance / target_distance < 0.70`
- any bridge `actual_distance / target_distance > 1.60`
- disconnected monomer-instance graph reconstructed from CIF bonding
- any nonbonded heavy-atom contact `< 1.05 A`
- `2D` cell area `< 10.0 A^2`
- `3D` cell volume `< 20.0 A^3`

Current hard-hard-invalid threshold:

- any actual bridge distance `>= 2.5 A`

Hydrogen atoms are ignored in the clash check. Warning-level bridge drift still goes through CIF-backed checks, so a candidate can be promoted from `warning` to `hard_invalid` if the exported structure shows a broken network or impossible heavy-atom contact. During generation, `valid` / `warning` / `hard_invalid` CIFs are written into separate subdirectories under `cifs/`, while `hard_hard_invalid` structures stay manifest-only with `cif_export_blocked = true`.

## Acknowledgements

`cofkit` is being built as an independent toolkit within the broader COF-software ecosystem.

PORMAKE has been an important source of inspiration and an occasional code-reading reference, especially around topology-oriented workflows, supporting utilities, and practical project ergonomics. That contribution is acknowledged here with appreciation.

`cofkit` maintains its own reaction-aware scientific core and is intended to grow into a broader toolkit rather than remain a point-for-point alternative to any single precursor project.
