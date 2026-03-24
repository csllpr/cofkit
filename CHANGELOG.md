# Changelog

All notable changes to `cofkit` are recorded here.

## 2026-03-24

### Added
- An installable unified `cofkit` CLI in `cofkit.cli`, including `single-pair`, `batch-binary-bridge`, `batch-all-binary-bridges`, `classify-output`, `build-default-library`, and template-listing entry points.
- `cofkit.batch_models` for neutral batch-facing summaries with `reactant_a_*` / `reactant_b_*` fields plus compatibility aliases for older imine-centric consumers.
- `cofkit.monomer_library` for extracted monomer-role resolution, detector-backed generic-library loading, and explicit-library regrouping outside the batch executor itself.
- CLI regression coverage, including direct `single-pair` CIF writing.

### Changed
- The example scripts under `examples/` now act as thin wrappers over the shared installable CLI instead of carrying separate command-line logic.
- Direct single-pair CLI generation now supports motif autodetection and default topology enumeration when explicit motif kinds or topology ids are omitted.
- The current full-batch CLI path has been verified against the wrapper path and reproduces the same per-template results for the practical binary-bridge workflows.

### Verification
- Local test suite passes (`116/116` tests at verification time).
- `cofkit single-pair` successfully wrote CIF output for TAPB + TFB on the autodetected imine path.
- The completed CLI full batch in `out/cli_full_batch_launch_complete_20260323` reproduced the same summary results as the wrapper run in `out/wrapper_full_batch_launch_complete_20260323`.

## 2026-03-23

### Added
- Registry-backed motif metadata in `cofkit.chem.motif_registry`, plus registry-backed RDKit construction and fallback-detection plumbing for future chemistry expansion.
- Practical binary-bridge chemistry support beyond imines for `hydrazone_bridge`, `keto_enamine_bridge`, `boronate_ester_bridge`, and `vinylene_bridge`.
- Atomistic reaction realization for those currently implemented binary-bridge products in CIF export.
- `examples/run_experimental_examples.py` covering literature-style smoke cases for COF-300, TpPa-1, COF-42, COF-5, and V-COF-1.
- Detector-backed library autodiscovery in the batch path, so generic `.txt` SMILES libraries can be grouped by detected role/connectivity without relying only on filename prefixes.
- `examples/build_default_monomers_library.py` plus the generated `examples/default_monomers_library/` snapshot, including `registry.jsonl`, `failures.jsonl`, and detector-regrouped `*_count_N.txt` libraries.
- `examples/run_available_binary_bridge_batches.py` for running one full batch per currently available binary-bridge template from a shared monomer-library directory.
- Regression coverage for the new chemistry registry path, reaction realization, RDKit motif detection, and autodetected library loading.

### Changed
- The generic binary-bridge API is now the practical path for the currently implemented non-imine chemistries, not just an imine compatibility shim.
- `examples/run_binary_bridge_generation.py` now accepts `--auto-detect-libraries` and documents generic-library input.
- The shipped example library now records role-detection mismatches explicitly, including hydrazides that were previously sitting in the amine fixture set and aldehydes that are better classified as keto aldehydes.
- Batch pair generation now uses a real process pool in the practical batch workflows, with a default `8`-worker budget.
- Batch CIF export now classifies written files directly into `cifs/valid`, `cifs/warning`, and `cifs/invalid`, while `hard_hard_invalid` structures are recorded in manifests but blocked from CIF export.
- Coarse validation now has a fourth class, `hard_hard_invalid`, triggered by extreme bridge distances (`>= 2.5 A` actual distance) that are too broken to export.

### Fixed
- `keto_enamine_bridge` export now removes the ortho-hydroxyl hydrogen, keeps the tautomerized oxygen, and shortens the realized local `C=O` bond so `TpPa-1`-style products no longer leave phenolic `O-H` atoms behind.
- `hydrazone_bridge` export now removes both terminal hydrazide hydrogens per linkage instead of only one, fixing the hydrogen count in `COF-42`-style products.

### Verification
- Local test suite passes (`113/113` tests at verification time).
- Experimental chemistry smoke examples now generate CIFs for COF-300 (`dia`), TpPa-1 (`hcb`), COF-42 (`hcb`), COF-5 (`hcb`), and V-COF-1 (`hcb`).
- Generated default monomer library contains `489` registered monomers and `23` failed or ambiguous autodetections.
- Fresh process-pool reruns completed for the two corrected linkages:
  - `hydrazone_bridge`: `176` pairs, `852` successful structures, `752` CIFs written
  - `keto_enamine_bridge`: `3406` pairs, `19662` successful structures, `13260` CIFs written

## 2026-03-22

### Added
- `cofkit.indexed_topology_layouts`, a builder-ready layout layer for chemistry-compatible indexed topologies beyond the handcrafted one-node families.
- Shared indexed-topology generation support in batch and direct single-pair flows, including explicit indexed `node-linker` and `node-node` construction through the normal topology-builder registry.
- `cofkit.validation` plus `examples/classify_batch_output.py` for coarse post-generation classification of exported structures.
- Quantitative triage classes `valid`, `warning`, and `hard_invalid`, with categorized output trees and manifest annotations for manual inspection.
- Regression coverage for indexed-topology generation and validation triage.

### Changed
- Default topology selection now includes a curated indexed-topology subset when the bundled chemistry metadata and current builder support agree, restoring practical defaults such as `sql`, `kgm`, `hxl`, `pts`, `ctn`, `bor`, `kgd`, `tbo`, `dia`, `pcu`, `acs`, `lon`, and `qtz`.
- Batch and direct single-pair generation no longer stop at the explicit one-node families; they can now route into the indexed-topology layout path for supported chemistry/topology combinations.
- Coarse validation now separates moderate bridge-geometry drift (`warning`) from clearly broken structures (`hard_invalid`) instead of using a single binary invalid bucket.

### Verification
- Local test suite passes (`88/88` tests at verification time).
- Full default-selector batch output was reclassified successfully into `74,472` `valid`, `2,562` `warning`, and `68,678` `hard_invalid` structures.

## 2026-03-19

### Added
- Chemistry-facing `two_monomer_*` topology metadata, recording allowed `node-node` and `node-linker` binary combinations such as `3+3`, `4+6`, and `3+2`.
- A generic `2D` symmetry-expansion path in `cofkit.topology_symmetry`, so asymmetric-unit `2D` topologies can be analyzed the same way as the broader `3D` catalog.
- Regression coverage for chemistry-mode topology metadata on representative bundled `2D` and `3D` topologies.

### Changed
- `zero_linker_*` metadata remains available, but it is now documented and treated as a lower-level periodic-graph diagnostic rather than the primary chemistry compatibility signal.
- Bundled topology browsing now exposes precomputed `two_monomer_*` metadata for bulk inspection without reparsing full topology definitions.
- The bundled topology report now focuses on chemistry-facing binary pairing modes instead of only raw bipartite status.

### Fixed
- Mixed-connectivity periodic nets such as `cpa` are no longer misread as chemically valid two-monomer node-node topologies just because the infinite periodic graph is bipartite in the gain-graph sense.

### Verification
- Local test suite passes (`79/79` tests at verification time).
- Bundled topology index regenerated successfully with `two_monomer_compatible = True` for `1844` topologies, `False` for `738`, and `None` for `15`.
- Remaining scan-unavailable residue is now `11` bundled `2D` entries and `4` bundled `3D` entries.

## 2026-03-18

### Added
- `cofkit.topology_analysis` for precomputed topology metadata, including raw `zero_linker_*` periodic-graph diagnostics, scan provenance, and builder-support flags.
- `cofkit.topology_symmetry`, an optional `gemmi`-backed space-group expansion layer that reconstructs indexed topologies into explicit `P1` node/edge graphs for metadata scans.
- Workspace index rebuild helpers in `cofkit.topology_data` plus the public convenience entry point `rebuild_imported_topology_index()`.
- Regression coverage for precomputed workspace-index browsing and generic `3D` symmetry-expansion scans on real bundled topologies.

### Changed
- Bundled `src/cofkit/data/topologies/index.json` now precomputes topology-scan metadata so bulk topology browsing does not need per-topology hydration.
- The zero-linker scan now versions its metadata so improved scan logic invalidates stale index rows automatically.
- Generic `3D` topology scanning now covers most of the bundled `3D` catalog through symmetry expansion instead of falling back to `unavailable`.
- Packaging now exposes an optional `crystal` extra for `gemmi`.

### Verification
- Local test suite passes (`76/76` tests at verification time).
- Bundled topology index regenerated successfully with scan-method counts: `gemmi-space-group-expanded-p1-scan 2397`, `expanded-p1-periodic-bipartite-scan 10`, `quotient-graph-parity-scan 14`, `unavailable 176`.
- Remaining `3D` `unavailable` entries dropped to `4`, with the residual uncovered set now dominated by older unsupported `2D` groups rather than the previously dominant `3D` families.

## 2026-03-17

### Added
- Symmetry-expanded one-node `2D` topology support in the batch pipeline for `hcb`, `hca`, `fes`, and `fxt`.
- Explicit batch builders for expanded `3+2` node-linker and bipartite expanded `3+3` node-node cases.
- One-node `3D` topology support for `dia` and `pcu`, including explicit primitive-cell builders for `4+4`, `4+2`, and `6+2` imine generation.
- Direct single-pair support for the shared expanded one-node `2D` / `3D` builders through `COFEngine` and `BatchStructureGenerator.generate_monomer_pair_candidate(s)`.
- Regression coverage for expanded-topology batch generation, topology selection, and direct single-pair engine support for the same one-node families.
- A linkage-profile registry in `cofkit.reactions` for binary-bridge pair-role ordering, target bridge distances, batch-library prefixes, and event-realizer dispatch.
- A shared `cofkit.topology_builders` registry so supported one-node topology construction is selected once and reused by batch and direct single-pair flows.
- Generic binary-bridge batch entry points plus a compatibility example script, `examples/run_binary_bridge_generation.py`, that reproduces the legacy imine batch outputs through the new registry-based API.

### Changed
- Batch runs now enumerate all applicable default one-node topologies per pair, spanning the supported `2D` and `3D` families.
- The batch CLI now writes CIFs by default and uses `--no-write-cif` / `--no-all-topologies` to opt out.
- The legacy two-node `hcb` node-node case now uses a self-contained pair builder instead of routing back through `COFEngine`, so the shared batch/single-pair topology path stays recursion-free.
- Documentation now describes the shared batch and single-pair topology path instead of a batch-only limitation.
- Default high-connectivity routing now prefers `3D` one-node nets (`dia` for `4`-connected cases and `pcu` for `6+2`) instead of forcing the uploaded nonplanar fixtures through the older `2D` families.
- Batch and direct single-pair generation now resolve supported binary-bridge pair order, topology-family dispatch, and target bridge distances through shared registries instead of hard-wired imine-specific branches.

### Fixed
- `fes` and `fxt` batch outputs no longer collapse onto the old `hcb` two-node placeholder embedding.
- `hca` is now treated as a topology-distinct expanded one-node net instead of being forced through an invalid `3+3` alternation.
- Supported expanded one-node topologies now generate directly from the normal single-pair engine path instead of failing or requiring the batch workflow.
- The new generic binary-bridge compatibility example reproduces the same manifest rows and CIF exports as the legacy imine batch example when run with `--template-id imine_bridge`.

### Verification
- Local test suite passes (`72/72` tests at verification time).
- Full all-topologies fixture batch over `examples/batch_test_monomers` completed with `23,963/23,963` successful pairs, `77,601/77,601` successful structures, and `77,601` CIFs written.
- Full-run topology counts: `hcb 18,500`, `fes 18,500`, `fxt 18,500`, `hca 16,638`, `dia 3,147`, `pcu 2,316`.

## 2026-03-15

### Added
- Batch imine-generation support for monomer libraries, including manifest/summary output and optional CIF export.
- Regression coverage for asymmetric `3+2` and `3+3` `hcb` placements plus a direct geometry rotation-mapping check.

### Changed
- `3+2` node-linker `hcb` embedding now uses generalized oblique 2D cells derived from the actual edge vectors instead of forcing a symmetric hexagonal metric.
- `3+3` single-node-bipartite `hcb` embedding now uses the same lower-symmetry cell construction when asymmetric trigonal nodes cannot satisfy a symmetric hexagonal metric.
- Batch linker placement now uses corrected frame-axis rotations and true reactive-site separations rather than projected linker spans.

### Fixed
- Frame-to-axis rotations now respect the matrix convention used by rigid-body placement.
- Asymmetric `3+2` and `3+3` fixture cases now close all bridge distances exactly instead of leaving disconnected reacting sites.

### Verification
- Local test suite passes (`40/40` tests at verification time).
- Full fixture batch over `examples/batch_test_monomers` completed with `18,500/18,500` successful pairs and `18,500` CIFs written.
- Manifest checks show `0` `3+3` pairs and `0` `3+2` pairs above `0.05 A` bridge-distance error.

## 2026-03-14

### Added
- Explicit CIF bond-loop export for atomistic candidates via `_geom_bond_*` records.
- Bond storage on `MonomerSpec` so atomistic monomer connectivity can survive through export.
- Tracking of removed atom ids during reaction realization so CIF export can omit bonds to removed leaving-group atoms.
- Regression coverage for atomistic bond-loop export and reacted RDKit imine CIF output.
- Project action checklist for tracking the next engineering steps.
- Regression coverage for a twisted imine-like bridge pose so the optimizer/scorer preserves a normal-misalignment residual signal.

### Changed
- RDKit-backed monomer construction now preserves atom-bond connectivity.
- Geometry-based motif detection now carries any available molecule bond information into `MonomerSpec`.
- CIF export now combines retained intramonomer bonds with realized inter-monomer reaction bonds instead of only exporting the realized bridge bonds.
- Example report text now describes the broader bond-loop export more accurately.
- Bridge-geometry scoring now includes a weighted normal-misalignment residual as a rigid-body torsion-style surrogate for planar or torsion-restricted bridge templates such as imines.
- Bridge score metadata now reports the normal-misalignment residual per bridge event, and optimizer acceptance now sees that residual through the aggregate bridge residual.
- The lightweight orientation cleanup now points each bridge participant toward its partner consistently during rigid-body refinement.

### Fixed
- Duplicate CIF bond records are filtered before export.
- Bond export now avoids referencing atoms removed during reaction realization.

### Verification
- Local test suite passes (`30/30` tests at verification time).
- Exported reacted artifact `out/tapb_tfb_hcb_rdkit_imine_product.cif` contains `_geom_bond_atom_site_label_1`, `_geom_bond_atom_site_label_2`, and `_geom_bond_distance` fields.
