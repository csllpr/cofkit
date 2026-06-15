# Agent Codebase Map

This document is for later coding/agent sessions that need to move quickly in `cofkit` without rebuilding the mental model from scratch.

## First orientation

Read these files first, in this order:

1. [README.md](../README.md)
2. [USER_MANUAL.md](../USER_MANUAL.md)
3. [ARCHITECTURE.md](../ARCHITECTURE.md)
4. [CHANGELOG.md](../CHANGELOG.md)

Then go straight to the module that matches the task.

## Main user entry points

- [src/cofkit/cli.py](../src/cofkit/cli.py)
  - Installed `cofkit` CLI.
  - Routes the top-level `build`, `analyze`, `calculate`, and `validate` namespaces.
- [src/cofkit/cli_build.py](../src/cofkit/cli_build.py)
  - Owns build-facing commands such as `cofkit build single-pair` and the batch workflows.
- [src/cofkit/cli_analyze.py](../src/cofkit/cli_analyze.py)
  - Owns analysis-facing commands such as `cofkit analyze classify-output`, `cofkit analyze decompose`, and `cofkit analyze zeopp`.
- [src/cofkit/cli_calculate.py](../src/cofkit/cli_calculate.py)
  - Owns external calculation commands such as `cofkit calculate lammps-optimize`, `graspa-widom`, `graspa-isotherm`, `graspa-mixture`, and `hybrid-mdmc`.
- [src/cofkit/engine.py](../src/cofkit/engine.py)
  - Direct project-style API via `COFEngine`.
- [src/cofkit/batch.py](../src/cofkit/batch.py)
  - Practical execution layer for topology-guided single-pair and batch generation.

## Core chemistry seams

- [src/cofkit/reactions.py](../src/cofkit/reactions.py)
  - Reaction templates and linkage profiles.
  - This is the first file to touch for a new linkage.
- [src/cofkit/chem/motif_registry.py](../src/cofkit/chem/motif_registry.py)
  - Motif-kind metadata.
  - Add new motif kinds here first.
- [src/cofkit/chem/rdkit.py](../src/cofkit/chem/rdkit.py)
  - Practical SMILES-to-`MonomerSpec` path.
  - New motif kinds normally need a match handler here.
- [src/cofkit/chem/detector.py](../src/cofkit/chem/detector.py)
  - Lightweight non-RDKit fallback detector.
  - Only some motif kinds are implemented here.
- [src/cofkit/reaction_realization.py](../src/cofkit/reaction_realization.py)
  - Atomistic bond/deletion realization for CIF export.
  - A new linkage is not really complete without this.

## Monomer-library and batch-input seams

- [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py)
  - `MonomerRoleResolver` for autodetection from SMILES.
  - `BinaryBridgeLibraryLoader` for explicit and autodetected batch libraries.
- [src/cofkit/batch_models.py](../src/cofkit/batch_models.py)
  - Neutral batch-facing data classes.
  - Use these instead of adding new ad hoc summary dictionaries.

## Topology seams

- [src/cofkit/topology_builders.py](../src/cofkit/topology_builders.py)
  - Shared dispatch for supported topology-family builders.
- [src/cofkit/single_node_topologies.py](../src/cofkit/single_node_topologies.py)
  - Handcrafted / expanded `2D` one-node families.
- [src/cofkit/single_node_topologies_3d.py](../src/cofkit/single_node_topologies_3d.py)
  - Supported `3D` one-node families.
- [src/cofkit/indexed_topology_layouts.py](../src/cofkit/indexed_topology_layouts.py)
  - Generic indexed-topology layout reconstruction.
- [src/cofkit/topology_analysis.py](../src/cofkit/topology_analysis.py)
  - Chemistry-facing `two_monomer_*` metadata and lower-level graph diagnostics.
- [src/cofkit/topology_symmetry.py](../src/cofkit/topology_symmetry.py)
  - Generic symmetry expansion for topology analysis.

## Geometry / scoring / validation

- [src/cofkit/embedding.py](../src/cofkit/embedding.py)
  - Initial periodic placement.
- [src/cofkit/optimizer.py](../src/cofkit/optimizer.py)
  - Lightweight continuous refinement.
- [src/cofkit/scoring.py](../src/cofkit/scoring.py)
  - Candidate scoring and bridge-geometry metrics.
- [src/cofkit/validation.py](../src/cofkit/validation.py)
  - `valid` / `warning` / `hard_invalid` / `hard_hard_invalid` triage.
- [src/cofkit/cif.py](../src/cofkit/cif.py)
  - CIF export, including realized inter-monomer bonds.
- [src/cofkit/decompose_cif.py](../src/cofkit/decompose_cif.py)
  - CIF atom-site and explicit-bond extraction for decomposition without adding ASE.
- [src/cofkit/decompose.py](../src/cofkit/decompose.py)
  - Explicit-bond binary-bridge decomposition into recovered monomers and COFid serialization.
  - This logic was adapted from the deCOFpose project at `https://github.com/r-fedorov/deCOFpose`.

## External calculation seams

- [src/cofkit/lammps.py](../src/cofkit/lammps.py)
  - UFF/DREIDING-backed explicit-bond CIF to LAMMPS data/input generation, EQeq charge staging, local minimization/MD orchestration, optional guest-restart atom/force-field merging for MD, and optimized or MD-updated CIF export.
- [src/cofkit/graspa.py](../src/cofkit/graspa.py)
  - EQeq to gRASPA/RASPA2 Widom, single-component isotherm, and mixture workflows; framework mixing-rule generation; simulation.input rendering; result parsing.
- [src/cofkit/hybrid_mdmc.py](../src/cofkit/hybrid_mdmc.py)
  - Cyclic LAMMPS MD plus gRASPA/RASPA2 GCMC workflow. Default framework exchange passes the MD-updated framework CIF between segments; opt-in guest-restart exchange parses final GCMC guest snapshots and passes massive guest atoms into the following LAMMPS MD segment.
- [src/cofkit/guest_restart.py](../src/cofkit/guest_restart.py)
  - GCMC movie/restart snapshot discovery, guest force-field asset synchronization from packaged/bundle RASPA rows into LAMMPS-ready guest sites/templates, binary guest parsing, and zero-mass pseudo-site rejection.
- [src/cofkit/guest_bundles.py](../src/cofkit/guest_bundles.py)
  - Shared external guest parameter-bundle contract for gRASPA/RASPA2 workflows. Bundles add RASPA molecule definitions, pseudo-atom rows, mixing-rule rows, aliases, rotatability, and a required non-empty `lammps` section for synchronized hybrid MD/MC guest restart.

## Typical execution paths

### Single pair from CLI

1. `cofkit build single-pair` in [src/cofkit/cli_build.py](../src/cofkit/cli_build.py)
2. motif-kind autodetection via [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py) if explicit kinds are omitted
3. monomer construction via [src/cofkit/chem/rdkit.py](../src/cofkit/chem/rdkit.py)
4. candidate generation via [src/cofkit/batch.py](../src/cofkit/batch.py)
5. CIF export and validation bucketing via [src/cofkit/cif.py](../src/cofkit/cif.py) and [src/cofkit/validation.py](../src/cofkit/validation.py)

### Full batch from CLI

1. `cofkit build batch-binary-bridge` or `cofkit build batch-all-binary-bridges`
2. library resolution via [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py)
3. pair enumeration and topology selection in [src/cofkit/batch.py](../src/cofkit/batch.py)
4. topology-family dispatch via [src/cofkit/topology_builders.py](../src/cofkit/topology_builders.py)
5. validation-aware CIF writing into `cifs/valid`, `cifs/warning`, or `cifs/invalid`

### CIF decomposition from CLI

1. `cofkit analyze decompose <CIF> --topology <TOKEN>` in [src/cofkit/cli_analyze.py](../src/cofkit/cli_analyze.py)
2. atom and bond extraction via [src/cofkit/decompose_cif.py](../src/cofkit/decompose_cif.py)
3. linkage-specific cutting and monomer repair via [src/cofkit/decompose.py](../src/cofkit/decompose.py)
4. COFid serialization through [src/cofkit/cofid.py](../src/cofkit/cofid.py)

### gRASPA/RASPA2 calculation from CLI

1. `cofkit calculate graspa-widom`, `graspa-isotherm`, or `graspa-mixture` in [src/cofkit/cli_calculate.py](../src/cofkit/cli_calculate.py)
2. optional external guest loading and alias canonicalization via [src/cofkit/guest_bundles.py](../src/cofkit/guest_bundles.py) and [src/cofkit/graspa.py](../src/cofkit/graspa.py)
3. EQeq framework charge staging in [src/cofkit/graspa.py](../src/cofkit/graspa.py)
4. backend-specific `simulation.input` and force-field asset materialization in [src/cofkit/graspa.py](../src/cofkit/graspa.py)
5. parsed JSON/CSV report writing from backend `Output/**/*.data`

### Hybrid MD/MC calculation from CLI

1. `cofkit calculate hybrid-mdmc` in [src/cofkit/cli_calculate.py](../src/cofkit/cli_calculate.py)
2. cycle orchestration in [src/cofkit/hybrid_mdmc.py](../src/cofkit/hybrid_mdmc.py)
3. per-cycle LAMMPS MD framework update through [src/cofkit/lammps.py](../src/cofkit/lammps.py)
4. per-cycle pure-component or mixture GCMC through [src/cofkit/graspa.py](../src/cofkit/graspa.py)
5. in `framework` mode, the next cycle starts from the LAMMPS MD output CIF only
6. in `guest_restart` mode, [src/cofkit/guest_restart.py](../src/cofkit/guest_restart.py) parses the final GCMC `Movies/System_0/result_*.data` snapshot and [src/cofkit/lammps.py](../src/cofkit/lammps.py) injects those guests into the next MD data file

## Current architectural constraints

These are important before editing:

- The practical generation path is still centered on one binary-bridge template per run.
- Batch autodetection still requires one resolved motif kind per monomer.
- RDKit is the practical chemistry path; the fallback detector is narrower.
- Non-binary / ring-forming chemistries are not yet fully wired into topology-guided generation.
- Validation thresholds are still globally configured, not fully per-linkage.

## Where to add tests

- [tests/test_cli.py](../tests/test_cli.py)
  - CLI coverage.
- [tests/test_batch.py](../tests/test_batch.py)
  - Batch and single-pair generation behavior.
- [tests/test_rdkit_monomer.py](../tests/test_rdkit_monomer.py)
  - RDKit motif detection and monomer building.
- [tests/test_reaction_realization.py](../tests/test_reaction_realization.py)
  - Atomistic reaction realization.
- [tests/test_decompose_cif.py](../tests/test_decompose_cif.py)
  - No-ASE CIF extraction and explicit-bond preservation.
- [tests/test_decompose.py](../tests/test_decompose.py)
  - CIF-to-COFid decomposition and generated hcb round trips across buildable binary-bridge linkages.
- [tests/test_lammps.py](../tests/test_lammps.py)
  - LAMMPS data/input generation, force-field parameter paths, optimization/MD orchestration, and CIF preservation behavior.
- [tests/test_graspa.py](../tests/test_graspa.py)
  - EQeq/gRASPA/RASPA2 workflow staging, CLI parsing, guest bundles, force-field asset generation, and parser behavior.
- [tests/test_hybrid_mdmc.py](../tests/test_hybrid_mdmc.py)
  - Hybrid MD/MC cycle orchestration, framework snapshot handoff behavior, and binary guest-restart propagation.
- [tests/test_guest_restart.py](../tests/test_guest_restart.py)
  - Guest restart force-field synchronization, GCMC movie snapshot parsing, and unsupported massless pseudo-site handling.
- [tests/test_core.py](../tests/test_core.py)
  - End-to-end project / template behavior.

## Practical advice for later sessions

- Start from the registry file, not the executor.
- If a new chemistry exports CIFs incorrectly, inspect [src/cofkit/reaction_realization.py](../src/cofkit/reaction_realization.py) before changing scoring or topology code.
- If autodetected libraries behave strangely, inspect [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py) before [src/cofkit/batch.py](../src/cofkit/batch.py).
- If a topology is present in metadata but never generated, inspect:
  - [src/cofkit/topology_analysis.py](../src/cofkit/topology_analysis.py)
  - [src/cofkit/topology_builders.py](../src/cofkit/topology_builders.py)
  - the default-selector logic in [src/cofkit/batch.py](../src/cofkit/batch.py)
