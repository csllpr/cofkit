# Adding Linkages And Monomer Types

This document is the practical checklist for extending `cofkit` chemistry.

Use it when adding:

- a new motif kind such as `hydrazide`, `catechol`, or a future `anhydride`
- a new linkage such as `hydrazone_bridge` or a future non-imine binary bridge

For general navigation, read [docs/AGENT_CODEBASE_MAP.md](../docs/AGENT_CODEBASE_MAP.md) first.

## 1. Decide what kind of chemistry it is

Before touching code, classify the new linkage:

- `binary_bridge`
  - two monomers, one linkage event between them
  - fits the current practical pipeline best
- not `binary_bridge`
  - ring-forming, hyperedge-like, or multi-step cyclization chemistry
  - do not force it through the current binary-bridge path unless the chemistry is genuinely reducible to that form

Current practical end-to-end support is strongest for `binary_bridge` templates.

## 2. If you are adding a new motif kind

### Required files

1. [src/cofkit/chem/motif_registry.py](../src/cofkit/chem/motif_registry.py)
2. [src/cofkit/chem/rdkit.py](../src/cofkit/chem/rdkit.py)

### Optional / situational files

3. [src/cofkit/chem/detector.py](../src/cofkit/chem/detector.py)
4. [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py)

### Steps

1. Register the motif kind in [src/cofkit/chem/motif_registry.py](../src/cofkit/chem/motif_registry.py)
   - choose `kind`
   - choose `id_prefix`
   - set a reasonable `cif_symbol`
   - list `allowed_reaction_templates`
   - add `rdkit_smarts` if RDKit autodetection should support it
2. Add an RDKit match handler in [src/cofkit/chem/rdkit.py](../src/cofkit/chem/rdkit.py)
   - the handler must return a `_DetectedMotif`
   - make sure `reactive_atom_id`, `anchor_atom_id`, and `atom_ids` are chemically meaningful
3. If needed, add a postprocess handler in [src/cofkit/chem/rdkit.py](../src/cofkit/chem/rdkit.py)
   - use this for site grouping or deduplication
   - current `keto_aldehyde` handling is the main example
4. Decide whether the lightweight fallback detector in [src/cofkit/chem/detector.py](../src/cofkit/chem/detector.py) should also support the motif
   - if not, document that RDKit is required for practical ingestion
5. If autodetect ambiguity needs chemistry-specific suppression, update [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py)
   - current example: generic `aldehyde` is suppressed when `keto_aldehyde` is also found

## 3. If you are adding a new linkage

### Required files

1. [src/cofkit/reactions.py](../src/cofkit/reactions.py)
2. [src/cofkit/reaction_realization.py](../src/cofkit/reaction_realization.py)

### Usually required files

3. [src/cofkit/chem/motif_registry.py](../src/cofkit/chem/motif_registry.py)
4. [src/cofkit/chem/rdkit.py](../src/cofkit/chem/rdkit.py)
5. [src/cofkit/batch.py](../src/cofkit/batch.py)
6. [src/cofkit/cli.py](../src/cofkit/cli.py)

### Steps

1. Add the `ReactionTemplate` in [src/cofkit/reactions.py](../src/cofkit/reactions.py)
   - set `id`
   - set `reactant_motif_kinds`
   - set `product_name`
   - set planarity / torsion priors honestly
2. Add the `ReactionLinkageProfile` in the same file
   - set `bridge_target_distance`
   - set the binary roles in chemically correct order
   - set `event_realizer`
   - choose `workflow_family`
   - choose `topology_assignment_mode`
   - choose `library_layout`
3. Register atomistic realization in [src/cofkit/reaction_realization.py](../src/cofkit/reaction_realization.py)
   - add a registry entry in `ReactionEventRealizationRegistry.builtin()`
   - implement the event handler
   - remove the right atoms
   - keep or reposition atoms if tautomerization/cyclization requires it
   - add explicit inter-monomer bonds
4. Check whether the linkage geometry can reuse the current coarse model
   - current bridge target distances propagate automatically from [src/cofkit/reactions.py](../src/cofkit/reactions.py) into embedding, scoring, and validation
   - if the linkage needs more than a target distance, you may need extra work in [src/cofkit/embedding.py](../src/cofkit/embedding.py), [src/cofkit/scoring.py](../src/cofkit/scoring.py), or [src/cofkit/validation.py](../src/cofkit/validation.py)
5. Make sure batch input can discover the right monomer libraries
   - explicit-library loading uses the linkage profile’s role `library_prefix`
   - autodetected libraries use motif detection plus [src/cofkit/monomer_library.py](../src/cofkit/monomer_library.py)
6. If the new linkage should be part of practical workflows, update user-facing docs and example scripts

## 4. Current constraints to respect

These are common sources of bad patches:

- The practical generation path still expects exactly one binary-bridge template per run.
- Batch autodetection still requires exactly one resolved motif kind per monomer.
- Non-binary / ring-forming chemistries are not yet fully integrated into topology-guided batch generation.
- Validation thresholds are still mostly global; do not assume they are linkage-specific.
- Topology compatibility and builder support are separate concepts.

## 5. Topology-related checklist for a new linkage

If the linkage is still a normal two-monomer topology-guided case:

1. Confirm the monomer connectivities map onto existing `two_monomer_*` topology metadata.
2. Confirm default topology selection in [src/cofkit/batch.py](../src/cofkit/batch.py) actually includes the intended topology family.
3. If not, inspect:
   - [src/cofkit/topology_analysis.py](../src/cofkit/topology_analysis.py)
   - [src/cofkit/topology_builders.py](../src/cofkit/topology_builders.py)
   - [src/cofkit/indexed_topology_layouts.py](../src/cofkit/indexed_topology_layouts.py)

If the linkage is not a clean two-monomer topology-guided case, do not fake it with topology metadata just to make it run.

## 6. Testing checklist

For a new motif kind:

- [tests/test_rdkit_monomer.py](../tests/test_rdkit_monomer.py)
- optionally [tests/test_batch.py](../tests/test_batch.py) if autodetected libraries should pick it up

For a new linkage:

- [tests/test_core.py](../tests/test_core.py)
- [tests/test_reaction_realization.py](../tests/test_reaction_realization.py)
- [tests/test_batch.py](../tests/test_batch.py)
- [tests/test_cli.py](../tests/test_cli.py)

If you add a literature-style smoke example:

- [examples/run_experimental_examples.py](../examples/run_experimental_examples.py)

## 7. What not to do

- Do not add a new linkage only in [src/cofkit/reactions.py](../src/cofkit/reactions.py) and call it done.
- Do not route chemically distinct products through the imine realizer.
- Do not silently coerce mixed-role monomers into one role in autodetect mode.
- Do not broaden topology defaults without checking current builder support.
- Do not use CIF export success alone as proof the chemistry is implemented correctly.

## 8. Fast sanity questions before you finish

- Can the new motif kind be detected from SMILES?
- Does the linkage appear in `cofkit list-templates`?
- Can `cofkit single-pair` build and write a CIF for one literature-style pair?
- Does batch loading work from both explicit and autodetected libraries, if both are supposed to be supported?
- Does the CIF atom/bond stoichiometry match the intended chemistry?
