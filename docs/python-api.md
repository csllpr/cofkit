# Python API

The CLI is the preferred human-facing interface. Use the Python API when you need direct programmatic control or want to integrate `cofkit` into another workflow.

## COFEngine

Use `COFEngine` when you already know the target topology and want a normal single-pair `CandidateEnsemble`.

```python
from cofkit import COFEngine, COFProject, build_rdkit_monomer, write_candidate_cif

tapb = build_rdkit_monomer(
    "tapb",
    "TAPB",
    "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N",
    "amine",
)
tfb = build_rdkit_monomer(
    "tfb",
    "TFB",
    "C1=C(C=C(C=C1C=O)C=O)C=O",
    "aldehyde",
)

project = COFProject(
    monomers=(tapb, tfb),
    allowed_reactions=("imine_bridge",),
    target_dimensionality="2D",
    target_topologies=("hcb",),
)

best = COFEngine().run(project).top(1)[0]
write_candidate_cif("out/tapb_tfb_hcb.cif", best, project.monomers)
```

`COFEngine` / `COFProject.stacking_mode` should remain `"disabled"`. For ring-forming projects, use `COFProject.stacking_ids=("AA", "AB")`; binary-bridge stacking remains the post-build registry path exposed by `BatchStructureGenerator` and CLI `--stacking`.

## RingFormingStructureGenerator

Use the dedicated builder for a single precursor whose functional groups cyclotrimerize into topology nodes:

```python
from cofkit import CIFWriter, RingFormingStructureGenerator, build_rdkit_monomer

precursor = build_rdkit_monomer(
    "ctf1_precursor",
    "terephthalonitrile",
    "N#Cc1ccc(C#N)cc1",
    "nitrile",
)
candidate = RingFormingStructureGenerator().generate(
    precursor,
    "triazine_trimerization",
    topology_id="hcb",
)
CIFWriter().write_candidate("out/ctf1.cif", candidate, {precursor.id: precursor})
```

`boroxine_trimerization` accepts boronic-acid motifs. The default hcb path requires a ditopic precursor; trigonal hcb and compatible indexed mixed-connectivity topologies are also supported. Inspect `candidate.metadata["ring_validation"]` before accepting a geometry.

Enumerate explicit bilayers with:

```python
from cofkit import RingFormationConfig, RingFormingStructureGenerator

generator = RingFormingStructureGenerator(
    RingFormationConfig(stacking_ids=("AA", "AB", "slipped"))
)
stacked_candidates = generator.generate_candidates(
    precursor,
    "triazine_trimerization",
    topology_id="hcb",
)
```

## BatchStructureGenerator

Use `BatchStructureGenerator` when you want the same topology-family-aware generation used by the batch pipeline.

```python
from cofkit import BatchGenerationConfig, BatchStructureGenerator, build_rdkit_monomer

tapb = build_rdkit_monomer(
    "tapb",
    "TAPB",
    "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N",
    "amine",
)
tpal = build_rdkit_monomer("tpal", "TPAL", "O=Cc1ccc(C=O)cc1", "aldehyde")

generator = BatchStructureGenerator(
    BatchGenerationConfig(
        write_cif=True,
        enumerate_all_topologies=True,
    )
)

summaries, candidates, attempted = generator.generate_monomer_pair_candidates(
    tapb,
    tpal,
    out_dir="out/single_pair_all",
)

print("attempted:", attempted)
for summary in summaries:
    print(summary.topology_id, summary.score, summary.cif_path)
```

Restrict topology selection with:

```python
BatchGenerationConfig(
    topology_ids=("hcb", "fes"),
    write_cif=True,
)
```

Opt into `2D` stacking variants with:

```python
BatchGenerationConfig(
    stacking_ids=("AA", "AB"),
    write_cif=True,
)
```

## Batch Libraries

```python
from cofkit import BatchGenerationConfig, BatchStructureGenerator

generator = BatchStructureGenerator(
    BatchGenerationConfig(
        rdkit_num_conformers=4,
        enumerate_all_topologies=True,
        write_cif=True,
    )
)

summary = generator.run_binary_bridge_batch(
    "examples/default_monomers_library",
    "out/binary_bridge_generation",
    template_id="imine_bridge",
)

print(summary.attempted_pairs)
print(summary.successful_structures)
print(summary.manifest_path)
```

The older `run_imine_batch(...)` convenience method remains available and delegates to the generic binary-bridge path with `template_id="imine_bridge"`.

## CIF Decomposition

Use `decompose_cif_to_cofid` for supported binary-bridge and boroxine/triazine ring-forming atomistic CIFs. Pass `topology=` when you know the topology, or omit it to use conservative topology detection from an embedded COFid comment or the recovered periodic linkage graph. The default `bond_mode="auto"` prefers explicit CIF bonds and falls back to distance inference; use `bond_mode="distance"` to force distance-inferred connectivity.

Nitrogen-containing C=N overlaps are resolved per bond, not per CIF. N-N environments use `azine > hydrazone > imine`, while keto-aldehyde-derived environments use the independent `bken > imine` branch. A more-specific bond is therefore unavailable to generic imine decomposition. Cross-branch matches are returned as ambiguous and exposed under `metadata["nitrogen_linkage_detection"]`.

Vinylene matching is also conservative: formal C=C bonds in five- or six-membered rings are excluded, both alkene endpoints must have carbon anchors, and exactly one endpoint must resolve as the more strongly activated-methylene-derived side. Recognized boronate-ester or boroxine chemistry in the same structure takes priority over an otherwise valid vinylene match. Detection and override counts are exposed under `metadata["vinylene_linkage_detection"]`.

Triazine matching yields to independently recoverable imine or vinylene chemistry. This uses complete monomer recovery rather than raw C=N detection, avoiding promotion of the alternating formal C=N bonds intrinsic to a Kekule triazine ring. Resolution details are exposed under `metadata["triazine_linkage_resolution"]`.

```python
from cofkit import decompose_cif_to_cofid, detect_cif_topology

detection = detect_cif_topology("out/tapb_tfb.cif", linkage="imine")
print(detection.status, detection.selected_topology, detection.confidence)

result = decompose_cif_to_cofid("out/tapb_tfb.cif", linkage="imine", bond_mode="auto")
if result.ok:
    print(result.cofid)
else:
    print(result.reason)
```

For ring-forming structures, decomposition recognizes complete B3O3 or C3N3 product rings, restores the precursor reactive groups, and rebuilds the periodic net through virtual three-connected ring nodes. Equivalent components in named stacked bilayers are reduced to one layer graph before topology ranking.

Topology detection only ranks topologies available in cofkit's local topology repository and remains scoped to the currently supported decomposition chemistries. Ambiguous cases return diagnostics instead of guessing; pass `topology="bex"` or another explicit topology when needed.
