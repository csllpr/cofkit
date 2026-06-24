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

`COFEngine` / `COFProject.stacking_mode` should remain `"disabled"`. Current stacking support is the post-build registry enumeration path exposed by `BatchStructureGenerator` and CLI `--stacking`.

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
