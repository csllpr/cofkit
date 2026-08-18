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
    # Summaries are returned best-first by mean bridge-geometry residual.
    # summary.score is None unless legacy scoring is enabled.
    print(summary.topology_id, summary.cif_path)
```

Candidates and summaries carry `metadata["scoring_mode"]` (`"residual"` by default). The deprecated event-count heuristic score can be restored with `BatchGenerationConfig(enable_legacy_scoring=True)` or `--legacy-scoring` on the CLI; it then populates `summary.score` / `candidate.score` and drives ranking again.

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

Use `decompose_cif_to_cofid` for supported binary-bridge and boroxine/triazine ring-forming atomistic CIFs expanded in `P1`. The default `linkage="auto"` evaluates every supported family and succeeds only when one produces a complete periodic precursor decomposition. Pass an explicit linkage to test one family. Pass `topology=` when you know the topology, or omit it to use conservative topology detection from a validated embedded COFid comment or the recovered periodic linkage graph. The default `bond_mode="auto"` prefers explicit CIF bonds and falls back to distance inference; use `bond_mode="distance"` to force distance-inferred connectivity. Distance inference is triclinic-safe, and explicit atom pairs may retain more than one periodic edge when their image gains differ.

The event/hypothesis engine is the default. It normalizes the graph once, detects atomic local linkage events, branches alternative site interpretations, and selects only globally validated reconstruction hypotheses. Pass `decomposition_mode="legacy"` to use the retained per-family compatibility engine. Both return the same `CifDecompositionResult` type; the default event result includes event and hypothesis diagnostics under `metadata`. See [event-decomposition.md](event-decomposition.md) for its status model, limitations, and benchmark results.

In default event mode, nitrogen-containing C-N/C=N overlaps are resolved locally before hypotheses are evaluated globally. N-N context recognizes azine and acylhydrazone products in both their conventional and keto-tautomerized representations before an overlapping beta-ketoenamine interpretation. The beta-ketoenamine detector accepts the conventional product tautomer with a C-N linkage bond and adjacent C=C bond. Event decisions and suppressed overlaps are exposed under `metadata["event_detection"]`; the explicit legacy engine retains its per-bond classifier diagnostics under `metadata["nitrogen_linkage_detection"]`.

Vinylene matching is also conservative: formal C=C bonds in five- or six-membered rings are excluded, both alkene endpoints must have carbon anchors, and candidate orientations are ranked by the activated-methylene-derived environment. Conjugated aza-aromatic and cyano-aromatic donors are supported while unactivated methylbenzene is rejected. Tied orientations branch into separate hypotheses for global validation.

Triazine is resolved globally. If another complete supported-family reconstruction retains every detected triazine ring inside a recovered monomer, event mode classifies that ring as a monomer motif instead of a linkage. This avoids treating the alternating formal C=N bonds intrinsic to a Kekule triazine ring as competing chemistry. Legacy-only resolution details remain available under `metadata["triazine_linkage_resolution"]` when that engine is requested.

```python
from cofkit import decompose_cif_to_cofid, detect_cif_topology

detection = detect_cif_topology("out/tapb_tfb.cif", linkage="imine")
print(detection.status, detection.selected_topology, detection.confidence)

result = decompose_cif_to_cofid("out/tapb_tfb.cif", linkage="auto", bond_mode="auto")
if result.ok:
    print(result.cofid)
else:
    print(result.status, result.reason)

legacy_result = decompose_cif_to_cofid(
    "out/tapb_tfb.cif",
    linkage="auto",
    decomposition_mode="legacy",
)
print(result.metadata["event_status"])
```

`status="skipped"` denotes unsupported, incomplete, chemically inconsistent, or nonperiodic recovery; `status="ambiguous"` means automatic linkage detection found more than one valid family. `status="error"` denotes malformed CIF, periodicity, bond, or topology input. A successful result has the complete linkage-specific precursor role set, exact motif-derived connectivity, a periodic rank compatible with its topology, and has passed the forward build-input validator.

For ring-forming structures, decomposition recognizes complete B3O3 or C3N3 product rings, restores the precursor reactive groups, and rebuilds the periodic net through virtual three-connected ring nodes. Equivalent components in named stacked bilayers are reduced to one layer graph before topology ranking.

Topology detection only ranks topologies available in cofkit's local topology repository and remains scoped to the currently supported decomposition chemistries. Ambiguous cases return diagnostics instead of guessing; pass `topology="bex"` or another explicit topology when needed.
