# Experimental Event-Based Decomposition

COFKit provides two CIF decomposition engines:

- `legacy` is the established per-family engine and remains the default.
- `event` is an experimental detector → event → reconstruction-hypothesis → global-validation engine.

The event engine is additional functionality. It does not replace or silently change the legacy result. No superiority claim is made until both engines have been evaluated against an independent labelled test set.

## Usage

CLI:

```bash
cofkit analyze decompose framework.cif \
  --decomposition-mode event \
  --linkage auto \
  --json
```

Python:

```python
from cofkit import decompose_cif_to_cofid

legacy_result = decompose_cif_to_cofid("framework.cif")
event_result = decompose_cif_to_cofid(
    "framework.cif",
    decomposition_mode="event",
)
```

The event engine accepts the same `topology`, `linkage`, and `bond_mode` arguments. Callers must opt in explicitly; omitting `decomposition_mode` continues to use `legacy`.

## Pipeline

Event mode normalizes the CIF graph once and then:

1. detects structured, immutable `LinkageEvent` objects;
2. applies only local chemical precedence;
3. groups alternative interpretations of the same site;
4. enumerates at most 256 hypotheses per linkage family;
5. cuts and repairs each event set atomically;
6. checks endpoint ownership and unexplained framework fragments;
7. reuses the established precursor-motif, periodic-rank, topology, COFid parse, and forward-build validators;
8. selects a result only from complete globally valid hypotheses.

Atomic events currently include paired `C=N-N=C` azines, paired five-membered `B-O-C-C-O` boronate esters, complete boroxine rings, and two alternating nitrile reconstructions for each triazine ring. Vinylene activation ranks candidate orientations; tied orientations branch, and structurally valid zero-score candidates remain low confidence instead of being discarded immediately. Canonical beta-ketoenamine events locally suppress only an overlapping vinylene event.

Triazine is handled globally. When another complete reconstruction leaves every structural triazine ring intact inside a recovered monomer, the ring is classified as a triazine-containing monomer motif rather than a triazine linkage. This policy considers every supported competing family, not only imine and vinylene.

## Result diagnostics

The normal `CifDecompositionResult` shape is preserved. Event-specific information is stored under `result.metadata`:

- `decomposition_mode`: always `event`;
- `event_status`: detailed final classification;
- `event_detection`: accepted and locally suppressed events;
- `hypothesis_generation`: site counts, bounded-enumeration diagnostics, and potential non-overlapping family combinations;
- `hypotheses`: every evaluated hypothesis, its events, repaired roles, validation status, and failure reasons;
- `benchmark_contract`: records that legacy remains the default and the independent benchmark is pending.

Detailed event statuses include `SUCCESS_COMPLETE`, `AMBIGUOUS_MULTIPLE_DECOMPOSITIONS`, `FAILED_CHEMICAL_VALIDATION`, `FAILED_ENDPOINT_ACCOUNTING`, `FAILED_TOPOLOGY_VALIDATION`, `FAILED_UNEXPLAINED_FRAMEWORK`, `UNSUPPORTED_LINKAGE`, and `SUPPRESSED_TRIAZINE_MOTIF`.

## Current limitations

- COFid currently serializes one linkage family. Event mode reports non-overlapping family combinations, but it does not serialize a mixed-linkage COFid.
- Hypothesis enumeration chooses one interpretation per detected site and is capped at 256 combinations per family. It does not yet search arbitrary subsets of high-confidence sites.
- Guest handling is conservative: disconnected components with no event atoms are ignored, while unexplained fragments from the event-bearing framework component invalidate a hypothesis.
- The same `P1`, bond-source, topology-repository, and supported-linkage restrictions as legacy mode still apply.
- The external labelled benchmark set is not available yet, so legacy remains the production default.

## Benchmark contract

When the labelled set is ready, run both engines on identical CIFs, topology inputs, linkage requests, and bond modes. Compare at least:

- exact linkage-label accuracy;
- exact reconstructed precursor identity and connectivity;
- complete-success, ambiguous, partial, unsupported, and failure rates;
- false-positive rate by linkage family;
- triazine motif-versus-linkage decisions;
- mixed-linkage cases;
- runtime and hypothesis truncation;
- per-family confusion matrices and failure-stage distributions.

Changing the default should require event mode to improve labelled chemical correctness without an unacceptable regression in coverage or runtime. Raw decomposition count alone is not sufficient.
