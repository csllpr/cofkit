# Event-Based Decomposition

COFKit provides two CIF decomposition engines:

- `event` is the default detector → event → reconstruction-hypothesis → global-validation engine.
- `legacy` is the retained per-family compatibility engine.

Event mode became the default after outperforming legacy on the largest available inspected linkage-label collection, the partially labelled CoRE-COFs 1242 table. This set was also used while improving the detector, so it is an operational benchmark rather than an independent holdout. Legacy mode remains available for compatibility and comparison.

## Usage

CLI:

```bash
cofkit analyze decompose framework.cif \
  --linkage auto \
  --json
```

Python:

```python
from cofkit import decompose_cif_to_cofid

event_result = decompose_cif_to_cofid("framework.cif")
legacy_result = decompose_cif_to_cofid(
    "framework.cif",
    decomposition_mode="legacy",
)
```

Both engines accept the same `topology`, `linkage`, and `bond_mode` arguments. Omitting `decomposition_mode` uses `event`; pass `legacy` explicitly for the compatibility implementation.

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

Atomic events currently include paired `C=N-N=C` azines, their keto-tautomerized `C-N-N-C` representation, acylhydrazones, paired five-membered `B-O-C-C-O` boronate esters, complete boroxine rings, and two alternating nitrile reconstructions for each triazine ring. N-N context takes precedence over an overlapping beta-ketoenamine interpretation. Vinylene activation ranks candidate orientations; tied orientations branch, and structurally valid zero-score candidates remain low confidence instead of being discarded immediately. Activated-methylene validation includes conjugated aza-aromatic and cyano-aromatic donors while continuing to reject an unactivated methylbenzene. Canonical beta-ketoenamine events locally suppress only an overlapping vinylene event.

Triazine is handled globally. When another complete reconstruction leaves every structural triazine ring intact inside a recovered monomer, the ring is classified as a triazine-containing monomer motif rather than a triazine linkage. This policy considers every supported competing family, not only imine and vinylene.

## Result diagnostics

The normal `CifDecompositionResult` shape is preserved. Event-specific information is stored under `result.metadata`:

- `decomposition_mode`: always `event`;
- `event_status`: detailed final classification;
- `event_detection`: accepted and locally suppressed events;
- `hypothesis_generation`: site counts, bounded-enumeration diagnostics, and potential non-overlapping family combinations;
- `hypotheses`: every evaluated hypothesis, its events, repaired roles, validation status, and failure reasons;
- `benchmark_contract`: records the default mode, retained legacy availability, and benchmark basis.

Detailed event statuses include `SUCCESS_COMPLETE`, `AMBIGUOUS_MULTIPLE_DECOMPOSITIONS`, `FAILED_CHEMICAL_VALIDATION`, `FAILED_ENDPOINT_ACCOUNTING`, `FAILED_TOPOLOGY_VALIDATION`, `FAILED_UNEXPLAINED_FRAMEWORK`, `UNSUPPORTED_LINKAGE`, and `SUPPRESSED_TRIAZINE_MOTIF`.

## Current limitations

- COFid currently serializes one linkage family. Event mode reports non-overlapping family combinations, but it does not serialize a mixed-linkage COFid.
- Hypothesis enumeration chooses one interpretation per detected site and is capped at 256 combinations per family. It does not yet search arbitrary subsets of high-confidence sites.
- Guest handling is conservative: disconnected components with no event atoms are ignored, while unexplained fragments from the event-bearing framework component invalidate a hypothesis.
- The same `P1`, bond-source, topology-repository, and supported-linkage restrictions as legacy mode still apply.
- The inspection table labels 920 of 1,242 structures and was consulted during detector improvement; reported scores are therefore not independent estimates of generalization.

## Benchmark basis

Both engines were run on identical CIFs, topology inputs, linkage requests, and bond modes from `CoRE-COFs_1242-v7.0/COF-linkage-inspection.csv`. Of 920 labelled structures, 797 have one supported linkage-family label, 6 are explicit multi-linkage labels, 4 are alternative labels, and 113 are outside the eight supported families.

On the 797 supported single-family labels, event mode achieved 65.2% exact accuracy, 69.0% unique coverage, 94.5% selective accuracy when unique, and 76.7% macro F1. Legacy achieved 46.8%, 56.2%, 83.3%, and 58.8%, respectively. Across all 920 labels, strict decision accuracy was 68.2% for event and 52.0% for legacy. In paired scoring, event gained 151 cases that legacy missed and lost 4 that legacy got right. Two alternating-order 16-worker timing runs averaged 5.851 seconds for event and 8.469 seconds for legacy over the 920 CIFs.

Future evaluations should continue to compare at least:

- exact linkage-label accuracy;
- exact reconstructed precursor identity and connectivity;
- complete-success, ambiguous, partial, unsupported, and failure rates;
- false-positive rate by linkage family;
- triazine motif-versus-linkage decisions;
- mixed-linkage cases;
- runtime and hypothesis truncation;
- per-family confusion matrices and failure-stage distributions.

Raw decomposition count alone is not sufficient: linkage correctness, unsupported-family false positives, ambiguity, and precursor reconstruction remain part of the release criterion.
