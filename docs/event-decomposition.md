# Event-Based Decomposition

COFKit provides two CIF decomposition engines:

- `event` is the default detector → event → reconstruction-hypothesis → global-validation engine.
- `legacy` is the retained per-family compatibility engine.

Event mode uses globally validated reconstruction hypotheses. Legacy mode remains available for compatibility and comparison. The default selection does not establish accuracy on a new dataset; validate linkage and precursor recovery for the structures you use.

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
4. enumerates at most 256 hypotheses per linkage family, in a canonical order derived from graph invariants (RDKit canonical atom symmetry classes) rather than input atom/label order, so permuting the CIF atom rows cannot change which hypotheses the cap explores;
5. cuts and repairs each event set atomically;
6. checks endpoint ownership and unexplained framework fragments;
7. reuses the established precursor-motif, periodic-rank, topology, COFid parse, and forward-build validators;
8. selects a result only from complete globally valid hypotheses.

When `bond_mode="auto"` uses an explicit CIF bond table but fails during endpoint, framework-accounting, or chemical validation, event mode may also evaluate an independently distance-inferred graph. The same fallback is allowed for a failed or rank-zero linkage-graph extraction only when at least one event-bearing source component independently has periodic rank 2 or 3. It is not used for ambiguous topology ranking or for a source framework whose own covalent graph has insufficient periodic rank. The fallback also remains disabled when several recovered same-role species are already individually buildable. The alternative replaces the explicit graph only if it produces exactly one globally validated decomposition; otherwise the original result and diagnostics are retained. Reconstructed fragments may also be merged across bond-order/H-placement variants, but only when molecular formula, element-labelled constitutional graph, role, connectivity, and net charge agree and at least one variant exposes the required number of buildable motifs. Constitutional isomers are never merged by this normalization, and a merge is never reported as exact identity: the result's `precursor_identity` block distinguishes `exact` recovery from `element_graph_ambiguous` canonicalization and preserves every input form with its fragment count.

Atomic events currently include paired `C=N-N=C` azines, their keto-tautomerized `C-N-N-C` representation, acylhydrazones, paired five-membered `B-O-C-C-O` boronate esters, complete boroxine rings, and two alternating nitrile reconstructions for each triazine ring. N-N context takes precedence over an overlapping beta-ketoenamine interpretation. Vinylene activation ranks candidate orientations; tied orientations branch, and structurally valid zero-score candidates remain low confidence instead of being discarded immediately. Activated-methylene validation includes conjugated aza-aromatic and cyano-aromatic donors plus methyl-bearing phenyl rings directly conjugated to an aromatic heterocycle, while continuing to reject unactivated methylbenzene and biaryl donors. Primary-amine precursor validation excludes resonance-deactivated amide, thioamide, sulfonamide, and thiourea nitrogens, while accepting terminal hydrazino sites only in a guanidinium-like three-nitrogen carbon environment. Beta-ketoenamine precursor/build validation supports both ortho-hydroxy aldehyde tautomerization and the explicit `CHO-CH2-C(=O)-C` beta-ketoenol route. Canonical beta-ketoenamine events locally suppress only an overlapping vinylene event.

Triazine is handled globally. When another complete reconstruction leaves every structural triazine ring intact inside a recovered monomer, the ring is classified as a triazine-containing monomer motif rather than a triazine linkage. This policy considers every supported competing family, not only imine and vinylene.

## Result diagnostics

The normal `CifDecompositionResult` shape is preserved. Event-specific information is stored under `result.metadata`:

- `decomposition_mode`: always `event`;
- `event_status`: detailed final classification;
- `search_status` / `search_coverage`: whether the bounded hypothesis search was `complete` or `truncated` at the per-family cap, with explored-versus-theoretical hypothesis counts and the truncated families; when a failed search was truncated, the result reason states explicitly that absence of a valid decomposition is not established;
- `hypothesis_status_counts`: how many considered hypotheses ended in each event status, so rejected alternatives keep their counts alongside the per-hypothesis reasons in `hypotheses`;
- `successful_alternatives`: every distinct complete decomposition (family, topology, COFid, score, hypothesis id) — a one-entry list on success, the full ambiguity set on `AMBIGUOUS_MULTIPLE_DECOMPOSITIONS`;
- `supporting_hypothesis_count`: on success, how many complete hypotheses converged on the selected COFid before deduplication;
- `precursor_identity`: the identity basis of the recovered precursors — `exact` when every recovered fragment of a reaction role had one identical canonical SMILES, or `element_graph_ambiguous` when chemically distinct bond-order/tautomer forms shared one molecular formula and element graph and were collapsed onto a deterministic buildable representative. The ambiguous case preserves every input form with its fragment count under `alternatives`; the selected form is a canonicalization, not an exact chemical identity;
- `atom_ledger`: the full atom accounting of the selected (or best failed) reconstruction — see the next section;
- `event_detection`: accepted and locally suppressed events;
- `hypothesis_generation`: site counts, bounded-enumeration diagnostics, and potential non-overlapping family combinations;
- `hypotheses`: every evaluated hypothesis, its events, repaired roles, validation status, and failure reasons;
- `bond_graph_fallback`: when attempted, the explicit-graph and distance-graph outcomes and whether the fully validated fallback was selected;
- `source_framework_periodicity`: periodic rank, connected-component sizes, and bond counts for the event-bearing source covalent graph, evaluated per component rather than by combining unrelated components;
- `fragment_identity_normalization`: formula/element-graph equivalent fragment forms and the selected buildable representative;
- `multispecies_precursor_recovery`: when present, every distinct recovered species, its declared and detected motif count, fragment count, and role-wise reactive-site balance;
- `defect_detection`: when present, the dominant precursor combination, fragment-agreement fraction, minority fragments, and structural-glitch evidence for a probable defective input;
- `benchmark_contract`: compatibility metadata recording the default mode, legacy availability, and selection unit; it does not certify predictive accuracy.

## Atom ledger and charge/stereo scope

Every event-mode cut/reconstruction produces a typed `AtomLedger` (`cofkit.decompose_ledger`), serialized in full as `metadata["atom_ledger"]` on the selected or best-failed result and as a compact `atom_ledger_summary` on each serialized hypothesis. The ledger assigns every input CIF atom to exactly one bucket:

- `recovered_fragments`: product atoms mapped to a fragment repaired into a precursor (per-fragment atom indices, element counts, role, recovered SMILES);
- `guests`: disconnected components with no linkage-event atoms. Each guest molecule is one ledger entry with its own atom count, Hill formula, and net formal charge; `guest_molecule_count` is deliberately distinct from `guest_atom_count`;
- `residue`: product-side atoms excluded from precursor recovery by the linkage chemistry (currently boroxine ring oxygens), with a recorded `residue_policy`;
- `unaccounted`: framework fragments not assignable to a precursor role. These still invalidate the hypothesis with `FAILED_UNEXPLAINED_FRAMEWORK`; the ledger itemizes their atoms rather than only counting them.

Because precursor repair restores condensation byproduct atoms and re-expresses product-state hydrogens, raw fragment-vs-input atom sums are not a valid conservation check. The ledger therefore records `reaction_additions` (atoms restored during repair, e.g. one O per imine aldehyde endpoint; two hydroxyl O per boronic-acid endpoint against three residue ring oxygens per boroxine ring — net 3 H2O per ring) and `reaction_deletions` (product-state explicit hydrogens re-expressed as implicit hydrogens), and verifies the explicit-atom balance per element:

    input + reaction_additions = recovered_precursor_explicit_atoms + reaction_deletions + residue + guests + unaccounted

Recovered-precursor implicit-hydrogen and total-atom counts are reported informationally alongside.

Charge and stereochemistry scope is conservative:

- `charge_scope`: formal charges from a CIF `_atom_site_pdbx_formal_charge` column are applied to the decomposition graph and verified per fragment against the recovered precursor SMILES (formal charge survives canonicalization and rides into the COFid). Statuses: `preserved` (charged input fragments recovered with the same net charge), `guest_localized` (charge only on reported guest molecules, which the COFid does not capture), `restored_by_heuristic` (charge created during repair by the valence-driven tetravalent-nitrogen rule on an input that declared none), `limited` (charge altered/lost on a recovered fragment, or charge on residue/unaccounted atoms — an explicit unsupported marker, never silent loss), or `not_present`. `charge_conserved` checks net charge across precursors + guests + residue + unaccounted.
- `stereo_scope`: decomposition never preserves stereochemistry — canonical SMILES and the COFid format use `isomericSmiles=False`. Input stereo evidence (tetrahedral atom chirality perceived from the 3D coordinates; bond E/Z perception is unavailable in this RDKit line) is recorded under `input_identity_evidence`, and `stereo_scope.status` is `unsupported` whenever evidence was detected, `not_detected` otherwise, `undetermined` if perception could not run. `preserved` is always `false`.

Atom-resolved formal-charge export into atomistic CIF charge columns, and ionic or stereochemical assembly, remain out of scope. The legacy decomposition engine does not compute an atom ledger.

Detailed event statuses include `SUCCESS_COMPLETE`, `AMBIGUOUS_MULTIPLE_DECOMPOSITIONS`, `DETECTED_PROBABLE_STRUCTURAL_DEFECT`, `UNSUPPORTED_MULTISPECIES_PRECURSORS`, `UNSUPPORTED_MIXED_LINKAGE_FAMILY`, `FAILED_CHEMICAL_VALIDATION`, `FAILED_ENDPOINT_ACCOUNTING`, `FAILED_TOPOLOGY_VALIDATION`, `FAILED_UNEXPLAINED_FRAMEWORK`, `FAILED_INTERNAL_ERROR`, `UNSUPPORTED_LINKAGE`, and `SUPPRESSED_TRIAZINE_MOTIF`. `FAILED_INTERNAL_ERROR` marks an unexpected internal failure (the original `TypeName: message` cause is preserved in `validation_errors` and an `internal_error` metadata block) and surfaces as result `status="error"`; it is never a chemical verdict. `UNSUPPORTED_MIXED_LINKAGE_FAMILY` marks a framework whose non-overlapping linkage families each reconstruct part of the structure but cannot be serialized under the one-family COFid contract — an explicit unsupported-representation verdict, not a chemical incompatibility verdict.

## Probable structural defects

For failed binary-linkage hypotheses, event mode can distinguish a mostly regular framework containing a small number of damaged fragments from an ordinary incomplete decomposition. Every reaction role must have one unique dominant monomer, the dominant combination must account for strictly more than 75% of recovered fragment instances, and the minority reconstruction must contain independent structural glitches such as role-connectivity loss, a reactive-site deficit, or unexplained framework fragments. Multiple same-role species alone are insufficient, so legitimate multivariate structures are not classified as defective merely because they contain more than two precursor species.

A probable defect remains a `skipped` decomposition and never produces a guessed COFid. JSON output includes the machine-readable `defect_detection` report and `probable_linkage`; normal CLI output prints the candidate linkage, agreement fraction, dominant monomers, and glitch summary before exiting unsuccessfully. Conversely, a role-complete set of multiple distinct species is reported as `UNSUPPORTED_MULTISPECIES_PRECURSORS` only when every species independently exposes its declared motif count and the role-wise reactive-site totals balance. This status is deliberately agnostic about whether the source represents an intentional multivariate material, disorder, or another data provenance issue.

## Current limitations

- COFid currently serializes one linkage family. When detected non-overlapping families each genuinely reconstruct part of the framework, event mode aborts with the explicit `UNSUPPORTED_MIXED_LINKAGE_FAMILY` verdict; it does not serialize a mixed-linkage COFid.
- Multivariate COFs with several chemically distinct precursors in the same reaction role are not yet supported as complete decompositions. Such reconstructions abstain with `UNSUPPORTED_MULTISPECIES_PRECURSORS`; external provenance is still required before calling a particular CIF intentionally multivariate. Distinct precursor identities are retained; they are not collapsed into a binary COFid or classified as defects solely because multiple species are present.
- Hypothesis enumeration chooses one interpretation per detected site and is capped at 256 combinations per family. It does not yet search arbitrary subsets of high-confidence sites. When the cap cuts off combinations, the result carries `search_status: "truncated"` with explored/theoretical coverage counts instead of presenting the search as exhaustive. Enumeration order is canonical with respect to input atom/label permutation, so the explored subset — and therefore the verdict under truncation — is independent of CIF row order. Assessment note: raising the cap from 256 to 4096 changed no verdict across the repository's decomposition fixtures, generated imine/vinylene round trips, and a 25-structure CoRE-COF sample (largest theoretical hypothesis count observed: 72), so the default cap is unchanged.
- Guest handling is conservative: disconnected components with no event atoms are treated as guests — now itemized in the atom ledger (per-molecule atom counts, formula, charge) and included in the atom balance — while unexplained fragments from the event-bearing framework component invalidate a hypothesis and are itemized as `unaccounted`.
- The same `P1`, bond-source, topology-repository, and supported-linkage restrictions as legacy mode still apply.
