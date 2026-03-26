# Sulfur-Enabled Imine Conversion Plan

## Status

This plan supersedes the earlier precursor-based benzothiazole plan that assumed a preinstalled `-SH` on the amine monomer.

Current implementation status:

- Stage 1 imine build remains the canonical graph representation.
- strict imine-monomer eligibility checks are implemented
- sulfur-enabled conversion eligibility annotation is implemented
- atomistic CIF realization now supports a requested benzothiazole-like external-sulfur annulation heuristic for eligible imine events
- the build graph and reaction-event graph still remain `imine_bridge`

The new rule is:

- build a normal imine COF first
- then mark that imine-linked framework as eligible for sulfur-enabled post-build conversion
- do not require sulfur-bearing functionality on the starting amine monomer

## Goal

Implement conservative rules for:

- detecting amine monomers that are valid imine-COF monomers
- building the imine COF from those monomers
- flagging the resulting imine COF as eligible for sulfur-enabled linkage conversion
- doing all of this without requiring a pre-existing `-SH` on the amine monomer

## Core Model

Treat this as a strict two-stage workflow.

### Stage 1: imine-COF eligibility and construction

- amine monomer + aldehyde monomer -> imine-linked COF

### Stage 2: post-build sulfur-conversion eligibility

- imine-linked COF + sulfur-enabled conversion conditions -> converted sulfur-containing linkage candidate

Do not collapse these two stages into a single monomer-detection rule.

## Chemistry Rules

### Valid amine monomer candidate

A valid amine monomer candidate must satisfy all of the following:

- contains at least 2 reactive primary amine sites
- each reactive site is an exocyclic neutral primary amine, `-NH2`
- each reactive amine nitrogen is attached directly to an `sp2` carbon on an aromatic or heteroaromatic scaffold
- the monomer can serve as a polyfunctional node for reticulation with a polyaldehyde monomer
- the reactive amine sites are not already consumed in another functionality
- the monomer scaffold is compatible with imine COF construction

Accept as reactive amine sites:

- aniline-type aryl amines
- heteroaryl primary amines directly attached to a conjugated scaffold
- polyarylamines such as triarylamine-like or triazine-centered aryl amines, if the reactive group itself is a neutral exocyclic `-NH2`
- multitopic monomers with `2`, `3`, `4`, or more primary amine sites

Reject as reactive amine sites:

- secondary amines
- tertiary amines
- amides
- imides
- ureas
- carbamates
- sulfonamides
- amidines
- guanidines
- hydrazides
- hydrazines unless the pipeline is explicitly in azine mode
- anilinium or protonated ammonium forms unless the code explicitly normalizes them to the free base
- amines whose nitrogen is already part of an aromatic ring as a pyridine-like `N` without exocyclic `NH2`
- protected amines

### Important non-requirement

Do not require any of the following on the amine monomer:

- ortho `-SH`
- ortho `-OH`
- ortho `-NH2` for ring-forming annulation
- preinstalled sulfur-containing substituents

Those requirements belong to other linkage templates and must not be imported into this workflow.

## Repo Constraints

The current codebase already has a strong Stage 1 path:

1. `imine_bridge` is a working binary-bridge template in [`src/cofkit/reactions.py`](../src/cofkit/reactions.py).
2. single-pair and batch imine building already run through [`src/cofkit/batch.py`](../src/cofkit/batch.py) and [`src/cofkit/engine.py`](../src/cofkit/engine.py).
3. amine and aldehyde RDKit motif building already exist in [`src/cofkit/chem/rdkit.py`](../src/cofkit/chem/rdkit.py).
4. atomistic imine realization already exists in [`src/cofkit/reaction_realization.py`](../src/cofkit/reaction_realization.py).

What the codebase does **not** currently have:

- a post-build conversion-candidate registry
- a candidate annotation layer for external reagent conversions
- an explicit distinction between "built linkage" and "eligible follow-up conversion"

Because of that, the safest plan is:

- reuse the existing imine build path for Stage 1
- add a separate post-build conversion-eligibility layer for Stage 2
- do not introduce a new build-time benzothiazole template in Stage 1

## Architectural Decision

### What should not be done

Do not implement this workflow as:

- a new build-time `benzothiazole_annulation` monomer-pair template
- a detector for prefunctionalized ortho-aminothiophenol monomers
- a fused-ring transform inside the current `binary_bridge` builder
- a replacement of built imine linkages during normal COF generation

That would violate the updated chemistry model and risk changing the current imine pipeline.

### Recommended architecture

Keep the architecture split:

- Stage 1 uses the existing `binary_bridge` workflow with `template_id="imine_bridge"`
- Stage 2 uses a new post-build conversion profile, separate from build-time linkage generation

Recommended new internal layer:

- module: `cofkit.post_build_conversions`
- purpose: annotate candidates that are eligible for external-reagent conversion workflows

This is not a new build workflow family. It is a second-stage annotation layer over built candidates.

### Why this is the stable approach

- the starting monomers are ordinary imine monomers
- the built framework should remain canonically imine-linked at Stage 1
- sulfur is external chemistry, not monomer provenance
- the converted sulfur-containing product should not be claimed atomistically until a separate second-stage transform is defined and validated

## Data Model Changes

### 1. Keep `imine_bridge` as the canonical built linkage

File: [`src/cofkit/reactions.py`](../src/cofkit/reactions.py)

Do not add a new Stage 1 build template for sulfur conversion.

Stage 1 should continue to use:

- `template_id="imine_bridge"`

Reason:

- the framework is first built as an imine COF
- the imine linkage must remain the actual built linkage stored in the graph and output metadata

### 2. Add a post-build conversion profile registry

Recommended new file:

- `src/cofkit/post_build_conversions.py`

Recommended dataclasses:

- `PostBuildConversionProfile`
- `PostBuildConversionRegistry`
- `PostBuildConversionAssessment`

Recommended first profile id:

- `sulfur_enabled_imine_conversion`

Suggested profile metadata:

- `required_built_template_ids=("imine_bridge",)`
- `required_external_conditions=("sulfur_enabled",)`
- `changes_built_graph=False` for the initial milestone
- `conversion_product_family="sulfur_containing_converted_linkage"`

### 3. Add candidate metadata for conversion eligibility

Files:

- [`src/cofkit/engine.py`](../src/cofkit/engine.py)
- [`src/cofkit/batch.py`](../src/cofkit/batch.py)

Recommended metadata shape:

```python
{
    "post_build_conversions": {
        "sulfur_enabled_imine_conversion": {
            "requested": True,
            "eligible": True,
            "reason_codes": ("built_from_imine_bridge", "valid_polyfunctional_amine_monomer", "external_sulfur_required"),
            "required_external_conditions": ("sulfur_enabled",),
            "canonical_built_linkage": "imine_bridge",
            "changes_applied": False,
        }
    }
}
```

Important rule:

- this annotation must not overwrite `candidate.events`
- this annotation must not overwrite `graph_summary["reaction_templates"]`
- this annotation must not make CIF export pretend the framework is already sulfur-converted

## Detection Plan

### 1. Shared strict aldehyde-site helper

Primary file: [`src/cofkit/chem/rdkit.py`](../src/cofkit/chem/rdkit.py)

Implement or centralize:

- `detect_aldehyde_sites(mol)`

This helper should be the single strict definition used by:

- the existing aldehyde motif builder
- any Stage 2 conversion-candidate checks that need aldehyde provenance

The conservative aldehyde rule remains:

- carbonyl carbon with one double bond to oxygen
- one hydrogen
- one scaffold substituent

Reject:

- ketones
- carboxylic acids
- esters
- amides
- acyl halides
- acetals
- hydrates

### 2. Strict imine-forming primary aryl amine helper

Primary file: [`src/cofkit/chem/rdkit.py`](../src/cofkit/chem/rdkit.py)

Add a helper such as:

- `detect_imine_eligible_primary_amine_sites(mol)`

This helper should enforce:

- neutral exocyclic primary amine
- exactly two hydrogens on the reactive nitrogen before condensation
- nitrogen attached to an aromatic or heteroaromatic `sp2` carbon
- site not embedded in excluded functional classes

Recommended behavior:

- use this helper to classify Stage 2 sulfur-conversion eligibility
- optionally reuse it for a future stricter imine auto-detect mode

### 3. Do not overload the generic `amine` motif kind without review

File: [`src/cofkit/chem/motif_registry.py`](../src/cofkit/chem/motif_registry.py)

The current `amine` motif kind is shared by:

- `imine_bridge`
- `keto_enamine_bridge`

If the generic `amine` detector is tightened globally, that may change current behavior.

Recommended conservative approach:

- keep the existing `amine` motif kind for Stage 1 building
- add stricter Stage 2 eligibility classification on top of it
- only tighten the shared detector globally after regression review

This avoids unintended changes to the current binary-bridge pipeline.

### 4. Monomer-level eligibility helpers

Recommended new helpers:

- `is_imine_cof_amine_monomer(monomer)`
- `is_imine_cof_aldehyde_monomer(monomer)`

Conservative amine-monomer rule:

- motif kind is `amine`
- at least `2` motifs
- every motif passes the strict imine-eligible primary aryl amine rule

Conservative aldehyde-monomer rule:

- motif kind is `aldehyde`
- at least `2` motifs
- every motif passes the strict aldehyde rule

## Stage 1 Build Plan

### Reuse the existing imine builder

Relevant files:

- [`src/cofkit/batch.py`](../src/cofkit/batch.py)
- [`src/cofkit/engine.py`](../src/cofkit/engine.py)
- [`src/cofkit/cli_build.py`](../src/cofkit/cli_build.py)

Stage 1 sequence should be:

1. detect polyaldehyde monomer sites
2. detect poly(primary aryl amine) monomer sites
3. pair compatible monomers by topic count and reticulation mode
4. build the imine-linked COF graph
5. store the imine linkage explicitly as the initial canonical linkage
6. annotate the built imine COF as eligible for sulfur-enabled conversion if the Stage 2 rules pass

### No new Stage 1 workflow family is needed

The imine build itself stays in the existing `binary_bridge` family.

The new logic belongs after candidate assembly, not in pair generation.

## Stage 2 Post-Build Conversion Eligibility Plan

### Conversion-candidate rule

A built imine COF is a sulfur-conversion candidate if:

- it was built from valid polyaldehyde and poly(primary aryl amine) monomers
- the imine linkage is present as the actual built linkage
- the framework does not rely on a different incompatible ring-forming linkage template
- the conversion chemistry is explicitly enabled as a second-stage transform
- sulfur is treated as an external reagent or condition

Do not require:

- pre-existing `-SH` on the amine monomer
- pre-existing sulfur on the imine nitrogen
- ortho-thiol substitution on the amine monomer

### Recommended assessment function

Suggested interface:

- `assess_sulfur_enabled_imine_conversion(candidate, monomer_specs, *, requested: bool)`

Inputs should examine:

- `candidate.events`
- `candidate.metadata["instance_to_monomer"]`
- monomer motif metadata
- built template ids

Outputs should include:

- `eligible: bool`
- `reason_codes`
- `required_external_conditions`
- `canonical_built_linkage`
- `changes_applied: False`

### Important representation rule

The canonical built representation remains `imine_bridge`.

Current implementation behavior:

- candidate graph construction stays on the imine path
- post-build conversion metadata records whether sulfur-enabled benzothiazole-like annulation is eligible
- when that profile is explicitly requested, atomistic CIF realization injects one external sulfur atom per eligible imine event and removes the aldehydic and ortho-aryl hydrogens needed for ring closure
- the graph summary and stored reaction-event template counts remain imine-based

The first milestone stops at eligibility assessment.

## Conversion Notes

The current atomistic conversion is intentionally conservative:

- sulfur is modeled as an external reagent atom added only during realization
- the conversion is heuristic and geometry-preserving, not a full local reoptimization
- the converted output is best treated as a chemically informed starting structure for later optimization

A future `calculate`-namespace workflow is still the right place for:

- local relaxation of the sulfur-converted linkage
- explicit external-tool optimization after conversion
- alternative sulfur-enabled conversion chemistries beyond the current benzothiazole-like annulation model

## CLI and API Plan

### Initial milestone

No new build-time template CLI is needed.

The existing build surface remains:

- `cofkit build single-pair --template-id imine_bridge ...`
- `cofkit build batch-binary-bridge --template-id imine_bridge ...`

Add only an opt-in annotation switch later, for example:

- `--annotate-post-build-conversions sulfur_enabled_imine_conversion`

### Future workflow split

- `build` should construct the imine COF
- `analyze` can report which outputs are eligible for sulfur-enabled conversion
- `calculate` is the better home for any future explicit sulfur-conversion transform

## Implementation Phases

### Phase 0: rewrite assumptions

Files:

- this plan document

Deliverables:

- remove the preinstalled `-SH` assumption
- define the workflow as Stage 1 imine build plus Stage 2 conversion eligibility

### Phase 1: strict eligibility helpers

Files:

- [`src/cofkit/chem/rdkit.py`](../src/cofkit/chem/rdkit.py)
- [`tests/test_rdkit_monomer.py`](../tests/test_rdkit_monomer.py)

Deliverables:

- shared strict aldehyde-site helper
- strict imine-eligible primary aryl amine helper
- positive and negative detector coverage

### Phase 2: monomer-level imine-COF eligibility

Files:

- new helper module, likely `src/cofkit/post_build_conversions.py`
- [`tests/test_chem.py`](../tests/test_chem.py) or a new dedicated test module

Deliverables:

- `is_imine_cof_amine_monomer`
- `is_imine_cof_aldehyde_monomer`
- clear rejection reasons for ineligible monomers

### Phase 3: post-build conversion annotation

Files:

- new `src/cofkit/post_build_conversions.py`
- [`src/cofkit/engine.py`](../src/cofkit/engine.py)
- [`src/cofkit/batch.py`](../src/cofkit/batch.py)
- new `tests/test_post_build_conversions.py`

Deliverables:

- conversion profile registry
- `sulfur_enabled_imine_conversion` profile
- candidate metadata annotation with explicit opt-in

### Phase 4: CLI and manifest integration

Files:

- [`src/cofkit/cli_build.py`](../src/cofkit/cli_build.py)
- [`tests/test_cli.py`](../tests/test_cli.py)
- batch summary/manifest writers in [`src/cofkit/batch.py`](../src/cofkit/batch.py)

Deliverables:

- optional CLI flag for conversion-candidate annotation
- manifest entries that record conversion eligibility
- no change to default CLI behavior when the flag is absent

### Phase 5: explicit conversion transform refinement

Files:

- likely new `calculate`-namespace code
- follow-up reaction-realization and CIF export refinements

Deliverables:

- local relaxation of the realized sulfur-containing linkage
- validation against targeted reference chemistry
- optional export modes that distinguish canonical imine graphs from converted atomistic realizations

## Test Matrix

### Detection tests

- aryl primary amine positive
- heteroaryl primary amine positive
- multitopic aryl amine positive
- secondary amine negative
- tertiary amine negative
- hydrazide negative
- sulfonamide negative
- protected amine negative
- aromatic ring `N` without exocyclic `NH2` negative
- aldehyde positive
- ketone negative

### Build-path tests

- valid amine + aldehyde pair still builds a normal imine COF
- `graph_summary["reaction_templates"]` remains `{"imine_bridge": ...}`
- CIF export still realizes imine, not sulfur-converted output

### Conversion-annotation tests

- valid imine COF is eligible only when conversion annotation is explicitly requested
- valid imine COF is not relabeled when annotation is absent
- non-imine candidate is ineligible
- candidate built from incompatible monomer classes is ineligible
- sulfur is recorded as an external condition, not a monomer feature

### Regression tests

- existing `imine_bridge` CLI behavior unchanged
- existing batch discovery unchanged
- existing binary-bridge summaries unchanged when annotation is disabled
- existing reaction realization unchanged for all current templates

## Recommended First Coding Milestone

The safest first coding milestone is:

1. implement strict imine-eligible primary aryl amine classification helpers
2. add a post-build conversion profile registry
3. annotate imine-built candidates as sulfur-conversion-eligible without changing the built graph

Do not start with a sulfur-containing transformed product.

The highest-value conservative milestone is proving that the code can:

- build the same imine COF it already builds today
- identify when that built framework is a valid sulfur-conversion candidate
- keep Stage 1 and Stage 2 cleanly separated
