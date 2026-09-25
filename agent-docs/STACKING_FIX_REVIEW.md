# Independent review of STACKING_FIX_PLAN.md

Reviewed 2026-09-24 against commit
`2e241794a21f0fefcb2aae0206467a7311ed9635`.

**Verdict:** the four main defects are real. The span, registry, and residual
fixes address demonstrated failures, but the plan needs revisions before it is
ready to implement. In particular, the proposed general clash rule produces
false positives, and measuring along c does not correctly define layer
thickness in the fitted cells currently produced by the code.

## Independent evidence

I generated new TAPB/TFB hcb and boroxine hcb candidates from SMILES through the
Python APIs, realized their products, exported fresh CIFs, and measured contacts
independently of the validator. I did not reuse the earlier audit's scripts,
CIFs, or optimization conclusions. Interlayer distances were computed from
exported fractional coordinates and the full cell matrix, enumerating periodic
translations in [-2, 2]^3 for these particular cells.

Reproduction commands (from the repository root):

```bash
uv run python out/stacking_independent_review/probe.py
uv run python out/stacking_independent_review/extra_probes.py
uv run pytest -q tests/test_ring_forming.py tests/test_validation.py tests/test_batch.py
```

The scripts, JSON results, and regenerated CIFs are local ignored evidence in
`out/stacking_independent_review/`. The existing suite above passed: **84 tests
in 9.25 s**. These passing tests do not cover the failures below.

| Claim | Independent result | Assessment |
|---|---|---|
| C1: missing span permits severe interlayer contacts | Fresh TAPB/TFB product span is 3.586133 Å; stacking records 0.0. AA has c = 6.8 Å and a 2.024706 Å interlayer C···C contact. | Confirmed. |
| C3: validator misses those contacts | AA, AB, and slipped all classify as `valid`, with no reasons and `min_nonbonded_heavy_distance = None`. Their minimum heavy interlayer contacts are 2.024706, 2.094295, and 2.171673 Å. | Confirmed. `None` means nothing found below 1.05 Å, not absence of clashes. |
| C2: 120° hexagonal registry uses the wrong shift | Fresh boroxine hcb has a = 15.117082 Å. Current AB shift has length 5.039065 Å; proposed shift has length 8.727851 Å. | Confirmed for the intended honeycomb vertex/pore registry. |
| C2/T6: fitted nearly hexagonal cell falls through | Fresh binary cell has gamma = 60.086662° and is labeled `oblique`; its AB shift consequently becomes (1/2, 1/2). A synthetic equal-length 60.09° cell also returns `oblique`. | Confirmed. |
| C4: residual normalization changes | Mean bridge residual changes from 0.005796032 to 0.002898016 upon stacking. | Confirmed: exactly a factor of two. |
| T1: nonzero fractional ring-center height breaks duplication | Translating a valid ring layer and its centers together by 0.2 c preserves its base geometry. Stacking introduces a 0.680380 Å ring planarity error. | Confirmed latent defect; default ring centers have w = 0. |
| Silent eligibility skip | Removing `net_plan` and requesting AA returns the original object, with no stderr diagnostic or added skip reason. | Confirmed. |
| T2: torsion copying | Code search finds passthrough copies, no active torsion lookup using instance/event keys. | Latent contract question; no demonstrated current numerical failure. |

Supplying the independently measured product span through the existing metadata
seam changes AA c to **13.972267 Å**. The minimum all-atom interlayer contacts
become **4.291525 / 4.188987 / 4.299418 Å** for AA/AB/slipped respectively. This
isolated experiment verifies the span correction; it does not implement the
registry correction at the same time. It also does not establish an equilibrium
spacing or validate the old optimization evidence.

## Changes required in the plan

### 1. W4.1 needs bond-graph-aware exclusions before widening clash detection

The current validator excludes direct bonded images only
(`validation.py:_min_nonbonded_heavy_distance_below_cutoff`). Applying the
proposed vdW rule to that pair set flags normal bonded-angle neighbors.

Independent counterexample: an ideal isolated benzene ring with 1.4 Å bonds has
1–3 C···C distances of **2.424871 Å**, below the proposed carbon threshold
**0.75 × (1.7 + 1.7) = 2.55 Å**. The existing search with that threshold finds
and would flag them, despite there being no intermolecular contact at all.

For the urgent stacking fix, use a layer-aware interlayer check. For general
validation, define periodic-image-aware 1–2 and 1–3 exclusions, an explicit 1–4
policy, and a separate severe-overlap check for excluded pairs. Specify the
radius source, supported-element/fallback policy, and a search radius sufficient
for the selected pair thresholds. Reporting hydrogen contacts separately is
reasonable. The existing `soft_relax.py` graph/image handling provides a useful
reference, although its force-field radii are not automatically interchangeable
with vdW radii. Add negative tests for ideal rings and ordinary angle neighbors;
a single clashing-pair positive test is insufficient. Do not schedule automatic
promotion to failure based solely on elapsed releases.

### 2. W1.2 and W4.4 should use the layer normal, or explicitly normalize the cell

The freshly generated binary cell already has c tilted **0.524955°** from the
normal to the ab plane. This is not exclusively a hypothetical malformed-input
case. Its a vector projects **-0.150668 Å** onto c-hat. Moving an atom to an
equivalent a-periodic image therefore changes its c-hat coordinate while leaving
its coordinate along the layer normal unchanged. A c-hat span is not invariant
to in-plane periodic representation in this cell.

Define `n = normalize(a × b)` and measure product coordinates along n. Either
construct the new c along n, documenting the cell change, or preserve its tilt
and solve for the required normal separation. In the latter case, the normal
repeat is `abs(c · n)/2`, not `|c|/2`; lateral components also affect successive
registries. A warning alone cannot make the advertised clearance derivation true.
Move this handling into M1, along with finite/degenerate-cell checks. Test rigid
rotations, tilted c, and equivalent in-plane images.

### 3. W2 fixes the reproduced honeycomb cases; limit its general claims

The proposed 120° shift passes a geometric check on the newly generated ring
net: one shifted vertex lies within **0.000087 Å** of a pore-center lattice
point, versus **5.039065 Å** for the current shift. The classifier tolerance
change also covers the measured 60.086662° case.

However, a metric family alone does not identify pore positions for every net
with that metric. Define the honeycomb Bernal and square hollow-site claims for
the topologies where they are tested. For other nets, describe AB as a specified
fractional translation or use topology-specific registry data. Explicitly handle
oblique, rectangular, unknown, and supercell settings rather than silently
claiming the square default has the same meaning. If claiming exact invariance
under basis changes, test transformed Cartesian geometry modulo lattice
translations; equal shift lengths alone are insufficient. The two symmetry-related
Bernal directions also need a stated convention.

### 4. W4.2 must check the realized periodic output, with a defined contact policy

`c2c > span` follows automatically from `c2c = positive_clearance + span`; it
cannot independently detect an underestimated span. The useful safety check is
on product atoms, using the same realization as export where possible, including
both galleries across the periodic boundary. In the span-corrected AA experiment,
the closest contact uses image **(0, 0, -1)**. Checking only the displayed pair of
layers would miss that minimum.

Use the same calibrated clash criterion as validation: “sub-vdW” alone is
ambiguous and could mean a stricter threshold than W4.1. Record whether each
minimum includes hydrogen, its atom labels and image, and whether it is an exact
minimum or a cutoff bound. A monomer-center broad phase needs conservative atom
extent bounds so it cannot discard contacting large monomers. Add a regression
whose clash exists only across the c boundary.

M1's claim that it “removes impossible outputs” is too strong if the specified
behavior is to emit them with warnings. State instead that the demonstrated span
error is fixed and remaining detected clashes are explicitly reported.

### 5. W3 must distinguish construction provenance from final artifact measurements

The richer schema and compatibility aliases are useful. But the existing batch
flow can run `soft_relax` after CIF export, changing positions without rebuilding
candidate stacking metadata. LAMMPS can change both coordinates and cell.
Construction span/contact values must not be presented as measurements of those
modified artifacts.

Give the block an explicit geometry stage, preserve initial construction values
as provenance, and recompute or invalidate current contact/spacing measurements
after modifications. `lammps.py:_render_optimized_cif` currently creates a fresh
header containing the COFid and optional warning; a second comment block does
not automatically survive. Test actual preservation/invalidation behavior, not
only continued parsing of the first COFid line.

### 6. W5.1 is sufficient for the demonstrated ranking bug; narrow the test wording

Doubling the bridge residual total restores **0.005796032**, exactly the base
mean. Compare normalized residual and the unreacted-motif tie-breaker, not the
entire `candidate_ranking_key` tuple: the last component intentionally changes
with the candidate ID. Also test ordering against an unexpanded candidate.

The error does not reorder candidates if all receive the same factor-of-two
change, so “critical ranking bias” overstates its effect within an exclusively
stacked binary candidate set. Mixed expanded/unexpanded comparisons and exported
residual values are affected. Separating topology ranks from registry variant
indices is sensible and should be explicit if those fields promise topology
counts.

### 7. W1.4, W6 and W7 need scope and acceptance-test adjustments

- A periodic one-layer ring model with c = 3.4 Å can represent an implicit AA
  repeat. Its difference from an 8 Å slab convention is not itself a geometry
  bug. Document distinct modes; make a default change a deliberate compatibility
  decision, separate from the urgent stacking repair. The span-corrected bilayer
  is a conservative starting model, not evidence of a physically correct bulk
  repeat.
- For T1, prefer converting the old Cartesian ring center into the new basis
  over asserting w = 0. The translated-layer experiment shows that w != 0 can
  describe an entirely valid geometry. Include existing center offsets in the
  conversion and ensure they are added exactly once. Recompute relevant geometry
  checks rather than trusting copied validation metadata.
- Define the torsion key contract before mechanically suffixing every key.
  Current code provides no evidence that arbitrary torsion keys are instance IDs.
- Handle realization failures, partial atom coverage, nonfinite coordinates,
  invalid clearances, and absent specs explicitly. Scope exception isolation to
  each candidate/variant: batch stacking expansion currently occurs in one tuple
  comprehension, so a new measurement/self-check exception can discard otherwise
  usable siblings if not handled there.
- Replace the proposed grep gate for `0.0` with behavioral tests: a genuinely
  planar measured layer legitimately has zero span. Assert measurement provenance
  and complete atom coverage. Distinguish coarse/mixed estimates from atomistic
  measurements.
- Turn regenerated evidence into self-contained fixtures or deterministic builders
  under `tests/`; CI cannot depend on ignored `out/stacking_audit_evidence` files.
  Add negative clash controls, periodic-boundary contacts, translated ring centers,
  tilted-cell cases, and failure-continuation tests.

## Suggested implementation order

1. Correct product-span measurement and the normal-axis contract, add a realized
   periodic interlayer diagnostic, and fix residual duplication. Include explicit
   skip/failure diagnostics and per-candidate exception isolation in this step.
2. Correct and test honeycomb registry settings, with supported-topology semantics.
3. Add artifact/provenance fields and graph-aware general validation, including
   negative controls and post-processing invalidation/recomputation.
4. Decide single-layer defaults and longer-term registry/relaxation exploration
   separately.

No runtime implementation was changed during this review. No LAMMPS calculations
were run; the claimed optimized gallery pairing and its physical interpretation
remain unverified here.
