# cofkit Action Checklist

_Last updated: 2026-03-14 (Asia/Shanghai)_

## Completed

- [x] Review PORMAKE and decide not to use it as the scientific core
- [x] Create `cofkit` and land the initial reaction-aware scaffold
- [x] Commit the first milestone baseline (`0eb185d`)
- [x] Add explicit bond-loop CIF export for atomistic outputs
- [x] Preserve monomer bond connectivity from RDKit / detected monomers into `MonomerSpec`
- [x] Merge retained intramonomer bonds with realized inter-monomer reaction bonds during CIF export
- [x] Deduplicate exported CIF bond records
- [x] Extend tests for atomistic CIF bond-loop export and reacted RDKit imine output
- [x] Verify current repo state still passes tests
- [x] Verify the new reacted CIF file can be sent out and checked externally

## Maintenance completed today

- [x] Recover yesterday's project state from archived transcript
- [x] Create a durable workspace recovery note outside the repo
- [x] Write/update project changelog for the bond-export milestone
- [x] Create a clean commit for the CIF bond-export milestone

## Next up

### Section: optimizer / refinement extension
- [x] Review the current lightweight optimizer/scoring loop for missing chemistry constraints
- [x] Extend bridge-geometry evaluation beyond raw bond distance/alignment where justified
- [x] Add at least one focused refinement metric or control tied to bond-angle / planarity behavior
- [x] Add regression tests for the new refinement behavior
- [x] Update docs once the optimizer extension lands

### Later backlog
- [ ] Better motif detection beyond current lightweight heuristics
- [ ] Richer reaction-geometry models, especially for ring-forming events
- [ ] Stronger validation/ranking (clash, strain, planarity, torsion)
- [ ] More downstream-friendly atomistic export metadata
- [ ] Keep stacking explicitly out of scope until the core chemistry pipeline is stronger
