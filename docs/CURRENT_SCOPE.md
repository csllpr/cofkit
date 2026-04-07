# Current Scope

The base package now installs `rdkit` and `gemmi` as mandatory runtime dependencies because the practical monomer-construction, CIF, and topology workflows depend on them. External binaries such as Zeo++, LAMMPS, EQeq, and gRASPA remain optional add-ons.

## Implemented so far

- geometry primitives and local frames
- core domain dataclasses
- builtin COF reaction library plus linkage-profile metadata for pair-role ordering, bridge targets, and realization hooks
- periodic product graph with motif usage validation
- topology repository/index support for RCSR CGD-style bundles, with bundled topology data preferred by default and builtin fallback hints
- a discrete assignment layer for motif-to-reaction event matching
- registry-backed motif metadata plus lightweight geometric motif detection for the currently supported fallback kinds
- RDKit/SMARTS-backed monomer construction for amine, aldehyde, hydrazine, hydrazide, boronic acid, catechol, keto aldehyde, and activated-methylene motifs, including conformer generation and bond retention in the standard install
- initial linkage geometry helpers for bridge-forming reactions
- initial periodic embedding for monomer instances from topology hints and motif/reaction heuristics, including oblique `hcb` cells for asymmetric `3+3` and `3+2` cases when a symmetric hexagonal metric is too restrictive
- a dependency-free continuous optimization pass for seed cell/pose refinement after embedding
- first-pass candidate scoring with event coverage, bridge geometry, topology bonuses, and unreacted penalties
- candidate metadata carrying embedding provenance, optimizer metrics, and score breakdowns
- legal P1 CIF export, with atomistic output when monomer coordinates are available and a coarse fallback otherwise, now including explicit `_ccdc_geom_bond_type` records for atomistic bond loops
- registry-backed atomistic reaction realization for the currently implemented binary-bridge products: imine, hydrazone, azine, beta-ketoenamine, boronate ester, and vinylene
- periodic-image-safe imine atomistic realization, with the same bent-linkage motif-origin correction now applied across the supported `2D` and `3D` builder paths
- reaction-aware batch binary-bridge generation over monomer libraries, with imine workflows as the primary documented path, including `3+3`, `3+2`, `4+4`, `4+2`, and `6+2` enumeration, manifest/summary writing, and CIF export enabled by default
- automatic monomer-role detection for batch library loading, so generic `.txt` SMILES libraries can be regrouped by detected role/connectivity instead of relying only on `*_count_N.txt` filenames
- an auto-generated example library under [`../examples/default_monomers_library`](../examples/default_monomers_library) with detector-scanned role/count registration metadata
- symmetry-expanded single-node generation across supported `2D` and `3D` one-node families, available through both batch runs and single-pair input, with `3+2` / `3+3` handling on `hcb` / `hca` / `fes` / `fxt`, plus `dia` and `pcu` builders for `4+2`, `4+4`, and `6+2` cases
- default topology-family routing that now sends higher-connectivity nonplanar inputs to `3D` one-node nets (`dia` for `4`-connected cases, `pcu` for `6+2`) while keeping older `2D` high-connectivity families opt-in
- shared topology-builder dispatch for the supported one-node families, reused by batch runs and direct single-pair generation
- indexed-topology layout generation for chemistry-compatible bundled topologies beyond the handcrafted one-node families, available in both batch and explicit single-pair generation
- compatibility-aware default topology selection that now includes curated indexed topologies such as `sql`, `kgm`, `hxl`, `pts`, `ctn`, `bor`, `kgd`, `tbo`, `dia`, `pcu`, `acs`, `lon`, and `qtz` when the current chemistry metadata and builder support permit them
- coarse post-generation validation with `valid` / `warning` / `hard_invalid` / `hard_hard_invalid` triage, categorized CIF output trees, and export blocking for obviously broken structures
- an initial `cofkit analyze zeopp` wrapper for Zeo++ pore analysis from CIF input, including a point-probe baseline from `-res` / `-resex` / `-chan 0` / `-sa 0 0` / `-vol 0 0` plus optional repeated probe-radius scans with accessibility summaries, with binary discovery through `COFKIT_ZEOPP_PATH`
- an initial `cofkit calculate lammps-optimize` wrapper for LAMMPS optimization of explicit-bond `P1` CIFs, including env-based binary discovery through `COFKIT_LMP_PATH`, nearest-image bond-length inference when CIF bond rows omit periodic image hints, required explicit bond-type reuse from CIF `_ccdc_geom_bond_type`, library-backed LAMMPS data-file generation, optional short restrained preruns before minimization, staged minimization support, optional final `fix box/relax`, and updated CIF export
- an initial `cofkit calculate graspa-widom` wrapper for staged `EQeq -> gRASPA` Widom insertion on one CIF, including env-based executable discovery through `COFKIT_EQEQ_PATH` and `COFKIT_GRASPA_PATH`, packaged Widom-template asset materialization, automatic `UnitCells` sizing from CIF cell lengths plus cutoff settings, parsed `results.csv` generation, and `graspa_widom_report.json` output
- a working UFF-backed LAMMPS parameterization path using Open Babel UFF atom typing, the bundled pinned Open Babel `UFF.prm` reference, `pymatgen` LAMMPS-data generation for bond, angle, dihedral, improper, and van der Waals terms, and default EQeq charge staging for charged `atom_style full` exports; this is the default and currently only implemented public LAMMPS force field, and local `spring/self` restraints now contribute to the minimized objective through `fix_modify energy yes`
- process-level batch pair generation with an `8`-worker default budget for the practical CLI workflows
- an installable `cofkit` CLI, organized under grouped `build` / `analyze` namespaces for direct `single-pair` generation plus unified batch, classification, and library-building entry points
- extracted batch-facing support layers (`cofkit.monomer_library`, `cofkit.batch_models`) so CLI, wrappers, and future linkage extensions share the same monomer-role resolution and summary schema
- tests covering core invariants

## Not implemented yet

- broader SMARTS/rule-based motif detection beyond the current binary-bridge set, especially for ring-forming and cyclization chemistries
- symmetry reduction beyond raw indexed topology filtering
- full general topology coverage beyond the current supported one-node families plus the present indexed-layout subset
- torsion-aware or force-field-backed optimization beyond the current lightweight pass
- stable ring-closure geometry models suitable for public CLI exposure
- a supported public benzothiazole conversion workflow; the current sulfur-enabled imine conversion prototype remains internal-only until its local geometry is more reliable
- Zeo++ PSD histograms, grid outputs, ray analyses, ZeoVis exports, and broader hidden or specialized workflows beyond the current pore-summary wrapper
- fuller LAMMPS force-field coverage beyond the current UFF-backed path, especially validated charge models and broader force-field families
- broader gRASPA workflow coverage beyond the current bundled Widom-template family, bundled component set, and Henry-coefficient / Widom-energy result parser
- non-`P1` bonded CIF optimization beyond the current conservative local cleanup wrapper
- chemically faithful atomistic CIF generation for arbitrary monomers without fallback/pseudo-sites
- semiempirical or quantum-chemistry cleanup beyond the current external-tool wrappers
- any stacking exploration, registry search, or stacking score terms
