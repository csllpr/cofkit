# Documentation

Use this directory as the canonical documentation set for `cofkit`.

## User Workflows

- [Jupyter tutorial series](../tutorials/README.md): executable CLI notebooks with commented Python API equivalents
- [Getting started](getting-started.md): repository setup, verification, first build, CLI shape
- [Building COFs](building.md): single-pair generation, batch generation, all binary bridges, default library rebuilds, topology and stacking options
- [Analysis and validation](analysis-validation.md): output classification, CIF-to-COFid decomposition, COFid validation, Zeo++ pore analysis
- [Calculations](calculations.md): LAMMPS optimization, EQeq charge staging, gRASPA/RASPA2 Widom/isotherm/mixture workflows, hybrid MD/MC
- [Python API](python-api.md): practical `COFEngine` and `BatchStructureGenerator` examples
- [External tools](external-tools.md): recommended Zeo++, LAMMPS, EQeq, gRASPA, and RASPA2 installs

## Reference

- [CURRENT_SCOPE.md](CURRENT_SCOPE.md): implemented capabilities and current limitations
- [TECHNICAL_OVERVIEW.md](TECHNICAL_OVERVIEW.md): design notes, topology repository details, and pipeline structure
- [COARSE_VALIDATION.md](COARSE_VALIDATION.md): validation buckets and current thresholds
- [COFid_Specification_v1.2.md](COFid_Specification_v1.2.md): COFid format reference
- [WRAPPER_SCRIPTS.md](WRAPPER_SCRIPTS.md): legacy notes for the thin scripts under `examples/`

## Related

- [../ARCHITECTURE.md](../ARCHITECTURE.md): broader architecture notes
- [../CHANGELOG.md](../CHANGELOG.md): release history
- [../agent-docs/README.md](../agent-docs/README.md): development and agent-facing guides
- [../skills/cofkit-navigator/SKILL.md](../skills/cofkit-navigator/SKILL.md): Codex workflow routing skill
