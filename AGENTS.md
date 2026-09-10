# AGENTS.md

Entrypoint for coding agents working in this repository. `cofkit` is a
reaction-aware periodic assembly toolkit for covalent organic frameworks
(COFs): monomer detection, reaction templates, topology-guided assembly,
CIF export, validation, and external-simulation wrappers.

## Documentation layout

- `docs/` — traditional user-facing documentation (install, workflows,
  CLI usage, tutorials entry points).
- `agent-docs/` — development documentation written for agents:
  - `agent-docs/AGENT_CODEBASE_MAP.md` — authoritative module/seam map;
    consult it before editing anything non-trivial.
  - `agent-docs/ADDING_LINKAGES_AND_MONOMERS.md` — chemistry-extension
    checklist.
- `ARCHITECTURE.md` — pipeline, batch flow, topology flow, design
  principles.
- `docs/CURRENT_SCOPE.md` — what is and is not implemented; do not claim
  or route users to unimplemented features.
- `skills/cofkit-navigator/SKILL.md` — for operating the CLI on behalf of
  users (intent → command routing, output interpretation).

Keep summaries short here; the files above carry the detail.

## Environment and commands

Use `uv`; the lockfile is canonical.

```bash
uv sync --locked --extra dev   # full dev environment
uv run pytest -q               # test suite (testpaths = tests)
uv run cofkit --help           # CLI smoke check
uv lock --check                # lockfile consistency (CI gate)
uvx ruff==0.12.12 check src --select E9,F63,F7,F82   # lint gate (CI)
uv build --wheel               # packaging check (CI)
```

- Python >= 3.10, `src/` layout, hatchling build backend.
- Runtime deps: `rdkit`, `gemmi`, `openbabel`, `pandas`, `pymatgen` — all
  mandatory; do not add ASE (decomposition intentionally uses `gemmi`).
- CI also runs `tests/engines` against a real LAMMPS wheel; locally these
  tests skip unless `COFKIT_TEST_LMP` / `COFKIT_REQUIRE_LAMMPS` are set.
- CI runs the test suite outside the checkout against the installed
  wheel, so never rely on cwd-relative imports of the source tree in
  tests.

## Where to make changes

- New linkage/reaction: start in `src/cofkit/reactions.py`, then
  `src/cofkit/chem/motif_registry.py`, `src/cofkit/chem/rdkit.py`, and
  `src/cofkit/reaction_realization.py`. Full checklist:
  `agent-docs/ADDING_LINKAGES_AND_MONOMERS.md`.
- CLI commands: grouped routers `cli_build.py` / `cli_analyze.py` /
  `cli_calculate.py` / `cli_validate.py` under `src/cofkit/`; legacy flat
  aliases must keep working with deprecation warnings.
- Batch-facing data: extend `src/cofkit/batch_models.py` dataclasses
  instead of introducing ad hoc summary dicts.
- Tests: mirror the existing `tests/test_*.py` layout
  (`agent-docs/AGENT_CODEBASE_MAP.md` maps seams to test files).

## Constraints and conventions

- The practical generation path is binary-bridge-first; ring-forming
  supports only boroxine/triazine. Check `cofkit build list-templates`
  (`supports_pair_generation`) before assuming a template can build.
- The benzothiazole conversion prototype is internal-only — never expose
  it through the public CLI.
- DREIDING is the recommended force field; UFF is experimental.
- Versioning is calendar-based `YYYY.M.D` (post-releases like
  `2026.4.7.post1` for same-day re-releases).
- Keep scoring and optimization separate even when they share metrics;
  keep seed assembly honest about not being physical relaxation.

## Repository hygiene

- `out/`, `files_for_reference/`, `.untracked/`,
  `reference_repositories/` are gitignored scratch/output trees — never
  commit generated artifacts.
- Never commit `.env`. The CLI auto-loads the nearest `.env` upward from
  the cwd; shell-set variables take precedence.
- External tools are optional and configured via env vars:
  `COFKIT_ZEOPP_PATH`, `COFKIT_LMP_PATH`, `COFKIT_EQEQ_PATH`,
  `COFKIT_GRASPA_PATH`, `COFKIT_RASPA2_PATH`. Build/validate workflows
  must keep working without them.

## When behavior changes

Update in the same change:

- `CHANGELOG.md`
- `docs/CURRENT_SCOPE.md` (implemented / not-implemented lists)
- `ARCHITECTURE.md` and/or `agent-docs/AGENT_CODEBASE_MAP.md` if module
  responsibilities or seams moved
- `skills/cofkit-navigator/SKILL.md` if CLI behavior, defaults, or
  output artifacts changed
