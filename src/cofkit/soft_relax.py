"""Prototype in-process clash-repair / strain-relief relaxer.

This is a deliberately rough, dependency-free optimizer intended to run on a
built structure before (or right after) CIF export so that graphs that are
topologically sound but geometrically clashing are not bucketed
``hard_invalid`` before the LAMMPS relaxation pipeline gets a chance at them.
It mirrors the staged ``pair_style soft`` pre-minimization already used by the
LAMMPS backend (``lammps.py``), but runs in-process with only gemmi:

- harmonic bond springs at the *initial* bond lengths (so realized linkage
  geometry from ``reaction_realization`` is preserved, not relaxed away),
- Urey-Bradley-style harmonic 1-3 distance springs at the initial 1-3
  distances, restraining angles against collapse without angle gradients,
- a ramped cosine soft-repulsion between non-bonded, non-1-3 atom pairs,
  with per-element radii from the DREIDING reference table,
- fixed unit cell; periodicity is handled through the explicit bond-image
  shifts in the ``_geom_bond`` loop plus ``periodic_geometry.images_within``.

It is NOT a force field and makes no claim of physical relaxation; the output
is a clash-relieved seed structure for downstream validation bucketing or a
real optimization backend.

Limitations of this prototype:

- CIF writing rewrites the ``_atom_site_fract_*`` columns of the atom-site
  loop textually and therefore expects the cofkit ``CIFWriter`` row format
  (one whitespace-separated row per atom with a unique label).
- Dihedral integrity is only preserved indirectly through the bond and 1-3
  restraint network; there are no torsion terms.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
import math
from pathlib import Path

import gemmi

from ._dreiding_reference import DREIDING_FRAMEWORK_TYPE_BY_ELEMENT, DREIDING_PARAMETERS
from .cif_checks import cif_value_str
from .periodic_geometry import images_within, p1_shift


@dataclass(frozen=True)
class SoftRelaxConfig:
    repulsion_scale: float = 0.8
    """r0_ij = repulsion_scale * (R_i + R_j); R from DREIDING LJ r0 / 2."""
    repulsion_coefficients: tuple[float, ...] = (1.0, 5.0, 20.0, 50.0, 100.0)
    """Soft-repulsion ramp, extending the LAMMPS soft pre-minimization staging
    with a final stiffer step for deep interlayer clashes."""
    bond_force_constant: float = 500.0
    """Harmonic bond spring (kcal/mol/A^2) at each initial bond length."""
    urey_bradley_force_constant: float = 50.0
    """Harmonic 1-3 distance spring (kcal/mol/A^2) at each initial 1-3 distance;
    restrains angles at their initial values without angle gradients."""
    neighbor_margin: float = 2.0
    """Extra pair-list margin (A) so the list can be built once."""
    repel_hydrogen: bool = True
    initial_step: float = 0.05
    max_displacement: float = 0.1
    """Per-step Cartesian displacement cap (A) for the adaptive descent."""
    max_steps_per_stage: int = 400
    pair_list_interval: int = 25
    """Rebuild the repulsion pair list this often; atoms can drift across cell
    boundaries into contacts that were not close when the list was built."""
    force_tolerance: float = 5.0
    """Stage-converged when max force falls below this (kcal/mol/A)."""
    fallback_radius: float = 2.0
    """Repulsion radius (A) for elements missing from the DREIDING table."""


@dataclass(frozen=True)
class SoftRelaxReport:
    input_path: str
    output_path: str
    n_atoms: int
    n_bonds: int
    n_pairs: int
    n_angle_restraints: int
    stage_energies: tuple[float, ...]
    converged: bool
    min_heavy_distance_before: float | None
    min_heavy_distance_after: float | None
    clashes_before: int
    clashes_after: int
    max_bond_drift: float
    warnings: tuple[str, ...] = field(default_factory=tuple)


@dataclass
class _System:
    labels: list[str]
    symbols: list[str]
    frac: list[list[float]]
    cell: gemmi.UnitCell
    orth: list[list[float]]
    # (i, j, shift) with i <= j; integer shift applied to j's fractional position.
    bonds: list[tuple[int, int, tuple[int, int, int]]]
    # 1-3 pairs (same encoding as bonds) used for Urey-Bradley angle restraints.
    angles13: list[tuple[int, int, tuple[int, int, int]]]
    excluded: set[tuple[int, int, tuple[int, int, int]]]
    bonded_label_images: set[tuple[str, str, tuple[int, int, int]]]


def _element_radius(symbol: str, config: SoftRelaxConfig) -> float:
    dtype = DREIDING_FRAMEWORK_TYPE_BY_ELEMENT.get(symbol)
    params = DREIDING_PARAMETERS.get(dtype) if dtype else None
    return params.r0 / 2.0 if params else config.fallback_radius


def _frac_to_cart(orth: list[list[float]], frac) -> list[float]:
    return [
        orth[0][0] * frac[0] + orth[0][1] * frac[1] + orth[0][2] * frac[2],
        orth[1][0] * frac[0] + orth[1][1] * frac[1] + orth[1][2] * frac[2],
        orth[2][0] * frac[0] + orth[2][1] * frac[1] + orth[2][2] * frac[2],
    ]


def _parse_system(block: gemmi.cif.Block, config: SoftRelaxConfig) -> tuple[_System, list[str]]:
    warnings: list[str] = []
    small = gemmi.make_small_structure_from_block(block)
    labels = [str(site.label) for site in small.sites]
    symbols = [str(site.element.name) for site in small.sites]
    frac = [[site.fract.x, site.fract.y, site.fract.z] for site in small.sites]
    index_by_label = {label: i for i, label in enumerate(labels)}

    bonds: list[tuple[int, int, tuple[int, int, int]]] = []
    bonded_label_images: set[tuple[str, str, tuple[int, int, int]]] = set()
    labels1 = block.find_loop("_geom_bond_atom_site_label_1")
    labels2 = block.find_loop("_geom_bond_atom_site_label_2")
    sym1 = block.find_loop("_geom_bond_site_symmetry_1")
    sym2 = block.find_loop("_geom_bond_site_symmetry_2")
    for k in range(min(len(labels1), len(labels2))):
        label_a = cif_value_str(labels1[k])
        label_b = cif_value_str(labels2[k])
        a = index_by_label.get(label_a)
        b = index_by_label.get(label_b)
        if a is None or b is None:
            warnings.append("geom_bond row references an unknown atom label; row skipped.")
            continue
        first = p1_shift(cif_value_str(sym1[k]) if len(sym1) else ".")
        second = p1_shift(cif_value_str(sym2[k]) if len(sym2) else ".")
        shift = tuple(sb - sa for sa, sb in zip(first, second))
        bonded_label_images.add((label_a, label_b, shift))
        bonded_label_images.add((label_b, label_a, tuple(-v for v in shift)))
        i, j = (a, b) if a <= b else (b, a)
        s = shift if a <= b else tuple(-v for v in shift)
        bonds.append((i, j, s))

    # Directed adjacency: adj[x] holds (y, t) meaning y at image t bonds x.
    adjacency: dict[int, set[tuple[int, tuple[int, int, int]]]] = {}
    for i, j, s in bonds:
        adjacency.setdefault(i, set()).add((j, s))
        adjacency.setdefault(j, set()).add((i, tuple(-v for v in s)))

    excluded: set[tuple[int, int, tuple[int, int, int]]] = set()
    for i, j, s in bonds:
        excluded.add((i, j, s))
        excluded.add((j, i, tuple(-v for v in s)))
    # Image-consistent 1-3 exclusion around every bonded center, kept both as
    # repulsion exclusions and as Urey-Bradley restraint pairs (deduplicated).
    angles13: list[tuple[int, int, tuple[int, int, int]]] = []
    seen13: set[tuple[int, int, tuple[int, int, int]]] = set()
    for center, neighbors in adjacency.items():
        ordered = sorted(neighbors)
        for x, (a, ta) in enumerate(ordered):
            for b, tb in ordered[x + 1:]:
                if a == b and ta == tb:
                    continue
                # Relative image of b seen from a is tb - ta.
                if (a, ta) <= (b, tb):
                    i, j = (a, b) if a <= b else (b, a)
                    shift = tuple(vb - va for va, vb in zip(ta, tb))
                    s = shift if a <= b else tuple(-v for v in shift)
                else:  # pragma: no cover - sorted() already orders them
                    continue
                if (i, j, s) in seen13:
                    continue
                seen13.add((i, j, s))
                excluded.add((i, j, s))
                excluded.add((j, i, tuple(-v for v in s)))
                if i != j or s != (0, 0, 0):
                    angles13.append((i, j, s))

    orth = [list(row) for row in small.cell.orth.mat.tolist()]
    system = _System(
        labels=labels,
        symbols=symbols,
        frac=frac,
        cell=small.cell,
        orth=orth,
        bonds=bonds,
        angles13=angles13,
        excluded=excluded,
        bonded_label_images=bonded_label_images,
    )
    return system, warnings


def _build_pair_list(
    system: _System, config: SoftRelaxConfig
) -> list[tuple[int, int, list[float], float]]:
    radii = [_element_radius(symbol, config) for symbol in system.symbols]
    if not config.repel_hydrogen:
        radii = [0.0 if s == "H" else r for s, r in zip(system.symbols, radii)]
    max_r0 = 2.0 * max(radii, default=config.fallback_radius) * config.repulsion_scale
    cutoff = max_r0 + config.neighbor_margin

    pairs: list[tuple[int, int, list[float], float]] = []
    n = len(system.frac)
    for i in range(n):
        for j in range(i, n):
            r0 = (radii[i] + radii[j]) * config.repulsion_scale
            if r0 <= 1e-9:
                continue
            for shift, _distance in images_within(system.cell, system.frac[i], system.frac[j], cutoff):
                if i == j and shift == (0, 0, 0):
                    continue
                if (i, j, shift) in system.excluded:
                    continue
                shift_cart = _frac_to_cart(system.orth, shift)
                pairs.append((i, j, shift_cart, r0))
                break  # nearest image only; the margin covers drift during the run
    return pairs


def _bond_distances(cart: list[list[float]], bonds_cart) -> list[float]:
    distances = []
    for (i, j), sc in bonds_cart:
        dx = cart[j][0] + sc[0] - cart[i][0]
        dy = cart[j][1] + sc[1] - cart[i][1]
        dz = cart[j][2] + sc[2] - cart[i][2]
        distances.append(math.sqrt(dx * dx + dy * dy + dz * dz))
    return distances


def _energy_and_forces(
    cart: list[list[float]],
    springs: list[tuple[tuple[int, int], list[float], float, float]],
    pairs: list[tuple[int, int, list[float], float]],
    repulsion_a: float,
    config: SoftRelaxConfig,
) -> tuple[float, list[list[float]]]:
    forces = [[0.0, 0.0, 0.0] for _ in cart]
    energy = 0.0

    for (i, j), sc, rest, k_spring in springs:
        bx = cart[j][0] + sc[0] - cart[i][0]
        by = cart[j][1] + sc[1] - cart[i][1]
        bz = cart[j][2] + sc[2] - cart[i][2]
        r = math.sqrt(bx * bx + by * by + bz * bz)
        if r < 1e-9:
            continue
        stretch = r - rest
        energy += 0.5 * k_spring * stretch * stretch
        fmag = k_spring * stretch / r
        fx, fy, fz = fmag * bx, fmag * by, fmag * bz
        forces[i][0] += fx
        forces[i][1] += fy
        forces[i][2] += fz
        forces[j][0] -= fx
        forces[j][1] -= fy
        forces[j][2] -= fz

    for i, j, sc, r0 in pairs:
        dx = cart[j][0] + sc[0] - cart[i][0]
        dy = cart[j][1] + sc[1] - cart[i][1]
        dz = cart[j][2] + sc[2] - cart[i][2]
        r2 = dx * dx + dy * dy + dz * dz
        if r2 >= r0 * r0:
            continue
        r = math.sqrt(r2)
        if r < 1e-6:
            # Exact overlap: deterministic pseudo-random push direction.
            h = (i * 73856093) ^ (j * 19349663)
            dx, dy, dz = float((h % 7) - 3), float(((h // 7) % 7) - 3), float(((h // 49) % 7) - 3)
            r = math.sqrt(dx * dx + dy * dy + dz * dz)
            fmag = repulsion_a * math.pi / (2.0 * r0) / r
            energy += repulsion_a
        else:
            x = math.pi * r / r0
            energy += repulsion_a * 0.5 * (1.0 + math.cos(x))
            fmag = repulsion_a * math.pi / (2.0 * r0) * math.sin(x) / r
        fx, fy, fz = fmag * dx, fmag * dy, fmag * dz
        forces[i][0] -= fx
        forces[i][1] -= fy
        forces[i][2] -= fz
        forces[j][0] += fx
        forces[j][1] += fy
        forces[j][2] += fz

    return energy, forces


def _small_structure(system: _System) -> gemmi.SmallStructure:
    small = gemmi.SmallStructure()
    small.cell = system.cell
    for label, symbol, frac in zip(system.labels, system.symbols, system.frac):
        site = gemmi.SmallStructure.Site()
        site.label = label
        site.element = gemmi.Element(symbol)
        site.fract = gemmi.Fractional(*frac)
        site.occ = 1.0
        small.sites.append(site)
    return small


def _min_heavy_distance_and_clashes(
    system: _System, clash_cutoff: float
) -> tuple[float | None, int]:
    """Mirror validation.py's clash check; also count sub-cutoff contacts."""
    cutoff = max(3.0, clash_cutoff)
    small = _small_structure(system)
    search = gemmi.NeighborSearch(small, cutoff).populate(include_h=False)
    minimum: float | None = None
    directed_clashes = 0
    for index, site in enumerate(small.sites):
        if site.element.is_hydrogen:
            continue
        candidates = {(index, 0)}
        candidates.update(
            (int(mark.atom_idx), int(mark.image_idx))
            for mark in search.find_site_neighbors(site, min_dist=0, max_dist=cutoff)
        )
        for other_index, image_index in candidates:
            other = small.sites[other_index]
            if other.element.is_hydrogen:
                continue
            position = other.fract
            if image_index:
                position = small.cell.images[image_index - 1].apply(position)
            for shift, distance in images_within(small.cell, site.fract, position, cutoff):
                if index == other_index and image_index == 0 and shift == (0, 0, 0):
                    continue
                if image_index == 0 and (site.label, other.label, shift) in system.bonded_label_images:
                    continue
                if distance < clash_cutoff:
                    directed_clashes += 1
                if minimum is None or distance < minimum:
                    minimum = distance
    return minimum, directed_clashes // 2


def _format_p1_shift(shift: tuple[int, int, int]) -> str:
    if shift == (0, 0, 0):
        return "."
    return "1_" + "".join(str(component + 5) for component in shift)


def _minimum_image_shift(
    cell: gemmi.UnitCell, frac_a, frac_b
) -> tuple[tuple[int, int, int], float]:
    """Shift s minimizing |orth(frac_b + s - frac_a)| over +-1 images."""
    best_shift = (0, 0, 0)
    best_distance = float("inf")
    pa = cell.orthogonalize(gemmi.Fractional(*frac_a))
    for dx, dy, dz in product((-1, 0, 1), repeat=3):
        pb = cell.orthogonalize(
            gemmi.Fractional(frac_b[0] + dx, frac_b[1] + dy, frac_b[2] + dz)
        )
        distance = pa.dist(pb)
        if distance + 1e-9 < best_distance:
            best_distance = distance
            best_shift = (dx, dy, dz)
    return best_shift, best_distance


def _rewrite_cif_coordinates(
    input_path: Path,
    output_path: Path,
    new_frac: list[list[float]],
    labels: list[str],
    cell: gemmi.UnitCell,
) -> None:
    new_by_label = {label: frac for label, frac in zip(labels, new_frac)}
    lines = input_path.read_text().splitlines(keepends=True)
    out_lines: list[str] = []
    atom_columns: list[str] = []
    bond_columns: list[str] = []
    mode = ""
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped == "loop_":
            j = i + 1
            tags: list[str] = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                tags.append(lines[j].strip().split()[0])
                j += 1
            if "_atom_site_fract_x" in tags and "_atom_site_label" in tags:
                mode = "atoms"
                atom_columns = tags
            elif "_geom_bond_atom_site_label_1" in tags and "_geom_bond_site_symmetry_2" in tags:
                mode = "bonds"
                bond_columns = tags
            else:
                mode = ""
            out_lines.extend(lines[i:j])
            i = j
            continue
        if mode and (not stripped or stripped.startswith(("#", "_", "loop_", "data_", "save_"))):
            mode = ""
        if mode == "atoms":
            tokens = stripped.split()
            label = cif_value_str(tokens[atom_columns.index("_atom_site_label")])
            frac = new_by_label.get(label)
            if frac is not None:
                wrapped = tuple(f - math.floor(f) for f in frac)
                tokens[atom_columns.index("_atom_site_fract_x")] = f"{wrapped[0]:.6f}"
                tokens[atom_columns.index("_atom_site_fract_y")] = f"{wrapped[1]:.6f}"
                tokens[atom_columns.index("_atom_site_fract_z")] = f"{wrapped[2]:.6f}"
                line = " ".join(tokens) + "\n"
        elif mode == "bonds":
            tokens = stripped.split()
            label_a = cif_value_str(tokens[bond_columns.index("_geom_bond_atom_site_label_1")])
            label_b = cif_value_str(tokens[bond_columns.index("_geom_bond_atom_site_label_2")])
            frac_a = new_by_label.get(label_a)
            frac_b = new_by_label.get(label_b)
            if frac_a is not None and frac_b is not None:
                # Coordinates were wrapped for output, so re-derive the bond
                # image shift from the wrapped positions (same convention as
                # CIFWriter._minimum_image_bond_geometry).
                wa = tuple(f - math.floor(f) for f in frac_a)
                wb = tuple(f - math.floor(f) for f in frac_b)
                shift, distance = _minimum_image_shift(cell, wa, wb)
                tokens[bond_columns.index("_geom_bond_site_symmetry_1")] = "."
                tokens[bond_columns.index("_geom_bond_site_symmetry_2")] = _format_p1_shift(shift)
                if "_geom_bond_distance" in bond_columns:
                    tokens[bond_columns.index("_geom_bond_distance")] = f"{distance:.3f}"
                line = " ".join(tokens) + "\n"
        out_lines.append(line)
        i += 1
    output_path.write_text("".join(out_lines))


def relax_cif_clashes(
    input_path: str | Path,
    output_path: str | Path | None = None,
    config: SoftRelaxConfig | None = None,
    *,
    clash_cutoff: float = 1.05,
) -> SoftRelaxReport:
    """Relieve non-bonded clashes in a cofkit-style explicit-bond P1 CIF.

    Runs the staged soft-repulsion / bond-spring descent and writes a new CIF
    with updated coordinates (cell, connectivity, and all non-coordinate CIF
    content preserved). Returns a report with before/after clash metrics.
    """
    config = config or SoftRelaxConfig()
    input_path = Path(input_path)
    if output_path is None:
        output_path = input_path.with_name(input_path.stem + "_softrelaxed.cif")
    output_path = Path(output_path)

    block = gemmi.cif.read_file(str(input_path)).sole_block()
    system, warnings = _parse_system(block, config)

    cart = [_frac_to_cart(system.orth, f) for f in system.frac]

    def _springs(bond_like, k_spring):
        indexed = [((i, j), _frac_to_cart(system.orth, s)) for i, j, s in bond_like]
        rest = _bond_distances(cart, indexed)
        return [((ij), sc, r0, k_spring) for (ij, sc), r0 in zip(indexed, rest)]

    springs = _springs(system.bonds, config.bond_force_constant)
    springs += _springs(system.angles13, config.urey_bradley_force_constant)
    bond_count = len(system.bonds)
    bond_rest = [s[2] for s in springs[:bond_count]]
    bonds_cart = [(s[0], s[1]) for s in springs[:bond_count]]

    min_before, clashes_before = _min_heavy_distance_and_clashes(system, clash_cutoff)
    pairs = _build_pair_list(system, config)

    stage_energies: list[float] = []
    converged = True
    for coefficient in config.repulsion_coefficients:
        energy, forces = _energy_and_forces(cart, springs, pairs, coefficient, config)
        step = config.initial_step
        stage_converged = False
        for iteration in range(config.max_steps_per_stage):
            if iteration > 0 and iteration % config.pair_list_interval == 0:
                system.frac = [
                    list(system.cell.fractionalize(gemmi.Position(*c))) for c in cart
                ]
                pairs = _build_pair_list(system, config)
                energy, forces = _energy_and_forces(cart, springs, pairs, coefficient, config)
            max_force = max(
                math.sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]) for f in forces
            ) if forces else 0.0
            if max_force < config.force_tolerance:
                stage_converged = True
                break
            proposal = [c[:] for c in cart]
            for idx, f in enumerate(forces):
                norm = math.sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2])
                if norm < 1e-12:
                    continue
                scale = step if norm * step <= config.max_displacement else config.max_displacement / norm
                proposal[idx][0] += f[0] * scale
                proposal[idx][1] += f[1] * scale
                proposal[idx][2] += f[2] * scale
            new_energy, new_forces = _energy_and_forces(
                proposal, springs, pairs, coefficient, config
            )
            if new_energy < energy:
                cart = proposal
                energy, forces = new_energy, new_forces
                step = min(step * 1.15, 0.5)
            else:
                step *= 0.5
                if step < 1e-5:
                    break
        stage_energies.append(energy)
        converged = converged and stage_converged

    new_frac = []
    for c in cart:
        frac = system.cell.fractionalize(gemmi.Position(*c))
        new_frac.append([frac.x, frac.y, frac.z])
    system.frac = [list(f) for f in new_frac]
    min_after, clashes_after = _min_heavy_distance_and_clashes(system, clash_cutoff)

    final_bonds = _bond_distances(cart, bonds_cart)
    max_bond_drift = max(
        (abs(a - b) for a, b in zip(final_bonds, bond_rest)), default=0.0
    )

    _rewrite_cif_coordinates(input_path, output_path, new_frac, system.labels, system.cell)

    return SoftRelaxReport(
        input_path=str(input_path),
        output_path=str(output_path),
        n_atoms=len(system.labels),
        n_bonds=len(system.bonds),
        n_pairs=len(pairs),
        n_angle_restraints=len(system.angles13),
        stage_energies=tuple(stage_energies),
        converged=converged,
        min_heavy_distance_before=min_before,
        min_heavy_distance_after=min_after,
        clashes_before=clashes_before,
        clashes_after=clashes_after,
        max_bond_drift=max_bond_drift,
        warnings=tuple(dict.fromkeys(warnings)),
    )
