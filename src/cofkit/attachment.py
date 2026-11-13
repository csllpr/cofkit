"""Versioned attachment metadata, independent of the compact COFid format.

Atom indices refer to the heavy atoms obtained by parsing the recipe's canonical
precursor SMILES. Instance ids and integer images name copies in the recorded
cell; they are references, not new material identifiers.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Mapping

from .geometry import Vec3, cross, dot
from .model import Image3

# Schema identity, not a numerical modeling parameter.
ATTACHMENT_SCHEMA_VERSION = 1
ATTACHMENT_SOURCE_ATOM_INDEX_PROP = "cofkit_event_original_atom_idx"


@dataclass(frozen=True)
class AttachmentPosition:
    atom_index: int
    position: Vec3
    source_label: str = ""


@dataclass(frozen=True)
class AttachmentInstance:
    id: str
    monomer_index: int
    seed_positions: tuple[AttachmentPosition, ...] = ()


@dataclass(frozen=True)
class AttachmentSite:
    instance_id: str
    atom_indices: tuple[int, ...]
    periodic_image: Image3 = (0, 0, 0)


@dataclass(frozen=True)
class AttachmentEvent:
    id: str
    template_id: str
    participants: tuple[AttachmentSite, ...]
    source_event_id: str = ""
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class AttachmentAtom:
    instance_id: str
    atom_index: int
    element: str

    @property
    def key(self) -> tuple[str, int]:
        return self.instance_id, self.atom_index


@dataclass(frozen=True)
class AttachmentBond:
    first: tuple[str, int]
    second: tuple[str, int]
    image: Image3 = (0, 0, 0)


@dataclass(frozen=True)
class AttachmentRecipe:
    cofid: str
    cell: tuple[Vec3, Vec3, Vec3]
    instances: tuple[AttachmentInstance, ...]
    events: tuple[AttachmentEvent, ...]
    reference_atoms: tuple[AttachmentAtom, ...] = ()
    reference_bonds: tuple[AttachmentBond, ...] = ()
    provenance: Mapping[str, object] = field(default_factory=dict)
    schema_version: int = ATTACHMENT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, object]:
        # Metadata can contain vectors represented by tuples or JSON lists.
        # Emit one JSON-native representation, including on an in-memory recipe.
        return json.loads(json.dumps(asdict(self), allow_nan=False))

    @property
    def fingerprint(self) -> str:
        """Recipe checksum, not a canonical periodic-structure identifier."""
        payload = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False
        )
        return hashlib.sha256(payload.encode()).hexdigest()

    def validate(self) -> None:
        from rdkit import Chem
        from .cofid import cofid_to_build_request
        from .reactions import ReactionLibrary

        _integer(self.schema_version, "schema version")
        if self.schema_version != ATTACHMENT_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported attachment schema version {self.schema_version!r}"
            )
        request = cofid_to_build_request(self.cofid)
        if len(self.cell) != 3 or any(
            len(v) != 3 or not all(math.isfinite(x) for x in v) for v in self.cell
        ):
            raise ValueError(
                "attachment cell must contain three finite Cartesian vectors"
            )
        if dot(self.cell[0], cross(self.cell[1], self.cell[2])) == 0:
            raise ValueError("attachment cell is singular")
        monomers = tuple(
            Chem.MolFromSmiles(m.canonical_smiles) for m in request.monomers
        )
        instances = {i.id: i for i in self.instances}
        if not instances or len(instances) != len(self.instances):
            raise ValueError("attachment instance ids must be nonempty and unique")
        for instance in self.instances:
            if not instance.id or any(c.isspace() for c in instance.id):
                raise ValueError(
                    "attachment instance ids must be nonempty and contain no whitespace"
                )
            _integer(instance.monomer_index, "monomer index")
            if not 0 <= instance.monomer_index < len(monomers):
                raise ValueError(
                    f"attachment instance {instance.id!r} has an unknown monomer index"
                )
            seen = set()
            for seed in instance.seed_positions:
                self._check_atom(instance.id, seed.atom_index, instances, monomers)
                if seed.atom_index in seen:
                    raise ValueError("duplicate attachment seed atom")
                seen.add(seed.atom_index)
                _vector(seed.position, "seed position")
        if not self.events or len({e.id for e in self.events}) != len(self.events):
            raise ValueError(
                "attachment event ids must be unique and events must not be empty"
            )
        template = ReactionLibrary.builtin().get(request.template_id)
        for event in self.events:
            if not event.id:
                raise ValueError("attachment event ids must not be empty")
            if event.template_id != request.template_id:
                raise ValueError(
                    "attachment event template disagrees with the COFid linkage"
                )
            if len(event.participants) != template.arity:
                raise ValueError(
                    f"attachment event {event.id!r} has the wrong participant count"
                )
            for site in event.participants:
                _image(site.periodic_image)
                if not site.atom_indices or len(set(site.atom_indices)) != len(
                    site.atom_indices
                ):
                    raise ValueError(
                        "attachment site atom indices must be nonempty and unique"
                    )
                for atom in site.atom_indices:
                    self._check_atom(site.instance_id, atom, instances, monomers)
            for key in (
                "ring_center_fractional",
                "ring_normal",
                "ring_center_offset_cartesian",
            ):
                if key in event.metadata:
                    vector = _vector(event.metadata[key], key)
                    if key == "ring_normal" and not any(vector):
                        raise ValueError("attachment ring normal must not be zero")
            if "attachment_retained_oxygen_ids" in event.metadata:
                selected = event.metadata["attachment_retained_oxygen_ids"]
                if request.linkage_code != "boroxine" or len(selected) != len(
                    event.participants
                ):
                    raise ValueError(
                        "attachment retained oxygen choices require one choice per boroxine participant"
                    )
                for site, index in zip(event.participants, selected):
                    mol = self._check_atom(site.instance_id, index, instances, monomers)
                    if mol.GetAtomWithIdx(index).GetSymbol() != "O":
                        raise ValueError(
                            "attachment retained oxygen choice must name oxygen"
                        )
        atoms = {a.key: a.element for a in self.reference_atoms}
        if len(atoms) != len(self.reference_atoms):
            raise ValueError("duplicate attachment reference atom")
        for atom in self.reference_atoms:
            mol = self._check_atom(
                atom.instance_id, atom.atom_index, instances, monomers
            )
            if (
                mol.GetAtomWithIdx(atom.atom_index).GetSymbol() != atom.element
                or atom.element == "H"
            ):
                raise ValueError(
                    "attachment reference atom element disagrees with its precursor"
                )
        seen_bonds = set()
        for bond in self.reference_bonds:
            _image(bond.image)
            if bond.first not in atoms or bond.second not in atoms:
                raise ValueError("attachment reference bond has an unknown endpoint")
            key = canonical_attachment_bond(bond.first, bond.second, bond.image)
            if key in seen_bonds:
                raise ValueError("duplicate attachment reference bond")
            seen_bonds.add(key)
        if self.reference_bonds and not self.reference_atoms:
            raise ValueError("attachment reference bonds require reference atoms")

    @staticmethod
    def _check_atom(instance_id, index, instances, monomers):
        _integer(index, "atom index")
        if instance_id not in instances:
            raise ValueError(f"unknown attachment instance {instance_id!r}")
        mol = monomers[instances[instance_id].monomer_index]
        if (
            index < 0
            or index >= mol.GetNumAtoms()
            or mol.GetAtomWithIdx(index).GetSymbol() == "H"
        ):
            raise ValueError(
                f"attachment atom index {index!r} is not a precursor heavy atom"
            )
        return mol

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> AttachmentRecipe:
        try:
            recipe = cls(
                cofid=str(value["cofid"]),
                cell=tuple(_vector(v, "cell vector") for v in value["cell"]),
                instances=tuple(
                    AttachmentInstance(
                        id=str(i["id"]),
                        monomer_index=_integer(i["monomer_index"], "monomer index"),
                        seed_positions=tuple(
                            AttachmentPosition(
                                _integer(p["atom_index"], "atom index"),
                                _vector(p["position"], "seed position"),
                                str(p.get("source_label", "")),
                            )
                            for p in i.get("seed_positions", ())
                        ),
                    )
                    for i in value["instances"]
                ),
                events=tuple(
                    AttachmentEvent(
                        id=str(e["id"]),
                        template_id=str(e["template_id"]),
                        participants=tuple(
                            AttachmentSite(
                                str(s["instance_id"]),
                                tuple(
                                    _integer(a, "atom index") for a in s["atom_indices"]
                                ),
                                _image(s.get("periodic_image", (0, 0, 0))),
                            )
                            for s in e["participants"]
                        ),
                        source_event_id=str(e.get("source_event_id", "")),
                        metadata=dict(e.get("metadata", {})),
                    )
                    for e in value["events"]
                ),
                reference_atoms=tuple(
                    AttachmentAtom(
                        str(a["instance_id"]),
                        _integer(a["atom_index"], "atom index"),
                        str(a["element"]),
                    )
                    for a in value.get("reference_atoms", ())
                ),
                reference_bonds=tuple(
                    AttachmentBond(
                        _atom_ref(b["first"]),
                        _atom_ref(b["second"]),
                        _image(b.get("image", (0, 0, 0))),
                    )
                    for b in value.get("reference_bonds", ())
                ),
                provenance=dict(value.get("provenance", {})),
                schema_version=_integer(value["schema_version"], "schema version"),
            )
            recipe.validate()
            return recipe
        except (KeyError, TypeError, OverflowError) as exc:
            raise ValueError(f"invalid attachment recipe: {exc}") from exc


def _integer(value, label):
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"attachment {label} must be an integer")
    return value


def _vector(value, label):
    if isinstance(value, (str, bytes)) or len(value) != 3:
        raise ValueError(f"attachment {label} must have three components")
    result = tuple(float(v) for v in value)
    if not all(math.isfinite(v) for v in result):
        raise ValueError(f"attachment {label} must be finite")
    return result


def _image(value):
    if len(value) != 3:
        raise ValueError("attachment periodic images must have three components")
    return tuple(_integer(v, "periodic image component") for v in value)


def _atom_ref(value):
    if len(value) != 2:
        raise ValueError(
            "attachment atom references require an instance id and atom index"
        )
    return str(value[0]), _integer(value[1], "atom index")


AttachmentInput = AttachmentRecipe | Mapping[str, object] | str | Path


def load_attachment_recipe(value: AttachmentInput) -> AttachmentRecipe:
    """Read a recipe, an attachment report, or full decomposition JSON."""
    if isinstance(value, AttachmentRecipe):
        value.validate()
        return value
    if isinstance(value, (str, Path)):
        value = json.loads(Path(value).read_text())
    if not isinstance(value, Mapping):
        raise ValueError("attachment recipe must be a JSON object")
    if "metadata" in value and "schema_version" not in value:
        metadata = value["metadata"]
        if not isinstance(metadata, Mapping):
            raise ValueError("decomposition metadata must be an object")
        value = metadata.get("attachment", {})
    if not isinstance(value, Mapping):
        raise ValueError("attachment metadata must be an object")
    if "recipe" in value:
        if value.get("status") != "available" or value["recipe"] is None:
            raise ValueError(
                "attachment metadata does not contain one available recipe"
            )
        value = value["recipe"]
    return AttachmentRecipe.from_dict(value)


def canonical_attachment_bond(first, second, image):
    image = tuple(image)
    if second < first or (second == first and tuple(-v for v in image) < image):
        first, second, image = second, first, tuple(-v for v in image)
    return first, second, image


def compare_attachment_graphs(
    expected_atoms: Mapping[tuple[str, int], str],
    expected_bonds: Iterable[AttachmentBond],
    actual_atoms: Mapping[tuple[str, int], str],
    actual_bonds: Iterable[AttachmentBond],
) -> dict[str, object]:
    """Verify a known atom correspondence, allowing per-atom image rebasing.

    This is an exact labelled periodic-graph check, not a general graph
    recognizer. For parallel edges the finite gain multiset determines a unique
    translation offset, so no isomorphism search or numerical tolerance is used.
    """
    report = {
        "status": "mismatch",
        "scope": "periodic_heavy_graph",
        "bond_orders_checked": False,
    }
    if dict(expected_atoms) != dict(actual_atoms):
        return {**report, "reason": "retained heavy atoms or their elements differ"}

    def pairs(bonds):
        result = defaultdict(Counter)
        for b in bonds:
            a, z, g = canonical_attachment_bond(b.first, b.second, b.image)
            result[a, z][g] += 1
        return result

    expected, actual = pairs(expected_bonds), pairs(actual_bonds)
    if expected.keys() != actual.keys():
        return {**report, "reason": "heavy-atom bond endpoints differ"}
    adjacency = defaultdict(list)
    for (a, b), gains in expected.items():
        observed = actual[a, b]
        delta = tuple(y - x for x, y in zip(min(gains), min(observed)))
        shifted = Counter(
            {tuple(x + d for x, d in zip(g, delta)): n for g, n in gains.items()}
        )
        if shifted != observed or (a == b and any(delta)):
            return {
                **report,
                "reason": "parallel or self-image bond translations differ",
            }
        adjacency[a].append((b, delta))
        adjacency[b].append((a, tuple(-v for v in delta)))
    potentials = {}
    for root in expected_atoms:
        if root in potentials:
            continue
        potentials[root] = (0, 0, 0)
        queue = deque([root])
        while queue:
            a = queue.popleft()
            for b, delta in adjacency[a]:
                required = tuple(x + d for x, d in zip(potentials[a], delta))
                if b in potentials:
                    if potentials[b] != required:
                        return {
                            **report,
                            "reason": "periodic cycle translations differ",
                        }
                else:
                    potentials[b] = required
                    queue.append(b)
    return {
        **report,
        "status": "verified",
        "reason": "all heavy atoms and periodic edges match under integer atom-image rebasing",
        "n_atoms": len(expected_atoms),
        "n_bonds": sum(sum(v.values()) for v in expected.values()),
    }
