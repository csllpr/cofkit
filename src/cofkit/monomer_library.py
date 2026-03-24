from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping

from .batch_models import BatchMonomerRecord
from .chem.rdkit import RDKitMotifBuilder
from .model import MonomerSpec, ReactionTemplate
from .reactions import BinaryBridgeRole, ReactionLibrary


_AUTO_DETECT_GENERIC_SUPPRESSION: Mapping[str, tuple[str, ...]] = {
    "aldehyde": ("keto_aldehyde",),
}


@dataclass(frozen=True)
class MonomerRoleResolver:
    motif_builder: RDKitMotifBuilder
    generic_suppression: Mapping[str, tuple[str, ...]] = field(
        default_factory=lambda: dict(_AUTO_DETECT_GENERIC_SUPPRESSION)
    )

    @classmethod
    def builtin(cls) -> "MonomerRoleResolver":
        return cls(motif_builder=RDKitMotifBuilder.builtin())

    def supported_motif_kinds(self) -> tuple[str, ...]:
        return self.motif_builder.supported_motif_kinds()

    def auto_detect_monomer_candidates(
        self,
        smiles: str,
        *,
        allowed_motif_kinds: Iterable[str] | None = None,
        num_conformers: int = 2,
        random_seed: int = 0xC0F,
    ) -> dict[str, MonomerSpec]:
        if allowed_motif_kinds is None:
            motif_kinds = self.motif_builder.supported_motif_kinds()
        else:
            allowed = set(allowed_motif_kinds)
            motif_kinds = tuple(
                motif_kind
                for motif_kind in self.motif_builder.supported_motif_kinds()
                if motif_kind in allowed
            )
        candidates: dict[str, MonomerSpec] = {}
        for motif_kind in motif_kinds:
            try:
                monomer = self.motif_builder.build_monomer(
                    "__auto_detect__",
                    "__auto_detect__",
                    smiles,
                    motif_kind,
                    num_conformers=max(1, min(2, num_conformers)),
                    random_seed=random_seed,
                )
            except Exception:
                continue
            candidates[motif_kind] = monomer
        return candidates

    def resolve_detected_motif_kind(
        self,
        candidates: Mapping[str, MonomerSpec],
    ) -> tuple[str, MonomerSpec]:
        if not candidates:
            allowed_text = tuple(self.motif_builder.supported_motif_kinds())
            raise ValueError(
                f"could not auto-detect a supported motif kind from SMILES against {allowed_text!r}"
            )
        remaining = dict(candidates)
        for generic_kind, specific_kinds in self.generic_suppression.items():
            if generic_kind not in remaining:
                continue
            if any(specific_kind in remaining for specific_kind in specific_kinds):
                remaining.pop(generic_kind, None)
        if len(remaining) != 1:
            details = tuple(
                f"{motif_kind}:{len(monomer.motifs)}"
                for motif_kind, monomer in sorted(remaining.items())
            )
            raise ValueError(
                f"auto-detected ambiguous motif kinds {details!r}; "
                "current batch generation requires one role per monomer"
            )
        motif_kind = next(iter(remaining))
        return motif_kind, remaining[motif_kind]

    def infer_record(
        self,
        smiles: str,
        *,
        record_id: str,
        name: str | None = None,
        source_path: str = "",
        source_line: int = 0,
        allowed_motif_kinds: Iterable[str] | None = None,
        library_stem: str = "",
        num_conformers: int = 2,
        random_seed: int = 0xC0F,
    ) -> BatchMonomerRecord:
        candidates = self.auto_detect_monomer_candidates(
            smiles,
            allowed_motif_kinds=allowed_motif_kinds,
            num_conformers=num_conformers,
            random_seed=random_seed,
        )
        motif_kind, monomer = self.resolve_detected_motif_kind(candidates)
        metadata: dict[str, object] = {
            "auto_detected": True,
            "detected_motif_kinds": tuple(sorted(candidates)),
            "detected_connectivities": {
                kind: len(candidate.motifs)
                for kind, candidate in sorted(candidates.items())
            },
        }
        if library_stem:
            metadata["library_stem"] = library_stem
        return BatchMonomerRecord(
            id=record_id,
            name=name or record_id,
            smiles=smiles,
            motif_kind=motif_kind,
            expected_connectivity=len(monomer.motifs),
            source_path=source_path,
            source_line=source_line,
            metadata=metadata,
        )


class BinaryBridgeLibraryLoader:
    def __init__(
        self,
        *,
        reaction_library: ReactionLibrary | None = None,
        role_resolver: MonomerRoleResolver | None = None,
    ) -> None:
        self.reaction_library = reaction_library or ReactionLibrary.builtin()
        self.role_resolver = role_resolver or MonomerRoleResolver.builtin()

    def selected_binary_bridge_template(
        self,
        *,
        allowed_reactions: Iterable[str],
        template_id: str | None = None,
    ) -> ReactionTemplate:
        allowed = tuple(allowed_reactions)
        if template_id is not None and template_id not in allowed:
            raise ValueError(
                f"requested template {template_id!r} is not enabled in allowed_reactions={allowed!r}"
            )
        requested_ids = (template_id,) if template_id is not None else allowed
        binary_bridge_templates = tuple(
            self.reaction_library.get(candidate_id)
            for candidate_id in requested_ids
            if self.reaction_library.supports_binary_bridge_pair_generation(candidate_id)
        )
        if len(binary_bridge_templates) != 1:
            raise ValueError(
                "batch pair generation currently requires exactly one supported binary bridge template; "
                f"got {tuple(template.id for template in binary_bridge_templates)!r} from {requested_ids!r}"
            )
        return binary_bridge_templates[0]

    def binary_bridge_roles(self, template: ReactionTemplate) -> tuple[BinaryBridgeRole, ...]:
        profile = self.reaction_library.linkage_profile(template)
        if profile is None or len(profile.binary_bridge_roles) != 2:
            raise ValueError(f"template {template.id!r} does not expose a two-role binary bridge profile")
        return profile.binary_bridge_roles

    def library_prefix_for_motif_kind(self, motif_kind: str) -> str:
        for profile in self.reaction_library.linkage_profiles.values():
            for role in profile.binary_bridge_roles:
                if role.motif_kind == motif_kind:
                    return role.library_prefix or f"{role.motif_kind}s"
        return f"{motif_kind}s"

    def load_smiles_library(
        self,
        path: str | Path,
        *,
        motif_kind: str,
        expected_connectivity: int,
        id_prefix: str | None = None,
    ) -> tuple[BatchMonomerRecord, ...]:
        library_path = Path(path)
        prefix = id_prefix or library_path.stem
        rows: list[BatchMonomerRecord] = []
        with library_path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line_number == 1 and line.lower() == "smiles":
                    continue
                index = len(rows) + 1
                rows.append(
                    BatchMonomerRecord(
                        id=f"{prefix}_{index:04d}",
                        name=f"{prefix}_{index:04d}",
                        smiles=line,
                        motif_kind=motif_kind,
                        expected_connectivity=expected_connectivity,
                        source_path=str(library_path),
                        source_line=line_number,
                        metadata={"library_stem": library_path.stem},
                    )
                )
        return tuple(rows)

    def load_smiles_library_auto(
        self,
        path: str | Path,
        *,
        allowed_motif_kinds: Iterable[str] | None = None,
        id_prefix: str | None = None,
        num_conformers: int = 2,
        random_seed: int = 0xC0F,
    ) -> tuple[BatchMonomerRecord, ...]:
        library_path = Path(path)
        prefix = id_prefix or library_path.stem
        rows: list[BatchMonomerRecord] = []
        with library_path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line_number == 1 and line.lower() == "smiles":
                    continue
                index = len(rows) + 1
                rows.append(
                    self.role_resolver.infer_record(
                        line,
                        record_id=f"{prefix}_{index:04d}",
                        name=f"{prefix}_{index:04d}",
                        source_path=str(library_path),
                        source_line=line_number,
                        allowed_motif_kinds=allowed_motif_kinds,
                        library_stem=library_path.stem,
                        num_conformers=num_conformers,
                        random_seed=random_seed,
                    )
                )
        return tuple(rows)

    def load_binary_bridge_test_set(
        self,
        root: str | Path,
        *,
        allowed_reactions: Iterable[str],
        template_id: str | None = None,
        auto_detect: bool = False,
        num_conformers: int = 2,
        random_seed: int = 0xC0F,
    ) -> dict[str, tuple[BatchMonomerRecord, ...]]:
        if auto_detect:
            return self.load_binary_bridge_test_set_auto(
                root,
                allowed_reactions=allowed_reactions,
                template_id=template_id,
                num_conformers=num_conformers,
                random_seed=random_seed,
            )
        base = Path(root)
        template = self.selected_binary_bridge_template(
            allowed_reactions=allowed_reactions,
            template_id=template_id,
        )
        libraries: dict[str, tuple[BatchMonomerRecord, ...]] = {}
        for role in self.binary_bridge_roles(template):
            prefix = role.library_prefix or f"{role.motif_kind}s"
            for path in sorted(
                base.glob(f"{prefix}_count_*.txt"),
                key=lambda candidate: int(candidate.stem.rsplit("_", 1)[1]),
            ):
                connectivity = int(path.stem.rsplit("_", 1)[1])
                libraries[path.stem] = self.load_smiles_library(
                    path,
                    motif_kind=role.motif_kind,
                    expected_connectivity=connectivity,
                )
        return libraries

    def load_binary_bridge_test_set_auto(
        self,
        root: str | Path,
        *,
        allowed_reactions: Iterable[str],
        template_id: str | None = None,
        num_conformers: int = 2,
        random_seed: int = 0xC0F,
    ) -> dict[str, tuple[BatchMonomerRecord, ...]]:
        base = Path(root)
        template = self.selected_binary_bridge_template(
            allowed_reactions=allowed_reactions,
            template_id=template_id,
        )
        allowed_motif_kinds = tuple(role.motif_kind for role in self.binary_bridge_roles(template))
        grouped_records: dict[str, list[BatchMonomerRecord]] = {}
        for path in sorted(base.glob("*.txt")):
            for record in self.load_smiles_library_auto(
                path,
                allowed_motif_kinds=allowed_motif_kinds,
                num_conformers=num_conformers,
                random_seed=random_seed,
            ):
                prefix = self.library_prefix_for_motif_kind(record.motif_kind)
                key = f"{prefix}_count_{record.expected_connectivity}"
                grouped_records.setdefault(key, []).append(record)
        return {key: tuple(records) for key, records in sorted(grouped_records.items())}

    def available_binary_bridge_template_ids(
        self,
        root: str | Path,
        *,
        allowed_reactions: Iterable[str] | None = None,
        auto_detect_libraries: bool = False,
        num_conformers: int = 2,
        random_seed: int = 0xC0F,
    ) -> tuple[str, ...]:
        base = Path(root)
        allowed = set(allowed_reactions or self.reaction_library.templates)
        available: list[str] = []
        for template_id, profile in sorted(self.reaction_library.linkage_profiles.items()):
            if template_id not in allowed or not profile.supports_binary_bridge_pair_generation:
                continue
            roles = tuple(profile.binary_bridge_roles)
            if auto_detect_libraries:
                allowed_motif_kinds = tuple(role.motif_kind for role in roles)
                grouped_records: dict[str, list[BatchMonomerRecord]] = {}
                for path in sorted(base.glob("*.txt")):
                    for record in self.load_smiles_library_auto(
                        path,
                        allowed_motif_kinds=allowed_motif_kinds,
                        num_conformers=num_conformers,
                        random_seed=random_seed,
                    ):
                        prefix = self.library_prefix_for_motif_kind(record.motif_kind)
                        key = f"{prefix}_count_{record.expected_connectivity}"
                        grouped_records.setdefault(key, []).append(record)
                libraries = {key: tuple(records) for key, records in grouped_records.items()}
            else:
                libraries: dict[str, tuple[BatchMonomerRecord, ...]] = {}
                for role in roles:
                    prefix = role.library_prefix or f"{role.motif_kind}s"
                    for path in sorted(
                        base.glob(f"{prefix}_count_*.txt"),
                        key=lambda candidate: int(candidate.stem.rsplit("_", 1)[1]),
                    ):
                        connectivity = int(path.stem.rsplit("_", 1)[1])
                        libraries[path.stem] = self.load_smiles_library(
                            path,
                            motif_kind=role.motif_kind,
                            expected_connectivity=connectivity,
                        )
            first_role, second_role = roles
            first_by_connectivity = self._records_by_connectivity(libraries, prefix=self.library_prefix_for_motif_kind(first_role.motif_kind))
            second_by_connectivity = self._records_by_connectivity(libraries, prefix=self.library_prefix_for_motif_kind(second_role.motif_kind))
            if self._has_candidate_pairs(first_by_connectivity, second_by_connectivity):
                available.append(template_id)
        return tuple(available)

    def _has_candidate_pairs(
        self,
        first_by_connectivity: Mapping[int, tuple[BatchMonomerRecord, ...]],
        second_by_connectivity: Mapping[int, tuple[BatchMonomerRecord, ...]],
    ) -> bool:
        if any(connectivity >= 3 and connectivity in second_by_connectivity for connectivity in first_by_connectivity):
            return True
        if 2 in second_by_connectivity and any(connectivity >= 3 for connectivity in first_by_connectivity):
            return True
        if 2 in first_by_connectivity and any(connectivity >= 3 for connectivity in second_by_connectivity):
            return True
        return False

    def _records_by_connectivity(
        self,
        libraries: Mapping[str, tuple[BatchMonomerRecord, ...]],
        *,
        prefix: str,
    ) -> dict[int, tuple[BatchMonomerRecord, ...]]:
        records_by_connectivity: dict[int, tuple[BatchMonomerRecord, ...]] = {}
        for key, records in libraries.items():
            if not key.startswith(prefix):
                continue
            connectivity = int(key.rsplit("_", 1)[1])
            records_by_connectivity[connectivity] = records
        return records_by_connectivity
