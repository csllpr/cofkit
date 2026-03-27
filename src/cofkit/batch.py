from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from itertools import permutations, product
from math import atan2, cos, pi, sin
from pathlib import Path
from threading import Lock
from typing import Iterable, Mapping

from .batch_models import BatchMonomerRecord, BatchPairSummary, BatchRunSummary, BuiltBatchMonomer
from .chem.rdkit import build_rdkit_monomer
from .cif import CIFWriter
from .embedding import EmbeddingConfig, EmbeddingResult, PeriodicEmbedder
from .engine import COFEngineConfig, COFProject
from .geometry import (
    Mat3,
    Vec3,
    Frame,
    add,
    cross,
    dot,
    matmul_vec,
    norm,
    normalize,
    rotation_from_frame_to_axes,
    scale,
    sub,
)
from .model import AssemblyState, Candidate, MonomerInstance, MonomerSpec, MotifRef, Pose, ReactionEvent, ReactionTemplate
from .monomer_library import BinaryBridgeLibraryLoader, MonomerRoleResolver
from .optimizer import ContinuousOptimizer, OptimizerConfig
from .planner import AssignmentPlan, NetPlan, NetPlanner
from .product_graph import PeriodicProductGraph
from .post_build_conversions import annotate_post_build_conversions
from .reactions import (
    BinaryBridgeRole,
    ReactionLibrary,
    bridge_target_distance,
)
from .indexed_topology_layouts import ExpandedIndexedNodeSite, ExpandedIndexedTopology, expand_indexed_topology
from .linkage_geometry import effective_motif_origin
from .scoring import CandidateScorer
from .search import AssignmentOutcome, AssignmentSolver
from .single_node_topologies import (
    ExpandedSingleNodeEdge,
    ExpandedSingleNodeNodeSite,
    ExpandedSingleNodeTopology,
    expand_single_node_topology,
    list_supported_single_node_topology_ids,
    resolve_single_node_topology_layout,
)
from .single_node_topologies_3d import (
    expand_three_d_single_node_topology,
    list_supported_3d_single_node_topology_ids,
    resolve_three_d_single_node_topology_layout,
)
from .topology_builders import PairTopologyBuildRequest, PairTopologyBuilder, PairTopologyBuilderRegistry
from .topologies import default_topology_repository, get_topology_hint
from .validation import CoarseStructureValidator, CoarseValidationReport, CoarseValidationThresholds

_PROCESS_BATCH_GENERATOR = None


_CURATED_DEFAULT_TOPOLOGY_IDS: tuple[str, ...] = (
    "hcb",
    "sql",
    "kgm",
    "hxl",
    "fes",
    "dia",
    "pts",
    "ctn",
    "bor",
    "srs",
    "kgd",
    "tth",
    "bex",
    "lon",
    "pcu",
    "acs",
    "qtz",
    "flu",
    "bcu",
    "tbo",
    "fof",
    "rht",
)


@dataclass(frozen=True)
class BatchGenerationConfig:
    allowed_reactions: tuple[str, ...] = ("imine_bridge",)
    target_dimensionality: str = "2D"
    hcb_topology_id: str = "hcb"
    topology_ids: tuple[str, ...] = ()
    single_node_topology_ids: tuple[str, ...] = ()
    use_indexed_topology_defaults: bool = False
    rdkit_num_conformers: int = 8
    rdkit_random_seed: int = 0xC0F
    retain_top_results: int = 25
    enumerate_all_topologies: bool = True
    post_build_conversions: tuple[str, ...] = ()
    write_cif: bool = True
    max_cif_exports: int | None = None
    max_workers: int = 8
    hard_hard_max_bridge_distance: float = 2.5
    separate_cif_outputs_by_validation: bool = True
    validation_thresholds: CoarseValidationThresholds = field(default_factory=CoarseValidationThresholds)
    embedding_config: EmbeddingConfig = field(default_factory=EmbeddingConfig)
    engine_config: COFEngineConfig = field(default_factory=COFEngineConfig)
    optimizer_config: OptimizerConfig = field(default_factory=OptimizerConfig)


@dataclass(frozen=True)
class _LinkerLayout:
    rotation: Mat3
    motif_ids: tuple[str, str]
    local_origins: tuple[Vec3, Vec3]
    offset_a: float
    offset_b: float
    span: float


@dataclass(frozen=True)
class _TopologyEvaluation:
    pair_mode: str
    candidates: tuple[Candidate, ...]
    available_topologies: tuple[str, ...]
    failed_topologies: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class _BinaryBridgePairContext:
    template: ReactionTemplate
    first: MonomerSpec
    second: MonomerSpec
    first_id: str
    second_id: str
    pair_id: str
    role_ids: tuple[str, str]


@dataclass(frozen=True)
class _NodeSitePlacement:
    rotation: Mat3
    motif_ids: tuple[str, ...]
    offsets: tuple[float, ...]
    motif_id_by_edge: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class _BatchPairTask:
    first_record: BatchMonomerRecord
    second_record: BatchMonomerRecord
    out_dir: str | Path | None
    write_cif: bool | None


def _init_batch_process_worker(
    config: BatchGenerationConfig,
    reaction_library: ReactionLibrary,
) -> None:
    global _PROCESS_BATCH_GENERATOR
    _PROCESS_BATCH_GENERATOR = BatchStructureGenerator(
        config=config,
        reaction_library=reaction_library,
    )


def _run_batch_pair_task_in_process(
    task: _BatchPairTask,
) -> tuple[tuple[BatchPairSummary, ...], int]:
    if _PROCESS_BATCH_GENERATOR is None:
        raise RuntimeError("batch process worker was not initialized")
    return _PROCESS_BATCH_GENERATOR._run_batch_pair_task(task)


class BatchStructureGenerator:
    def __init__(
        self,
        config: BatchGenerationConfig | None = None,
        reaction_library: ReactionLibrary | None = None,
        smiles_monomer_builder=None,
    ):
        self.config = config or BatchGenerationConfig()
        self.reaction_library = reaction_library or ReactionLibrary.builtin()
        self.smiles_monomer_builder = smiles_monomer_builder or build_rdkit_monomer
        self._role_resolver = MonomerRoleResolver.builtin()
        self._library_loader = BinaryBridgeLibraryLoader(
            reaction_library=self.reaction_library,
            role_resolver=self._role_resolver,
        )
        self.net_planner = NetPlanner()
        self.assignment_solver = AssignmentSolver()
        self.embedder = PeriodicEmbedder(self.config.embedding_config)
        self.scorer = CandidateScorer()
        self.optimizer = ContinuousOptimizer(self.config.optimizer_config, scorer=self.scorer)
        self.cif_writer = CIFWriter()
        self.structure_validator = CoarseStructureValidator(thresholds=self.config.validation_thresholds)
        self._monomer_cache: dict[str, BuiltBatchMonomer] = {}
        self._monomer_cache_lock = Lock()
        self._spatial_rotation_cache: dict[tuple[str, tuple[tuple[float, float, float], ...]], tuple[Mat3, tuple[str, ...], tuple[float, ...]]] = {}
        self._spatial_rotation_cache_lock = Lock()
        self._topology_id_cache: dict[tuple[tuple[int, int], str], tuple[str, ...]] = {}
        self._topology_id_cache_lock = Lock()
        self._topology_repository = default_topology_repository()
        self._pair_topology_builders = PairTopologyBuilderRegistry(
            (
                PairTopologyBuilder(
                    id="two-node-node-node",
                    description="legacy two-node 2D single-node node-node builder",
                    supports=self._supports_two_node_node_node_builder,
                    build=self._build_two_node_node_node_request,
                ),
                PairTopologyBuilder(
                    id="two-node-node-linker",
                    description="legacy two-node 2D single-node node-linker builder",
                    supports=self._supports_two_node_node_linker_builder,
                    build=self._build_two_node_node_linker_request,
                ),
                PairTopologyBuilder(
                    id="expanded-node-node",
                    description="expanded single-node node-node builder",
                    supports=self._supports_expanded_node_node_builder,
                    build=self._build_expanded_node_node_request,
                ),
                PairTopologyBuilder(
                    id="expanded-node-linker",
                    description="expanded single-node node-linker builder",
                    supports=self._supports_expanded_node_linker_builder,
                    build=self._build_expanded_node_linker_request,
                ),
                PairTopologyBuilder(
                    id="indexed-node-node",
                    description="generic indexed topology node-node builder",
                    supports=self._supports_indexed_node_node_builder,
                    build=self._build_indexed_node_node_request,
                ),
                PairTopologyBuilder(
                    id="indexed-node-linker",
                    description="generic indexed topology node-linker builder",
                    supports=self._supports_indexed_node_linker_builder,
                    build=self._build_indexed_node_linker_request,
                ),
            )
        )

    def load_smiles_library(
        self,
        path: str | Path,
        *,
        motif_kind: str,
        expected_connectivity: int,
        id_prefix: str | None = None,
    ) -> tuple[BatchMonomerRecord, ...]:
        return self._library_loader.load_smiles_library(
            path,
            motif_kind=motif_kind,
            expected_connectivity=expected_connectivity,
            id_prefix=id_prefix,
        )

    def load_smiles_library_auto(
        self,
        path: str | Path,
        *,
        allowed_motif_kinds: Iterable[str] | None = None,
        id_prefix: str | None = None,
    ) -> tuple[BatchMonomerRecord, ...]:
        return self._library_loader.load_smiles_library_auto(
            path,
            allowed_motif_kinds=allowed_motif_kinds,
            id_prefix=id_prefix,
            num_conformers=self._autodetect_num_conformers(),
            random_seed=self.config.rdkit_random_seed,
        )

    def infer_monomer_record(
        self,
        smiles: str,
        *,
        record_id: str,
        name: str | None = None,
        source_path: str = "",
        source_line: int = 0,
        allowed_motif_kinds: Iterable[str] | None = None,
        library_stem: str = "",
    ) -> BatchMonomerRecord:
        return self._library_loader.role_resolver.infer_record(
            smiles,
            record_id=record_id,
            name=name,
            source_path=source_path,
            source_line=source_line,
            allowed_motif_kinds=allowed_motif_kinds,
            library_stem=library_stem,
            num_conformers=self._autodetect_num_conformers(),
            random_seed=self.config.rdkit_random_seed,
        )

    def _auto_detect_monomer_candidates(
        self,
        smiles: str,
        *,
        allowed_motif_kinds: Iterable[str] | None = None,
    ) -> dict[str, MonomerSpec]:
        return self._role_resolver.auto_detect_monomer_candidates(
            smiles,
            allowed_motif_kinds=allowed_motif_kinds,
            num_conformers=self._autodetect_num_conformers(),
            random_seed=self.config.rdkit_random_seed,
        )

    def _resolve_detected_motif_kind(
        self,
        candidates: Mapping[str, MonomerSpec],
    ) -> tuple[str, MonomerSpec]:
        return self._role_resolver.resolve_detected_motif_kind(candidates)

    def load_binary_bridge_test_set(
        self,
        root: str | Path,
        *,
        template_id: str | None = None,
        auto_detect: bool = False,
    ) -> dict[str, tuple[BatchMonomerRecord, ...]]:
        return self._library_loader.load_binary_bridge_test_set(
            root,
            allowed_reactions=self.config.allowed_reactions,
            template_id=template_id,
            auto_detect=auto_detect,
            num_conformers=self._autodetect_num_conformers(),
            random_seed=self.config.rdkit_random_seed,
        )

    def load_binary_bridge_test_set_auto(
        self,
        root: str | Path,
        *,
        template_id: str | None = None,
    ) -> dict[str, tuple[BatchMonomerRecord, ...]]:
        return self._library_loader.load_binary_bridge_test_set_auto(
            root,
            allowed_reactions=self.config.allowed_reactions,
            template_id=template_id,
            num_conformers=self._autodetect_num_conformers(),
            random_seed=self.config.rdkit_random_seed,
        )

    def load_imine_test_set(self, root: str | Path) -> dict[str, tuple[BatchMonomerRecord, ...]]:
        return self.load_binary_bridge_test_set(root, template_id="imine_bridge")

    def available_binary_bridge_template_ids(
        self,
        root: str | Path,
        *,
        auto_detect_libraries: bool = False,
    ) -> tuple[str, ...]:
        available: list[str] = []
        for template_id, profile in sorted(self.reaction_library.linkage_profiles.items()):
            if not profile.supports_binary_bridge_pair_generation:
                continue
            libraries = self._library_loader.load_binary_bridge_test_set(
                root,
                allowed_reactions=(template_id,),
                template_id=template_id,
                auto_detect=auto_detect_libraries,
                num_conformers=self._autodetect_num_conformers(),
                random_seed=self.config.rdkit_random_seed,
            )
            template = self.reaction_library.get(template_id)
            pair_iterables = self._batch_pair_iterables(libraries, template=template)
            if any(next(iter(iterable), None) is not None for iterable in pair_iterables):
                available.append(template_id)
        return tuple(available)

    def _library_prefix_for_motif_kind(self, motif_kind: str) -> str:
        return self._library_loader.library_prefix_for_motif_kind(motif_kind)

    def build_monomer(self, record: BatchMonomerRecord) -> BuiltBatchMonomer:
        with self._monomer_cache_lock:
            cached = self._monomer_cache.get(record.id)
        if cached is not None:
            return cached

        try:
            monomer = self.smiles_monomer_builder(
                record.id,
                record.name,
                record.smiles,
                record.motif_kind,
                num_conformers=self.config.rdkit_num_conformers,
                random_seed=self.config.rdkit_random_seed,
            )
            actual_connectivity = len(monomer.motifs)
            if actual_connectivity != record.expected_connectivity:
                raise ValueError(
                    f"expected {record.expected_connectivity} motifs from {record.motif_kind!r} detection, "
                    f"got {actual_connectivity}"
                )
            result = BuiltBatchMonomer(record=record, monomer=monomer)
        except Exception as exc:  # pragma: no cover - exercised against large fixture libraries
            result = BuiltBatchMonomer(record=record, error=f"{type(exc).__name__}: {exc}")

        with self._monomer_cache_lock:
            existing = self._monomer_cache.get(record.id)
            if existing is not None:
                return existing
            self._monomer_cache[record.id] = result
        return result

    def _autodetect_num_conformers(self) -> int:
        return max(1, min(2, self.config.rdkit_num_conformers))

    def _selected_binary_bridge_template(self, *, template_id: str | None = None) -> ReactionTemplate:
        return self._library_loader.selected_binary_bridge_template(
            allowed_reactions=self.config.allowed_reactions,
            template_id=template_id,
        )

    def _binary_bridge_roles(self, template: ReactionTemplate) -> tuple[BinaryBridgeRole, ...]:
        return self._library_loader.binary_bridge_roles(template)

    def _resolve_binary_bridge_pair_context(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
        *,
        first_id: str,
        second_id: str,
        pair_id: str | None = None,
    ) -> _BinaryBridgePairContext:
        template = self._selected_binary_bridge_template()
        order = self.reaction_library.resolve_binary_bridge_pair_order(
            template,
            self._single_motif_kind(first),
            self._single_motif_kind(second),
        )
        monomers = (first, second)
        record_ids = (first_id, second_id)
        ordered_first, ordered_second = (monomers[index] for index in order.ordered_indices)
        ordered_first_id, ordered_second_id = (record_ids[index] for index in order.ordered_indices)
        resolved_pair_id = pair_id or f"{ordered_first_id}__{ordered_second_id}"
        return _BinaryBridgePairContext(
            template=template,
            first=ordered_first,
            second=ordered_second,
            first_id=ordered_first_id,
            second_id=ordered_second_id,
            pair_id=resolved_pair_id,
            role_ids=tuple(role.role_id for role in order.roles),
        )

    def generate_pair_candidate(
        self,
        amine_record: BatchMonomerRecord,
        aldehyde_record: BatchMonomerRecord,
        *,
        out_dir: str | Path | None = None,
        write_cif: bool | None = None,
        cif_export_index: int = 0,
    ) -> tuple[BatchPairSummary, Candidate | None]:
        summaries, candidates, _attempted_structures = self.generate_pair_candidates(
            amine_record,
            aldehyde_record,
            out_dir=None,
            write_cif=False,
            cif_export_start_index=0,
        )
        if not summaries:
            raise ValueError("pair evaluation produced no summaries")

        if candidates:
            best_candidate = max(candidates, key=lambda candidate: candidate.score)
            best_topology = best_candidate.metadata["net_plan"]["topology"]
            best_index = next(
                index
                for index, candidate in enumerate(candidates)
                if candidate.metadata["net_plan"]["topology"] == best_topology
            )
            summary = summaries[best_index]
            effective_write_cif = self.config.write_cif if write_cif is None else write_cif
            cif_path = None
            validation_metadata = None
            export_blocked = bool(summary.metadata.get("hard_hard_invalid_reasons"))
            if effective_write_cif and out_dir is not None and not export_blocked:
                amine = self.build_monomer(amine_record)
                aldehyde = self.build_monomer(aldehyde_record)
                assert amine.monomer is not None
                assert aldehyde.monomer is not None
                cif_path, validation_metadata = self._write_classified_candidate_cif(
                    out_dir=out_dir,
                    structure_id=summary.structure_id,
                    candidate=best_candidate,
                    monomer_specs=(amine.monomer, aldehyde.monomer),
                    provisional_summary=summary,
                )
            elif export_blocked:
                hard_hard_metrics = summary.metadata.get("hard_hard_invalid_metrics", {})
                validation_metadata = self._hard_hard_validation_metadata(
                    reasons=tuple(summary.metadata.get("hard_hard_invalid_reasons", ())),
                    metrics=hard_hard_metrics if isinstance(hard_hard_metrics, Mapping) else {},
                )
            if cif_path is not None:
                summary = replace(
                    summary,
                    cif_path=cif_path,
                    metadata={
                        **dict(summary.metadata),
                        **({"validation": validation_metadata} if validation_metadata is not None else {}),
                        "cif_export_index": cif_export_index if cif_export_index > 0 else 1,
                    },
                )
            elif validation_metadata is not None:
                summary = replace(
                    summary,
                    metadata={
                        **dict(summary.metadata),
                        "validation": validation_metadata,
                    },
                )
            return summary, best_candidate

        return summaries[0], None

    def generate_pair_candidates(
        self,
        first_record: BatchMonomerRecord,
        second_record: BatchMonomerRecord,
        *,
        out_dir: str | Path | None = None,
        write_cif: bool | None = None,
        cif_export_start_index: int = 0,
    ) -> tuple[tuple[BatchPairSummary, ...], tuple[Candidate, ...], int]:
        first = self.build_monomer(first_record)
        second = self.build_monomer(second_record)

        if not first.ok or not second.ok:
            errors = {}
            if first.error is not None:
                errors["first_error"] = first.error
            if second.error is not None:
                errors["second_error"] = second.error
            return (
                (
                    BatchPairSummary(
                        structure_id=f"{first_record.id}__{second_record.id}",
                        pair_id=f"{first_record.id}__{second_record.id}",
                        pair_mode="unresolved",
                        status="monomer-build-failed",
                        reactant_a_record_id=first_record.id,
                        reactant_b_record_id=second_record.id,
                        metadata=errors,
                    ),
                ),
                (),
                0,
            )

        assert first.monomer is not None
        assert second.monomer is not None
        pair = self._resolve_binary_bridge_pair_context(
            first.monomer,
            second.monomer,
            first_id=first_record.id,
            second_id=second_record.id,
        )
        return self._generate_pair_candidates_from_monomers(
            pair=pair,
            out_dir=out_dir,
            write_cif=write_cif,
            cif_export_start_index=cif_export_start_index,
        )

    def generate_monomer_pair_candidate(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
        *,
        pair_id: str | None = None,
        out_dir: str | Path | None = None,
        write_cif: bool | None = None,
        cif_export_index: int = 0,
    ) -> tuple[BatchPairSummary, Candidate | None]:
        summaries, candidates, _attempted_structures = self.generate_monomer_pair_candidates(
            first,
            second,
            pair_id=pair_id,
            out_dir=None,
            write_cif=False,
            cif_export_start_index=0,
        )
        if not summaries:
            raise ValueError("pair evaluation produced no summaries")

        if candidates:
            best_candidate = max(candidates, key=lambda candidate: candidate.score)
            best_topology = best_candidate.metadata["net_plan"]["topology"]
            best_index = next(
                index
                for index, candidate in enumerate(candidates)
                if candidate.metadata["net_plan"]["topology"] == best_topology
            )
            summary = summaries[best_index]
            effective_write_cif = self.config.write_cif if write_cif is None else write_cif
            cif_path = None
            validation_metadata = None
            export_blocked = bool(summary.metadata.get("hard_hard_invalid_reasons"))
            if effective_write_cif and out_dir is not None and not export_blocked:
                cif_path, validation_metadata = self._write_classified_candidate_cif(
                    out_dir=out_dir,
                    structure_id=summary.structure_id,
                    candidate=best_candidate,
                    monomer_specs=(first, second),
                    provisional_summary=summary,
                )
            elif export_blocked:
                hard_hard_metrics = summary.metadata.get("hard_hard_invalid_metrics", {})
                validation_metadata = self._hard_hard_validation_metadata(
                    reasons=tuple(summary.metadata.get("hard_hard_invalid_reasons", ())),
                    metrics=hard_hard_metrics if isinstance(hard_hard_metrics, Mapping) else {},
                )
            if cif_path is not None:
                summary = replace(
                    summary,
                    cif_path=cif_path,
                    metadata={
                        **dict(summary.metadata),
                        **({"validation": validation_metadata} if validation_metadata is not None else {}),
                        "cif_export_index": cif_export_index if cif_export_index > 0 else 1,
                    },
                )
            elif validation_metadata is not None:
                summary = replace(
                    summary,
                    metadata={
                        **dict(summary.metadata),
                        "validation": validation_metadata,
                    },
                )
            return summary, best_candidate

        return summaries[0], None

    def generate_monomer_pair_candidates(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
        *,
        pair_id: str | None = None,
        out_dir: str | Path | None = None,
        write_cif: bool | None = None,
        cif_export_start_index: int = 0,
    ) -> tuple[tuple[BatchPairSummary, ...], tuple[Candidate, ...], int]:
        pair = self._resolve_binary_bridge_pair_context(
            first,
            second,
            first_id=first.id,
            second_id=second.id,
            pair_id=pair_id,
        )
        return self._generate_pair_candidates_from_monomers(
            pair=pair,
            out_dir=out_dir,
            write_cif=write_cif,
            cif_export_start_index=cif_export_start_index,
        )

    def run_binary_bridge_batch(
        self,
        input_dir: str | Path,
        output_dir: str | Path,
        *,
        template_id: str | None = None,
        max_pairs: int | None = None,
        write_cif: bool | None = None,
        auto_detect_libraries: bool = False,
    ) -> BatchRunSummary:
        template = self._selected_binary_bridge_template(template_id=template_id)
        libraries = self.load_binary_bridge_test_set(
            input_dir,
            template_id=template.id,
            auto_detect=auto_detect_libraries,
        )
        output_root = Path(output_dir)
        output_root.mkdir(parents=True, exist_ok=True)
        manifest_path = output_root / "manifest.jsonl"

        all_records = tuple(record for records in libraries.values() for record in records)
        if self.config.max_workers > 1 and len(all_records) > 1:
            with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                built_results = tuple(executor.map(self.build_monomer, all_records))
            built = {result.record.id: result for result in built_results}
        else:
            built = {record.id: self.build_monomer(record) for record in all_records}
        build_failures = {
            record_id: result.error
            for record_id, result in built.items()
            if result.error is not None
        }
        valid = {
            key: tuple(result.record for result in (built[record.id] for record in records) if result.ok)
            for key, records in libraries.items()
        }

        pair_iterables = self._batch_pair_iterables(valid, template=template)
        pair_tasks = tuple(
            _BatchPairTask(
                first_record=first_record,
                second_record=second_record,
                out_dir=output_root / "cifs",
                write_cif=write_cif,
            )
            for first_record, second_record in self._iter_batch_pairs(pair_iterables, max_pairs=max_pairs)
        )

        attempted_pairs = len(pair_tasks)
        successful_pairs = 0
        attempted_structures = 0
        successful_structures = 0
        cifs_written = 0
        mode_counts: dict[str, int] = {}
        topology_counts: dict[str, int] = {}
        top_results: list[BatchPairSummary] = []
        effective_write_cif = self.config.write_cif if write_cif is None else write_cif
        parallelize_pairs = (
            self.config.max_workers > 1
            and len(pair_tasks) > 1
            and (not effective_write_cif or self.config.max_cif_exports is None)
        )

        with manifest_path.open("w", encoding="utf-8") as manifest:
            if parallelize_pairs:
                if self._supports_process_pair_pool():
                    try:
                        with ProcessPoolExecutor(
                            max_workers=self.config.max_workers,
                            initializer=_init_batch_process_worker,
                            initargs=(self.config, self.reaction_library),
                        ) as executor:
                            pair_results = executor.map(_run_batch_pair_task_in_process, pair_tasks)
                            for summaries, pair_attempted_structures in pair_results:
                                successful_pairs += self._record_batch_pair_results(
                                    manifest=manifest,
                                    summaries=summaries,
                                    top_results=top_results,
                                    mode_counts=mode_counts,
                                    topology_counts=topology_counts,
                                )
                                attempted_structures += pair_attempted_structures
                                successful_structures += sum(1 for summary in summaries if summary.status == "ok")
                                cifs_written += sum(1 for summary in summaries if summary.cif_path is not None)
                    except (NotImplementedError, PermissionError, OSError):
                        with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                            pair_results = executor.map(self._run_batch_pair_task, pair_tasks)
                            for summaries, pair_attempted_structures in pair_results:
                                successful_pairs += self._record_batch_pair_results(
                                    manifest=manifest,
                                    summaries=summaries,
                                    top_results=top_results,
                                    mode_counts=mode_counts,
                                    topology_counts=topology_counts,
                                )
                                attempted_structures += pair_attempted_structures
                                successful_structures += sum(1 for summary in summaries if summary.status == "ok")
                                cifs_written += sum(1 for summary in summaries if summary.cif_path is not None)
                else:
                    with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
                        pair_results = executor.map(self._run_batch_pair_task, pair_tasks)
                        for summaries, pair_attempted_structures in pair_results:
                            successful_pairs += self._record_batch_pair_results(
                                manifest=manifest,
                                summaries=summaries,
                                top_results=top_results,
                                mode_counts=mode_counts,
                                topology_counts=topology_counts,
                            )
                            attempted_structures += pair_attempted_structures
                            successful_structures += sum(1 for summary in summaries if summary.status == "ok")
                            cifs_written += sum(1 for summary in summaries if summary.cif_path is not None)
            else:
                for pair_task in pair_tasks:
                    summaries, pair_attempted_structures = self._run_batch_pair_task(
                        pair_task,
                        cif_export_start_index=cifs_written,
                    )
                    successful_pairs += self._record_batch_pair_results(
                        manifest=manifest,
                        summaries=summaries,
                        top_results=top_results,
                        mode_counts=mode_counts,
                        topology_counts=topology_counts,
                    )
                    attempted_structures += pair_attempted_structures
                    successful_structures += sum(1 for summary in summaries if summary.status == "ok")
                    cifs_written += sum(1 for summary in summaries if summary.cif_path is not None)

        summary = BatchRunSummary(
            input_dir=str(Path(input_dir)),
            output_dir=str(output_root),
            attempted_pairs=attempted_pairs,
            successful_pairs=successful_pairs,
            attempted_structures=attempted_structures,
            successful_structures=successful_structures,
            cifs_written=cifs_written,
            built_monomers=sum(1 for result in built.values() if result.ok),
            failed_monomers=len(build_failures),
            build_failures={key: value for key, value in build_failures.items() if value is not None},
            mode_counts=mode_counts,
            topology_counts=topology_counts,
            manifest_path=str(manifest_path),
            top_results=tuple(top_results),
        )
        staging_dir = output_root / "cifs" / ".staging"
        if staging_dir.is_dir():
            try:
                staging_dir.rmdir()
            except OSError:
                pass
        self._write_summary_report(output_root / "summary.md", summary)
        return summary

    def _supports_process_pair_pool(self) -> bool:
        return self.smiles_monomer_builder is build_rdkit_monomer

    def _iter_batch_pairs(
        self,
        pair_iterables: Iterable[Iterable[tuple[BatchMonomerRecord, BatchMonomerRecord]]],
        *,
        max_pairs: int | None,
    ) -> Iterable[tuple[BatchMonomerRecord, BatchMonomerRecord]]:
        emitted = 0
        for pairs in pair_iterables:
            for first_record, second_record in pairs:
                if max_pairs is not None and emitted >= max_pairs:
                    return
                yield first_record, second_record
                emitted += 1

    def _run_batch_pair_task(
        self,
        task: _BatchPairTask,
        *,
        cif_export_start_index: int = 0,
    ) -> tuple[tuple[BatchPairSummary, ...], int]:
        if self.config.enumerate_all_topologies:
            summaries, _candidates, pair_attempted_structures = self.generate_pair_candidates(
                task.first_record,
                task.second_record,
                out_dir=task.out_dir,
                write_cif=task.write_cif,
                cif_export_start_index=cif_export_start_index,
            )
            return summaries, pair_attempted_structures

        summary, _candidate = self.generate_pair_candidate(
            task.first_record,
            task.second_record,
            out_dir=task.out_dir,
            write_cif=task.write_cif,
            cif_export_index=cif_export_start_index + 1,
        )
        pair_attempted_structures = 1 if summary.status == "ok" else 0
        return (summary,), pair_attempted_structures

    def _record_batch_pair_results(
        self,
        *,
        manifest,
        summaries: tuple[BatchPairSummary, ...],
        top_results: list[BatchPairSummary],
        mode_counts: dict[str, int],
        topology_counts: dict[str, int],
    ) -> int:
        pair_had_success = False
        for summary in summaries:
            if summary.status == "ok":
                pair_had_success = True
                mode_counts[summary.pair_mode] = mode_counts.get(summary.pair_mode, 0) + 1
                if summary.topology_id is not None:
                    topology_counts[summary.topology_id] = topology_counts.get(summary.topology_id, 0) + 1
                top_results.append(summary)
                top_results.sort(
                    key=lambda item: item.score if item.score is not None else float("-inf"),
                    reverse=True,
                )
                del top_results[self.config.retain_top_results :]
            manifest.write(json.dumps(self._json_safe(self._summary_to_dict(summary)), sort_keys=True) + "\n")
        return 1 if pair_had_success else 0

    def run_imine_batch(
        self,
        input_dir: str | Path,
        output_dir: str | Path,
        *,
        max_pairs: int | None = None,
        write_cif: bool | None = None,
        auto_detect_libraries: bool = False,
    ) -> BatchRunSummary:
        return self.run_binary_bridge_batch(
            input_dir,
            output_dir,
            template_id="imine_bridge",
            max_pairs=max_pairs,
            write_cif=write_cif,
            auto_detect_libraries=auto_detect_libraries,
        )

    def _evaluate_pair_topologies(
        self,
        pair: _BinaryBridgePairContext,
    ) -> _TopologyEvaluation:
        first = pair.first
        second = pair.second
        first_connectivity = len(first.motifs)
        second_connectivity = len(second.motifs)
        if first_connectivity >= 3 and second_connectivity >= 3:
            return self._evaluate_node_node_topologies(
                first,
                second,
                template=pair.template,
                first_connectivity=first_connectivity,
                second_connectivity=second_connectivity,
            )
        if min(first_connectivity, second_connectivity) == 2 and max(first_connectivity, second_connectivity) >= 3:
            return self._evaluate_node_linker_topologies(
                first,
                second,
                template=pair.template,
                node_connectivity=max(first_connectivity, second_connectivity),
            )
        raise ValueError(
            f"unsupported monomer connectivity combination {(first_connectivity, second_connectivity)!r}"
        )

    def _run_three_plus_three_pair(self, first: MonomerSpec, second: MonomerSpec) -> Candidate:
        pair = self._resolve_binary_bridge_pair_context(
            first,
            second,
            first_id=first.id,
            second_id=second.id,
        )
        evaluation = self._evaluate_same_connectivity_topologies(
            pair.first,
            pair.second,
            template=pair.template,
            connectivity=3,
        )
        if not evaluation.candidates:
            raise ValueError(f"single-node node-node generation failed: {evaluation.failed_topologies}")
        return max(evaluation.candidates, key=lambda candidate: candidate.score)

    def _run_three_plus_two_pair(self, first: MonomerSpec, second: MonomerSpec) -> Candidate:
        pair = self._resolve_binary_bridge_pair_context(
            first,
            second,
            first_id=first.id,
            second_id=second.id,
        )
        evaluation = self._evaluate_node_linker_topologies(
            pair.first,
            pair.second,
            template=pair.template,
            node_connectivity=3,
        )
        if not evaluation.candidates:
            raise ValueError(f"single-node node-linker generation failed: {evaluation.failed_topologies}")
        return max(evaluation.candidates, key=lambda candidate: candidate.score)

    def _evaluate_same_connectivity_topologies(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
        *,
        template: ReactionTemplate,
        connectivity: int,
    ) -> _TopologyEvaluation:
        return self._evaluate_node_node_topologies(
            first,
            second,
            template=template,
            first_connectivity=connectivity,
            second_connectivity=connectivity,
        )

    def _evaluate_node_node_topologies(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
        *,
        template: ReactionTemplate,
        first_connectivity: int,
        second_connectivity: int,
    ) -> _TopologyEvaluation:
        pair_mode = self._pair_mode_from_connectivities(first_connectivity, second_connectivity)
        topology_ids = self._topology_ids_for_pair(
            connectivities=(first_connectivity, second_connectivity),
            pair_mode="node-node",
        )
        if not topology_ids:
            return _TopologyEvaluation(
                pair_mode=pair_mode,
                candidates=(),
                available_topologies=(),
                failed_topologies=self._topology_unavailable_errors(
                    connectivities=(first_connectivity, second_connectivity),
                    pair_mode="node-node",
                    generic_message=(
                        "no topology-guided node-node builders are available for "
                        f"{min(first_connectivity, second_connectivity)}+{max(first_connectivity, second_connectivity)} "
                        "batch generation"
                    ),
                ),
            )

        raw_candidates: list[tuple[str, Candidate]] = []
        topology_errors: dict[str, str] = {}
        for topology_id in topology_ids:
            topology = get_topology_hint(topology_id)
            if not template.allows_dimensionality(topology.dimensionality):
                topology_errors[topology_id] = (
                    f"template {template.id!r} does not allow {topology.dimensionality!r} generation"
                )
                continue
            try:
                request = self._topology_build_request(
                    first=first,
                    second=second,
                    template=template,
                    topology_id=topology_id,
                    pair_mode="node-node",
                    connectivities=(first_connectivity, second_connectivity),
                )
                raw_candidates.append((topology_id, self._pair_topology_builders.build(request)))
            except Exception as exc:
                topology_errors[topology_id] = f"{type(exc).__name__}: {exc}"

        candidates = tuple(
            self._annotate_topology_selection(
                candidate,
                pair_mode=pair_mode,
                available_topologies=topology_ids,
                failed_topologies=topology_errors,
                selected_topology=topology_id,
                evaluation_mode="all-topologies",
            )
            for topology_id, candidate in raw_candidates
        )
        return _TopologyEvaluation(
            pair_mode=pair_mode,
            candidates=candidates,
            available_topologies=topology_ids,
            failed_topologies=topology_errors,
        )

    def _evaluate_three_plus_three_topologies(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
    ) -> _TopologyEvaluation:
        pair = self._resolve_binary_bridge_pair_context(
            first,
            second,
            first_id=first.id,
            second_id=second.id,
        )
        return self._evaluate_same_connectivity_topologies(
            pair.first,
            pair.second,
            template=pair.template,
            connectivity=3,
        )

    def _evaluate_node_linker_topologies(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
        *,
        template: ReactionTemplate,
        node_connectivity: int,
    ) -> _TopologyEvaluation:
        pair_mode = self._pair_mode_from_connectivities(len(first.motifs), len(second.motifs))
        topology_ids = self._topology_ids_for_pair(
            connectivities=(len(first.motifs), len(second.motifs)),
            pair_mode="node-linker",
        )
        if not topology_ids:
            return _TopologyEvaluation(
                pair_mode=pair_mode,
                candidates=(),
                available_topologies=(),
                failed_topologies=self._topology_unavailable_errors(
                    connectivities=(len(first.motifs), len(second.motifs)),
                    pair_mode="node-linker",
                    generic_message=(
                        "no topology-guided node-linker builders are available for "
                        f"{node_connectivity}+2 batch generation"
                    ),
                ),
            )

        raw_candidates: list[tuple[str, Candidate]] = []
        topology_errors: dict[str, str] = {}
        for topology_id in topology_ids:
            topology = get_topology_hint(topology_id)
            if not template.allows_dimensionality(topology.dimensionality):
                topology_errors[topology_id] = (
                    f"template {template.id!r} does not allow {topology.dimensionality!r} generation"
                )
                continue
            try:
                request = self._topology_build_request(
                    first=first,
                    second=second,
                    template=template,
                    topology_id=topology_id,
                    pair_mode="node-linker",
                    connectivities=(len(first.motifs), len(second.motifs)),
                )
                raw_candidates.append((topology_id, self._pair_topology_builders.build(request)))
            except Exception as exc:
                topology_errors[topology_id] = f"{type(exc).__name__}: {exc}"

        candidates = tuple(
            self._annotate_topology_selection(
                candidate,
                pair_mode=pair_mode,
                available_topologies=topology_ids,
                failed_topologies=topology_errors,
                selected_topology=topology_id,
                evaluation_mode="all-topologies",
            )
            for topology_id, candidate in raw_candidates
        )
        return _TopologyEvaluation(
            pair_mode=pair_mode,
            candidates=candidates,
            available_topologies=topology_ids,
            failed_topologies=topology_errors,
        )

    def _evaluate_three_plus_two_topologies(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
    ) -> _TopologyEvaluation:
        pair = self._resolve_binary_bridge_pair_context(
            first,
            second,
            first_id=first.id,
            second_id=second.id,
        )
        return self._evaluate_node_linker_topologies(
            pair.first,
            pair.second,
            template=pair.template,
            node_connectivity=3,
        )

    def _topology_build_request(
        self,
        *,
        first: MonomerSpec,
        second: MonomerSpec,
        template: ReactionTemplate,
        topology_id: str,
        pair_mode: str,
        connectivities: tuple[int, int],
    ) -> PairTopologyBuildRequest:
        topology = get_topology_hint(topology_id)
        try:
            expanded = self._expand_supported_single_node_topology(topology_id)
            builder_family = "single-node"
        except (KeyError, ValueError):
            expanded = expand_indexed_topology(topology_id)
            builder_family = "indexed"
        return PairTopologyBuildRequest(
            project=COFProject(
                monomers=(first, second),
                allowed_reactions=(template.id,),
                target_dimensionality=topology.dimensionality,
                target_topologies=(topology_id,),
            ),
            monomers=(first, second),
            template=template,
            topology=topology,
            expanded=expanded,
            builder_family=builder_family,
            pair_mode=pair_mode,
            connectivities=connectivities,
        )

    def _supports_two_node_node_node_builder(self, request: PairTopologyBuildRequest) -> bool:
        return (
            request.builder_family == "single-node"
            and request.pair_mode == "node-node"
            and request.topology.dimensionality == "2D"
            and len(request.expanded.node_sites) == 2
        )

    def _supports_two_node_node_linker_builder(self, request: PairTopologyBuildRequest) -> bool:
        return (
            request.builder_family == "single-node"
            and request.pair_mode == "node-linker"
            and request.topology.dimensionality == "2D"
            and len(request.expanded.node_sites) == 2
        )

    def _supports_expanded_node_node_builder(self, request: PairTopologyBuildRequest) -> bool:
        return request.builder_family == "single-node" and request.pair_mode == "node-node"

    def _supports_expanded_node_linker_builder(self, request: PairTopologyBuildRequest) -> bool:
        return request.builder_family == "single-node" and request.pair_mode == "node-linker"

    def _supports_indexed_node_node_builder(self, request: PairTopologyBuildRequest) -> bool:
        return request.builder_family == "indexed" and request.pair_mode == "node-node"

    def _supports_indexed_node_linker_builder(self, request: PairTopologyBuildRequest) -> bool:
        return request.builder_family == "indexed" and request.pair_mode == "node-linker"

    def _build_two_node_node_node_request(self, request: PairTopologyBuildRequest) -> Candidate:
        first, second = request.monomers
        return self._build_two_node_single_node_node_node_candidate(
            first,
            second,
            request.template,
            request.project,
        )

    def _build_two_node_node_linker_request(self, request: PairTopologyBuildRequest) -> Candidate:
        first, second = request.monomers
        node_spec, linker_spec = self._resolve_node_linker_specs(first, second, max(request.connectivities))
        outcome = self._build_single_node_node_linker_outcome(node_spec, linker_spec, request.template, request.topology)
        embedding = self._embed_single_node_node_linker(
            outcome,
            node_spec,
            linker_spec,
            request.template,
            request.topology,
        )
        return self._assemble_single_node_candidate(
            candidate_id=f"{first.id}__{second.id}__{request.topology.id}__single_node_node_linker",
            project=request.project,
            outcome=outcome,
            embedding=embedding,
            monomer_specs={first.id: first, second.id: second},
            template=request.template,
            extra_flags=("batch_single_node_node_linker",),
        )

    def _build_expanded_node_node_request(self, request: PairTopologyBuildRequest) -> Candidate:
        first, second = request.monomers
        expanded = (
            self._node_node_supercell_topology(request.expanded)
            if request.topology.dimensionality == "2D" and len(request.expanded.node_sites) != 2
            else request.expanded
        )
        outcome = self._build_expanded_single_node_node_node_outcome(
            first,
            second,
            request.template,
            request.topology,
            expanded,
        )
        embedding = self._embed_expanded_single_node_node_node(
            outcome,
            first,
            second,
            request.template,
            request.topology,
            expanded,
        )
        return self._assemble_single_node_candidate(
            candidate_id=f"{first.id}__{second.id}__{request.topology.id}__single_node_node_node",
            project=request.project,
            outcome=outcome,
            embedding=embedding,
            monomer_specs={first.id: first, second.id: second},
            template=request.template,
            extra_flags=("batch_single_node_node_node",),
        )

    def _build_expanded_node_linker_request(self, request: PairTopologyBuildRequest) -> Candidate:
        first, second = request.monomers
        node_spec, linker_spec = self._resolve_node_linker_specs(first, second, max(request.connectivities))
        outcome = self._build_expanded_single_node_node_linker_outcome(
            node_spec,
            linker_spec,
            request.template,
            request.topology,
            request.expanded,
        )
        embedding = self._embed_expanded_single_node_node_linker(
            outcome,
            node_spec,
            linker_spec,
            request.template,
            request.topology,
            request.expanded,
        )
        return self._assemble_single_node_candidate(
            candidate_id=f"{first.id}__{second.id}__{request.topology.id}__single_node_node_linker",
            project=request.project,
            outcome=outcome,
            embedding=embedding,
            monomer_specs={first.id: first, second.id: second},
            template=request.template,
            extra_flags=("batch_single_node_node_linker",),
        )

    def _build_indexed_node_linker_request(self, request: PairTopologyBuildRequest) -> Candidate:
        first, second = request.monomers
        expanded = request.expanded
        assert isinstance(expanded, ExpandedIndexedTopology)
        node_spec, linker_spec = self._resolve_node_linker_specs(first, second, max(request.connectivities))
        outcome = self._build_indexed_node_linker_outcome(
            node_spec,
            linker_spec,
            request.template,
            request.topology,
            expanded,
        )
        embedding = self._embed_indexed_node_linker(
            outcome,
            node_spec,
            linker_spec,
            request.template,
            request.topology,
            expanded,
        )
        return self._assemble_single_node_candidate(
            candidate_id=f"{first.id}__{second.id}__{request.topology.id}__indexed_node_linker",
            project=request.project,
            outcome=outcome,
            embedding=embedding,
            monomer_specs={first.id: first, second.id: second},
            template=request.template,
            extra_flags=("batch_indexed_node_linker",),
        )

    def _build_indexed_node_node_request(self, request: PairTopologyBuildRequest) -> Candidate:
        first, second = request.monomers
        expanded = request.expanded
        assert isinstance(expanded, ExpandedIndexedTopology)
        candidates: list[Candidate] = []
        for variant_index, spec_by_node_id in enumerate(
            self._indexed_node_node_spec_assignments(expanded, first, second),
            start=1,
        ):
            outcome = self._build_indexed_node_node_outcome(
                first,
                second,
                request.template,
                request.topology,
                expanded,
                spec_by_node_id,
            )
            embedding = self._embed_indexed_node_node(
                outcome,
                first,
                second,
                request.template,
                request.topology,
                expanded,
                spec_by_node_id,
            )
            candidate = self._assemble_single_node_candidate(
                candidate_id=f"{first.id}__{second.id}__{request.topology.id}__indexed_node_node_{variant_index}",
                project=request.project,
                outcome=outcome,
                embedding=embedding,
                monomer_specs={first.id: first, second.id: second},
                template=request.template,
                extra_flags=("batch_indexed_node_node",),
            )
            candidate_metadata = dict(candidate.metadata)
            candidate_metadata["indexed_node_assignment_variant"] = variant_index
            candidates.append(replace(candidate, metadata=candidate_metadata))
        if not candidates:
            raise ValueError(f"topology {request.topology.id!r} did not produce an indexed node-node assignment")
        return max(candidates, key=lambda candidate: candidate.score)

    def _resolve_node_linker_specs(
        self,
        first: MonomerSpec,
        second: MonomerSpec,
        node_connectivity: int,
    ) -> tuple[MonomerSpec, MonomerSpec]:
        if len(first.motifs) == node_connectivity and len(second.motifs) == 2:
            return first, second
        if len(second.motifs) == node_connectivity and len(first.motifs) == 2:
            return second, first
        raise ValueError(
            f"expected a {node_connectivity}+2 node-linker pair, got {(len(first.motifs), len(second.motifs))!r}"
        )

    def _assemble_single_node_candidate(
        self,
        *,
        candidate_id: str,
        project: COFProject,
        outcome: AssignmentOutcome,
        embedding: EmbeddingResult,
        monomer_specs: Mapping[str, MonomerSpec],
        template: ReactionTemplate,
        extra_flags: tuple[str, ...],
    ) -> Candidate:
        graph = PeriodicProductGraph()
        for instance in outcome.monomer_instances:
            graph.add_monomer(instance)
        template_map = {template.id: template}
        for event in outcome.events:
            graph.add_reaction_event(event, monomer_specs, template_map)
        return self._assemble_candidate(
            candidate_id=candidate_id,
            project=project,
            outcome=outcome,
            graph=graph,
            embedding=embedding,
            monomer_specs=monomer_specs,
            templates=template_map,
            extra_flags=extra_flags,
        )

    def _build_single_node_node_linker_outcome(
        self,
        node_spec: MonomerSpec,
        linker_spec: MonomerSpec,
        template,
        topology,
    ) -> AssignmentOutcome:
        layout = resolve_single_node_topology_layout(topology.id)
        effective_template_id = template.id
        node_motif_ids = self._ordered_planar_motif_ids(node_spec)
        edge_directions = layout.directions
        image_sequence = ((0, 0, 0), (-1, 0, 0), (0, -1, 0))
        linker_orders = [
            self._linker_layout_for_direction(linker_spec, direction, template_id=effective_template_id).motif_ids
            for direction in edge_directions
        ]

        net_plan = NetPlan(
            topology=topology,
            monomer_ids=(node_spec.id, linker_spec.id),
            reaction_ids=(template.id,),
            metadata={
                "planning_mode": "batch-single-node-node-linker",
                "node_monomer_id": node_spec.id,
                "linker_monomer_id": linker_spec.id,
                "connectivities": (len(node_spec.motifs), len(linker_spec.motifs)),
                "topology_metadata": dict(topology.metadata),
                "topology_family": "single-node-2d",
                "topology_inference_mode": layout.inference_mode,
            },
        )
        assignment_plan = AssignmentPlan(
            net_plan=net_plan,
            slot_to_monomer={
                "node_a": node_spec.id,
                "node_b": node_spec.id,
                "linker_1": linker_spec.id,
                "linker_2": linker_spec.id,
                "linker_3": linker_spec.id,
            },
            slot_to_conformer={
                "node_a": node_spec.conformer_ids[0],
                "node_b": node_spec.conformer_ids[0],
                "linker_1": linker_spec.conformer_ids[0],
                "linker_2": linker_spec.conformer_ids[0],
                "linker_3": linker_spec.conformer_ids[0],
            }
            if node_spec.conformer_ids and linker_spec.conformer_ids
            else {},
            metadata={"assignment_mode": "batch-single-node-node-linker"},
        )
        instances = (
            MonomerInstance(id="m1", monomer_id=node_spec.id, conformer_id=node_spec.conformer_ids[0] if node_spec.conformer_ids else None, metadata={"slot_id": "node_a", "role": "node", "sublattice": "A"}),
            MonomerInstance(id="m2", monomer_id=node_spec.id, conformer_id=node_spec.conformer_ids[0] if node_spec.conformer_ids else None, metadata={"slot_id": "node_b", "role": "node", "sublattice": "B"}),
            MonomerInstance(id="m3", monomer_id=linker_spec.id, conformer_id=linker_spec.conformer_ids[0] if linker_spec.conformer_ids else None, metadata={"slot_id": "linker_1", "role": "linker", "edge_index": 1}),
            MonomerInstance(id="m4", monomer_id=linker_spec.id, conformer_id=linker_spec.conformer_ids[0] if linker_spec.conformer_ids else None, metadata={"slot_id": "linker_2", "role": "linker", "edge_index": 2}),
            MonomerInstance(id="m5", monomer_id=linker_spec.id, conformer_id=linker_spec.conformer_ids[0] if linker_spec.conformer_ids else None, metadata={"slot_id": "linker_3", "role": "linker", "edge_index": 3}),
        )

        events: list[ReactionEvent] = []
        for edge_index, (image, linker_motif_ids) in enumerate(zip(image_sequence, linker_orders), start=1):
            node_motif_id = node_motif_ids[edge_index - 1]
            linker_instance_id = f"m{edge_index + 2}"
            linker_near_a, linker_near_b = linker_motif_ids
            events.append(
                ReactionEvent(
                    id=f"rxn{len(events) + 1}",
                    template_id=template.id,
                    participants=(
                        MotifRef(monomer_instance_id="m1", monomer_id=node_spec.id, motif_id=node_motif_id),
                        MotifRef(monomer_instance_id=linker_instance_id, monomer_id=linker_spec.id, motif_id=linker_near_a),
                    ),
                    metadata={"edge_index": edge_index, "endpoint": "node_a"},
                )
            )
            events.append(
                ReactionEvent(
                    id=f"rxn{len(events) + 1}",
                    template_id=template.id,
                    participants=(
                        MotifRef(monomer_instance_id="m2", monomer_id=node_spec.id, motif_id=node_motif_id, periodic_image=image),
                        MotifRef(monomer_instance_id=linker_instance_id, monomer_id=linker_spec.id, motif_id=linker_near_b),
                    ),
                    metadata={"edge_index": edge_index, "endpoint": "node_b", "node_image": image},
                )
            )

        return AssignmentOutcome(
            assignment_plan=assignment_plan,
            monomer_instances=instances,
            events=tuple(events),
            unreacted_motifs=(),
            consumed_count=sum(len(event.participants) for event in events),
        )

    def _embed_single_node_node_linker(
        self,
        outcome: AssignmentOutcome,
        node_spec: MonomerSpec,
        linker_spec: MonomerSpec,
        template,
        topology,
    ) -> EmbeddingResult:
        layout = resolve_single_node_topology_layout(topology.id)
        effective_template_id = template.id
        bridge_target = self._template_target_distance(template)
        canonical_a = layout.directions
        canonical_b = tuple(scale(direction, -1.0) for direction in canonical_a)
        node_rotation_a, node_motif_ids, node_offsets_a = self._rotation_for_planar_motifs(
            node_spec,
            canonical_a,
            template_id=effective_template_id,
        )
        node_rotation_b, _node_motif_ids_b, node_offsets_b = self._rotation_for_planar_motifs(
            node_spec,
            canonical_b,
            template_id=effective_template_id,
        )
        node_motifs = {motif.id: motif for motif in node_spec.motifs}
        edge_directions: list[Vec3] = []
        linker_layouts: list[_LinkerLayout] = []
        edge_vectors: list[Vec3] = []
        reactive_site_distances: list[float] = []

        for node_motif_id in node_motif_ids:
            node_motif = node_motifs[node_motif_id]
            direction = self._node_linker_edge_direction(
                motif=node_motif,
                monomer=node_spec,
                node_rotation_a=node_rotation_a,
                node_rotation_b=node_rotation_b,
                template_id=effective_template_id,
            )
            linker_layout = self._linker_layout_for_direction(linker_spec, direction, template_id=effective_template_id)
            reactive_site_span = bridge_target + linker_layout.span + bridge_target
            edge_vector = self._generalized_single_node_edge_vector(
                node_motif,
                monomer=node_spec,
                node_rotation_a=node_rotation_a,
                node_rotation_b=node_rotation_b,
                direction=direction,
                reactive_site_span=reactive_site_span,
                template_id=effective_template_id,
            )
            edge_directions.append(direction)
            linker_layouts.append(linker_layout)
            edge_vectors.append(edge_vector)
            reactive_site_distances.append(norm(edge_vector))

        cell = self._generalized_single_node_cell(tuple(edge_vectors))
        center_a = scale(add(cell[0], cell[1]), 1.0 / 3.0)
        center_b = add(center_a, edge_vectors[0])
        monomer_poses: dict[str, Pose] = {
            "m1": Pose(translation=center_a, rotation_matrix=node_rotation_a),
            "m2": Pose(translation=center_b, rotation_matrix=node_rotation_b),
        }
        pose_details: dict[str, object] = {
            "m1": {
                "slot_id": "node_a",
                "translation": center_a,
                "rotation_matrix": node_rotation_a,
                "motif_count": len(node_spec.motifs),
                "role": "node",
                "sublattice": "A",
                "radial_offsets": tuple(round(value, 6) for value in node_offsets_a),
            },
            "m2": {
                "slot_id": "node_b",
                "translation": center_b,
                "rotation_matrix": node_rotation_b,
                "motif_count": len(node_spec.motifs),
                "role": "node",
                "sublattice": "B",
                "radial_offsets": tuple(round(value, 6) for value in node_offsets_b),
            },
        }

        node_images = ((0, 0, 0), (-1, 0, 0), (0, -1, 0))
        for edge_index, (direction, image, linker_layout) in enumerate(zip(edge_directions, node_images, linker_layouts), start=1):
            node_motif_id = node_motif_ids[edge_index - 1]
            node_a_origin = self._world_motif_origin(
                cell,
                monomer_poses["m1"],
                node_spec,
                node_motifs[node_motif_id],
                (0, 0, 0),
                template_id=effective_template_id,
            )
            node_b_origin = self._world_motif_origin(
                cell,
                monomer_poses["m2"],
                node_spec,
                node_motifs[node_motif_id],
                image,
                template_id=effective_template_id,
            )
            bridge_direction = self._safe_normalize(sub(node_b_origin, node_a_origin))
            desired_a = add(node_a_origin, scale(bridge_direction, bridge_target))
            desired_b = add(node_b_origin, scale(bridge_direction, -bridge_target))
            linker_local_a, linker_local_b = linker_layout.local_origins
            translation_a = sub(desired_a, matmul_vec(linker_layout.rotation, linker_local_a))
            translation_b = sub(desired_b, matmul_vec(linker_layout.rotation, linker_local_b))
            linker_translation = scale(add(translation_a, translation_b), 0.5)
            linker_instance_id = f"m{edge_index + 2}"
            monomer_poses[linker_instance_id] = Pose(
                translation=linker_translation,
                rotation_matrix=linker_layout.rotation,
            )
            pose_details[linker_instance_id] = {
                "slot_id": f"linker_{edge_index}",
                "translation": linker_translation,
                "rotation_matrix": linker_layout.rotation,
                "motif_count": len(linker_spec.motifs),
                "role": "linker",
                "edge_index": edge_index,
                "node_b_image": image,
                "edge_direction": tuple(round(value, 6) for value in direction),
                "layout_span": linker_layout.span,
                "layout_offsets": (round(linker_layout.offset_a, 6), round(linker_layout.offset_b, 6)),
                "translation_mismatch": round(norm(sub(translation_b, translation_a)), 6),
            }

        return EmbeddingResult(
            state=AssemblyState(cell=cell, monomer_poses=monomer_poses, stacking_state="disabled"),
            metadata={
                "mode": "topology-guided",
                "topology": topology.id,
                "target_distance": bridge_target,
                "reactive_site_distance": sum(reactive_site_distances) / len(reactive_site_distances),
                "edge_reactive_site_distances": tuple(round(value, 6) for value in reactive_site_distances),
                "cell_kind": self._single_node_cell_kind(cell),
                "stacking_enabled": False,
                "placement_mode": "node-linker-single-node",
                "node_instance_count": 2,
                "linker_instance_count": 3,
                "topology_family": "single-node-2d",
                "topology_inference_mode": layout.inference_mode,
                "poses": pose_details,
            },
        )

    def _build_expanded_single_node_node_linker_outcome(
        self,
        node_spec: MonomerSpec,
        linker_spec: MonomerSpec,
        template,
        topology,
        expanded: ExpandedSingleNodeTopology,
    ) -> AssignmentOutcome:
        if len(linker_spec.motifs) != 2:
            raise ValueError("single-node node-linker placement requires a ditopic linker")
        effective_template_id = template.id

        node_placements = self._expanded_node_site_placements(
            expanded.node_sites,
            {node_site.id: node_spec for node_site in expanded.node_sites},
            dimensionality=topology.dimensionality,
            template_id=effective_template_id,
        )
        linker_instance_ids = {
            edge.id: f"l{index + 1}" for index, edge in enumerate(expanded.edge_sites)
        }

        slot_to_monomer: dict[str, str] = {}
        slot_to_conformer: dict[str, str] = {}
        instances: list[MonomerInstance] = []
        for node_site in expanded.node_sites:
            slot_id = f"node_{node_site.id}"
            slot_to_monomer[slot_id] = node_spec.id
            if node_spec.conformer_ids:
                slot_to_conformer[slot_id] = node_spec.conformer_ids[0]
            instances.append(
                MonomerInstance(
                    id=node_site.id,
                    monomer_id=node_spec.id,
                    conformer_id=node_spec.conformer_ids[0] if node_spec.conformer_ids else None,
                    metadata={
                        "slot_id": slot_id,
                        "role": "node",
                        "topology_node_id": node_site.id,
                        "sublattice": node_site.sublattice,
                        "fractional_position": node_site.fractional_position,
                    },
                )
            )
        for edge in expanded.edge_sites:
            linker_instance_id = linker_instance_ids[edge.id]
            slot_id = f"linker_{edge.id}"
            slot_to_monomer[slot_id] = linker_spec.id
            if linker_spec.conformer_ids:
                slot_to_conformer[slot_id] = linker_spec.conformer_ids[0]
            instances.append(
                MonomerInstance(
                    id=linker_instance_id,
                    monomer_id=linker_spec.id,
                    conformer_id=linker_spec.conformer_ids[0] if linker_spec.conformer_ids else None,
                    metadata={
                        "slot_id": slot_id,
                        "role": "linker",
                        "topology_edge_id": edge.id,
                        "end_image": edge.end_image,
                        "center_fractional": edge.center_fractional,
                    },
                )
            )

        net_plan = NetPlan(
            topology=topology,
            monomer_ids=(node_spec.id, linker_spec.id),
            reaction_ids=(template.id,),
            metadata={
                "planning_mode": "batch-single-node-node-linker-expanded",
                "node_monomer_id": node_spec.id,
                "linker_monomer_id": linker_spec.id,
                "connectivities": (len(node_spec.motifs), len(linker_spec.motifs)),
                "topology_metadata": dict(topology.metadata),
                "topology_family": "single-node-2d",
                "expanded_node_count": len(expanded.node_sites),
                "expanded_edge_count": len(expanded.edge_sites),
            },
        )
        assignment_plan = AssignmentPlan(
            net_plan=net_plan,
            slot_to_monomer=slot_to_monomer,
            slot_to_conformer=slot_to_conformer,
            metadata={"assignment_mode": "batch-single-node-node-linker-expanded"},
        )

        events: list[ReactionEvent] = []
        for edge in expanded.edge_sites:
            direction = self._safe_normalize(edge.base_vector)
            linker_layout = self._linker_layout_for_direction(linker_spec, direction, template_id=effective_template_id)
            linker_instance_id = linker_instance_ids[edge.id]
            start_motif_id = str(
                node_placements[edge.start_node_id].motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")]
            )
            end_motif_id = str(
                node_placements[edge.end_node_id].motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")]
            )
            linker_near_start, linker_near_end = linker_layout.motif_ids
            events.append(
                ReactionEvent(
                    id=f"rxn{len(events) + 1}",
                    template_id=template.id,
                    participants=(
                        MotifRef(
                            monomer_instance_id=edge.start_node_id,
                            monomer_id=node_spec.id,
                            motif_id=start_motif_id,
                        ),
                        MotifRef(
                            monomer_instance_id=linker_instance_id,
                            monomer_id=linker_spec.id,
                            motif_id=linker_near_start,
                        ),
                    ),
                    metadata={"topology_edge_id": edge.id, "endpoint": "start"},
                )
            )
            events.append(
                ReactionEvent(
                    id=f"rxn{len(events) + 1}",
                    template_id=template.id,
                    participants=(
                        MotifRef(
                            monomer_instance_id=edge.end_node_id,
                            monomer_id=node_spec.id,
                            motif_id=end_motif_id,
                            periodic_image=edge.end_image,
                        ),
                        MotifRef(
                            monomer_instance_id=linker_instance_id,
                            monomer_id=linker_spec.id,
                            motif_id=linker_near_end,
                        ),
                    ),
                    metadata={"topology_edge_id": edge.id, "endpoint": "end", "node_image": edge.end_image},
                )
            )

        return AssignmentOutcome(
            assignment_plan=assignment_plan,
            monomer_instances=tuple(instances),
            events=tuple(events),
            unreacted_motifs=(),
            consumed_count=sum(len(event.participants) for event in events),
        )

    def _embed_expanded_single_node_node_linker(
        self,
        outcome: AssignmentOutcome,
        node_spec: MonomerSpec,
        linker_spec: MonomerSpec,
        template,
        topology,
        expanded: ExpandedSingleNodeTopology,
    ) -> EmbeddingResult:
        effective_template_id = template.id
        bridge_target = self._template_target_distance(template)
        node_placements = self._expanded_node_site_placements(
            expanded.node_sites,
            {node_site.id: node_spec for node_site in expanded.node_sites},
            dimensionality=topology.dimensionality,
            template_id=effective_template_id,
        )
        node_motifs = {motif.id: motif for motif in node_spec.motifs}
        linker_layouts: dict[str, _LinkerLayout] = {}

        target_edge_vectors: dict[str, Vec3] = {}
        reactive_site_distances: list[float] = []
        for edge in expanded.edge_sites:
            direction = self._safe_normalize(edge.base_vector)
            start_spec = node_spec
            end_spec = node_spec
            start_placement = node_placements[edge.start_node_id]
            end_placement = node_placements[edge.end_node_id]
            start_motif = node_motifs[
                str(start_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")])
            ]
            end_motif = node_motifs[
                str(end_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")])
            ]
            if topology.dimensionality == "3D":
                start_normal = matmul_vec(start_placement.rotation, start_motif.frame.normal)
                end_normal = matmul_vec(end_placement.rotation, end_motif.frame.normal)
                target_normal = self._bridge_plane_normal(direction, start_normal, end_normal)
                linker_layout = self._linker_layout_for_axes(
                    linker_spec,
                    direction,
                    target_normal,
                    template_id=effective_template_id,
                )
            else:
                linker_layout = self._linker_layout_for_direction(linker_spec, direction, template_id=effective_template_id)
            linker_layouts[edge.id] = linker_layout
            reactive_site_span = bridge_target + linker_layout.span + bridge_target
            edge_vector = add(
                scale(direction, reactive_site_span),
                sub(
                    matmul_vec(
                        start_placement.rotation,
                        effective_motif_origin(effective_template_id, node_spec, start_motif),
                    ),
                    matmul_vec(
                        end_placement.rotation,
                        effective_motif_origin(effective_template_id, node_spec, end_motif),
                    ),
                ),
            )
            target_edge_vectors[edge.id] = edge_vector
            reactive_site_distances.append(norm(edge_vector))

        cell = self._fit_expanded_single_node_cell(
            expanded,
            target_edge_vectors,
            dimensionality=topology.dimensionality,
        )
        monomer_poses: dict[str, Pose] = {}
        pose_details: dict[str, object] = {}
        for node_site in expanded.node_sites:
            placement = node_placements[node_site.id]
            translation = self._fractional_position_in_cell(cell, node_site.fractional_position)
            monomer_poses[node_site.id] = Pose(translation=translation, rotation_matrix=placement.rotation)
            pose_details[node_site.id] = {
                "slot_id": f"node_{node_site.id}",
                "translation": translation,
                "rotation_matrix": placement.rotation,
                "motif_count": len(node_spec.motifs),
                "role": "node",
                "topology_node_id": node_site.id,
                "sublattice": node_site.sublattice,
                "fractional_position": node_site.fractional_position,
                "radial_offsets": tuple(round(value, 6) for value in placement.offsets),
            }

        for edge in expanded.edge_sites:
            linker_layout = linker_layouts[edge.id]
            linker_instance_id = next(
                instance.id
                for instance in outcome.monomer_instances
                if instance.metadata.get("topology_edge_id") == edge.id
            )
            start_spec = node_spec
            end_spec = node_spec
            start_placement = node_placements[edge.start_node_id]
            end_placement = node_placements[edge.end_node_id]
            start_motif = node_motifs[
                str(start_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")])
            ]
            end_motif = node_motifs[
                str(end_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")])
            ]
            start_origin = self._world_motif_origin(
                cell,
                monomer_poses[edge.start_node_id],
                node_spec,
                start_motif,
                (0, 0, 0),
                template_id=effective_template_id,
            )
            end_origin = self._world_motif_origin(
                cell,
                monomer_poses[edge.end_node_id],
                node_spec,
                end_motif,
                edge.end_image,
                template_id=effective_template_id,
            )
            bridge_direction = self._safe_normalize(sub(end_origin, start_origin))
            desired_start = add(start_origin, scale(bridge_direction, bridge_target))
            desired_end = add(end_origin, scale(bridge_direction, -bridge_target))
            linker_local_a, linker_local_b = linker_layout.local_origins
            translation_a = sub(desired_start, matmul_vec(linker_layout.rotation, linker_local_a))
            translation_b = sub(desired_end, matmul_vec(linker_layout.rotation, linker_local_b))
            linker_translation = scale(add(translation_a, translation_b), 0.5)
            monomer_poses[linker_instance_id] = Pose(
                translation=linker_translation,
                rotation_matrix=linker_layout.rotation,
            )
            pose_details[linker_instance_id] = {
                "slot_id": f"linker_{edge.id}",
                "translation": linker_translation,
                "rotation_matrix": linker_layout.rotation,
                "motif_count": len(linker_spec.motifs),
                "role": "linker",
                "topology_edge_id": edge.id,
                "node_image": edge.end_image,
                "edge_direction": tuple(round(value, 6) for value in self._safe_normalize(edge.base_vector)),
                "layout_span": linker_layout.span,
                "layout_offsets": (round(linker_layout.offset_a, 6), round(linker_layout.offset_b, 6)),
                "translation_mismatch": round(norm(sub(translation_b, translation_a)), 6),
            }

        return EmbeddingResult(
            state=AssemblyState(cell=cell, monomer_poses=monomer_poses, stacking_state="disabled"),
            metadata={
                "mode": "topology-guided",
                "topology": topology.id,
                "target_distance": bridge_target,
                "reactive_site_distance": sum(reactive_site_distances) / len(reactive_site_distances),
                "edge_reactive_site_distances": tuple(round(value, 6) for value in reactive_site_distances),
                "cell_kind": self._topology_cell_kind(cell, topology.dimensionality),
                "stacking_enabled": False,
                "placement_mode": (
                    "node-linker-single-node-3d" if topology.dimensionality == "3D" else "node-linker-single-node-expanded"
                ),
                "node_instance_count": len(expanded.node_sites),
                "linker_instance_count": len(expanded.edge_sites),
                "topology_family": f"single-node-{topology.dimensionality.lower()}",
                "expanded_node_orbit_size": len(expanded.node_sites),
                "expanded_edge_count": len(expanded.edge_sites),
                "poses": pose_details,
            },
        )

    def _build_expanded_single_node_node_node_outcome(
        self,
        amine_spec: MonomerSpec,
        aldehyde_spec: MonomerSpec,
        template,
        topology,
        expanded: ExpandedSingleNodeTopology,
    ) -> AssignmentOutcome:
        if not expanded.is_bipartite:
            raise ValueError(
                f"topology {topology.id!r} is not bipartite and cannot host a single-node node-node alternation"
            )
        effective_template_id = template.id

        spec_by_node_id = {
            node_site.id: amine_spec if node_site.sublattice == 0 else aldehyde_spec
            for node_site in expanded.node_sites
        }
        placements = self._expanded_node_site_placements(
            expanded.node_sites,
            spec_by_node_id,
            dimensionality=topology.dimensionality,
            template_id=effective_template_id,
        )

        slot_to_monomer: dict[str, str] = {}
        slot_to_conformer: dict[str, str] = {}
        instances: list[MonomerInstance] = []
        for node_site in expanded.node_sites:
            spec = spec_by_node_id[node_site.id]
            slot_id = f"node_{node_site.id}"
            slot_to_monomer[slot_id] = spec.id
            if spec.conformer_ids:
                slot_to_conformer[slot_id] = spec.conformer_ids[0]
            instances.append(
                MonomerInstance(
                    id=node_site.id,
                    monomer_id=spec.id,
                    conformer_id=spec.conformer_ids[0] if spec.conformer_ids else None,
                    metadata={
                        "slot_id": slot_id,
                        "role": "node",
                        "topology_node_id": node_site.id,
                        "sublattice": node_site.sublattice,
                        "fractional_position": node_site.fractional_position,
                    },
                )
            )

        net_plan = NetPlan(
            topology=topology,
            monomer_ids=(amine_spec.id, aldehyde_spec.id),
            reaction_ids=(template.id,),
            metadata={
                "planning_mode": "batch-single-node-node-node-expanded",
                "connectivities": (len(amine_spec.motifs), len(aldehyde_spec.motifs)),
                "topology_metadata": dict(topology.metadata),
                "topology_family": "single-node-2d",
                "expanded_node_count": len(expanded.node_sites),
                "expanded_edge_count": len(expanded.edge_sites),
            },
        )
        assignment_plan = AssignmentPlan(
            net_plan=net_plan,
            slot_to_monomer=slot_to_monomer,
            slot_to_conformer=slot_to_conformer,
            metadata={"assignment_mode": "batch-single-node-node-node-expanded"},
        )

        events: list[ReactionEvent] = []
        for edge in expanded.edge_sites:
            start_spec = spec_by_node_id[edge.start_node_id]
            end_spec = spec_by_node_id[edge.end_node_id]
            start_motif_id = str(
                placements[edge.start_node_id].motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")]
            )
            end_motif_id = str(
                placements[edge.end_node_id].motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")]
            )
            events.append(
                ReactionEvent(
                    id=f"rxn{len(events) + 1}",
                    template_id=template.id,
                    participants=(
                        MotifRef(
                            monomer_instance_id=edge.start_node_id,
                            monomer_id=start_spec.id,
                            motif_id=start_motif_id,
                        ),
                        MotifRef(
                            monomer_instance_id=edge.end_node_id,
                            monomer_id=end_spec.id,
                            motif_id=end_motif_id,
                            periodic_image=edge.end_image,
                        ),
                    ),
                    metadata={"topology_edge_id": edge.id, "node_image": edge.end_image},
                )
            )

        return AssignmentOutcome(
            assignment_plan=assignment_plan,
            monomer_instances=tuple(instances),
            events=tuple(events),
            unreacted_motifs=(),
            consumed_count=sum(len(event.participants) for event in events),
        )

    def _embed_expanded_single_node_node_node(
        self,
        outcome: AssignmentOutcome,
        amine_spec: MonomerSpec,
        aldehyde_spec: MonomerSpec,
        template,
        topology,
        expanded: ExpandedSingleNodeTopology,
    ) -> EmbeddingResult:
        effective_template_id = template.id
        bridge_target = self._template_target_distance(template)
        spec_by_node_id = {
            node_site.id: amine_spec if node_site.sublattice == 0 else aldehyde_spec
            for node_site in expanded.node_sites
        }
        placements = self._expanded_node_site_placements(
            expanded.node_sites,
            spec_by_node_id,
            dimensionality=topology.dimensionality,
            template_id=effective_template_id,
        )
        motifs_by_spec = {
            amine_spec.id: {motif.id: motif for motif in amine_spec.motifs},
            aldehyde_spec.id: {motif.id: motif for motif in aldehyde_spec.motifs},
        }

        target_edge_vectors: dict[str, Vec3] = {}
        reactive_site_distances: list[float] = []
        for edge in expanded.edge_sites:
            direction = self._safe_normalize(edge.base_vector)
            start_spec = spec_by_node_id[edge.start_node_id]
            end_spec = spec_by_node_id[edge.end_node_id]
            start_placement = placements[edge.start_node_id]
            end_placement = placements[edge.end_node_id]
            start_motif = motifs_by_spec[start_spec.id][
                str(start_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")])
            ]
            end_motif = motifs_by_spec[end_spec.id][
                str(end_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")])
            ]
            edge_vector = add(
                scale(direction, bridge_target),
                sub(
                    matmul_vec(
                        start_placement.rotation,
                        effective_motif_origin(effective_template_id, start_spec, start_motif),
                    ),
                    matmul_vec(
                        end_placement.rotation,
                        effective_motif_origin(effective_template_id, end_spec, end_motif),
                    ),
                ),
            )
            target_edge_vectors[edge.id] = edge_vector
            reactive_site_distances.append(norm(edge_vector))

        cell = self._fit_expanded_single_node_cell(
            expanded,
            target_edge_vectors,
            dimensionality=topology.dimensionality,
        )
        monomer_poses: dict[str, Pose] = {}
        pose_details: dict[str, object] = {}
        for node_site in expanded.node_sites:
            spec = spec_by_node_id[node_site.id]
            placement = placements[node_site.id]
            translation = self._fractional_position_in_cell(cell, node_site.fractional_position)
            monomer_poses[node_site.id] = Pose(translation=translation, rotation_matrix=placement.rotation)
            pose_details[node_site.id] = {
                "slot_id": f"node_{node_site.id}",
                "translation": translation,
                "rotation_matrix": placement.rotation,
                "motif_count": len(spec.motifs),
                "role": "node",
                "topology_node_id": node_site.id,
                "sublattice": node_site.sublattice,
                "fractional_position": node_site.fractional_position,
                "radial_offsets": tuple(round(value, 6) for value in placement.offsets),
            }

        return EmbeddingResult(
            state=AssemblyState(cell=cell, monomer_poses=monomer_poses, stacking_state="disabled"),
            metadata={
                "mode": "topology-guided",
                "topology": topology.id,
                "target_distance": bridge_target,
                "reactive_site_distance": sum(reactive_site_distances) / len(reactive_site_distances),
                "edge_reactive_site_distances": tuple(round(value, 6) for value in reactive_site_distances),
                "cell_kind": self._topology_cell_kind(cell, topology.dimensionality),
                "stacking_enabled": False,
                "placement_mode": (
                    "single-node-bipartite"
                    if topology.dimensionality == "2D" and len(expanded.node_sites) == 2
                    else (
                        "single-node-3d-node-node"
                        if topology.dimensionality == "3D"
                        else "single-node-expanded-node-node"
                    )
                ),
                "node_instance_count": len(expanded.node_sites),
                "topology_family": f"single-node-{topology.dimensionality.lower()}",
                "expanded_node_orbit_size": len(expanded.node_sites),
                "expanded_edge_count": len(expanded.edge_sites),
                "poses": pose_details,
            },
        )

    def _indexed_node_node_spec_assignments(
        self,
        expanded: ExpandedIndexedTopology,
        first: MonomerSpec,
        second: MonomerSpec,
    ) -> tuple[Mapping[str, MonomerSpec], ...]:
        first_connectivity = len(first.motifs)
        second_connectivity = len(second.motifs)
        node_by_id = {node_site.id: node_site for node_site in expanded.node_sites}

        if first_connectivity != second_connectivity:
            spec_by_node_id: dict[str, MonomerSpec] = {}
            for node_site in expanded.node_sites:
                if node_site.connectivity == first_connectivity:
                    spec_by_node_id[node_site.id] = first
                elif node_site.connectivity == second_connectivity:
                    spec_by_node_id[node_site.id] = second
                else:
                    raise ValueError(
                        f"topology {expanded.topology_id!r} contains node connectivity {node_site.connectivity}, "
                        f"which does not match pair connectivities {(first_connectivity, second_connectivity)!r}"
                    )
            self._validate_indexed_node_node_assignment(expanded, spec_by_node_id)
            return (spec_by_node_id,)

        source_indices = tuple(sorted({int(node_site.source_node_index) for node_site in expanded.node_sites}))
        if len(source_indices) != 2:
            raise ValueError(
                f"equal-connectivity indexed node-node generation for topology {expanded.topology_id!r} currently "
                "requires exactly two source node orbits"
            )
        for edge in expanded.edge_sites:
            start_source = int(node_by_id[edge.start_node_id].source_node_index)
            end_source = int(node_by_id[edge.end_node_id].source_node_index)
            if start_source == end_source:
                raise ValueError(
                    f"topology {expanded.topology_id!r} requires a more general bipartite role coloring than the "
                    "current indexed node-node builder provides"
                )

        assignments: list[Mapping[str, MonomerSpec]] = []
        for first_source, second_source in ((source_indices[0], source_indices[1]), (source_indices[1], source_indices[0])):
            spec_by_node_id = {
                node_site.id: first if int(node_site.source_node_index) == first_source else second
                for node_site in expanded.node_sites
            }
            self._validate_indexed_node_node_assignment(expanded, spec_by_node_id)
            assignments.append(spec_by_node_id)
        return tuple(assignments)

    def _validate_indexed_node_node_assignment(
        self,
        expanded: ExpandedIndexedTopology,
        spec_by_node_id: Mapping[str, MonomerSpec],
    ) -> None:
        for edge in expanded.edge_sites:
            start_spec = spec_by_node_id[edge.start_node_id]
            end_spec = spec_by_node_id[edge.end_node_id]
            if start_spec.id == end_spec.id:
                raise ValueError(
                    f"topology {expanded.topology_id!r} edge {edge.id!r} would connect identical node roles in "
                    "a node-node pairing"
                )

    def _build_indexed_node_linker_outcome(
        self,
        node_spec: MonomerSpec,
        linker_spec: MonomerSpec,
        template,
        topology,
        expanded: ExpandedIndexedTopology,
    ) -> AssignmentOutcome:
        if len(linker_spec.motifs) != 2:
            raise ValueError("indexed node-linker placement requires a ditopic linker")
        effective_template_id = template.id

        node_placements = self._expanded_node_site_placements(
            expanded.node_sites,
            {node_site.id: node_spec for node_site in expanded.node_sites},
            dimensionality=topology.dimensionality,
            template_id=effective_template_id,
        )
        linker_instance_ids = {
            edge.id: f"l{index + 1}" for index, edge in enumerate(expanded.edge_sites)
        }

        slot_to_monomer: dict[str, str] = {}
        slot_to_conformer: dict[str, str] = {}
        instances: list[MonomerInstance] = []
        for node_site in expanded.node_sites:
            slot_id = f"node_{node_site.id}"
            slot_to_monomer[slot_id] = node_spec.id
            if node_spec.conformer_ids:
                slot_to_conformer[slot_id] = node_spec.conformer_ids[0]
            instances.append(
                MonomerInstance(
                    id=node_site.id,
                    monomer_id=node_spec.id,
                    conformer_id=node_spec.conformer_ids[0] if node_spec.conformer_ids else None,
                    metadata={
                        "slot_id": slot_id,
                        "role": "node",
                        "topology_node_id": node_site.id,
                        "source_node_index": node_site.source_node_index,
                        "source_label": node_site.source_label,
                        "connectivity": node_site.connectivity,
                        "fractional_position": node_site.fractional_position,
                    },
                )
            )
        for edge in expanded.edge_sites:
            linker_instance_id = linker_instance_ids[edge.id]
            slot_id = f"linker_{edge.id}"
            slot_to_monomer[slot_id] = linker_spec.id
            if linker_spec.conformer_ids:
                slot_to_conformer[slot_id] = linker_spec.conformer_ids[0]
            instances.append(
                MonomerInstance(
                    id=linker_instance_id,
                    monomer_id=linker_spec.id,
                    conformer_id=linker_spec.conformer_ids[0] if linker_spec.conformer_ids else None,
                    metadata={
                        "slot_id": slot_id,
                        "role": "linker",
                        "topology_edge_id": edge.id,
                        "end_image": edge.end_image,
                        "center_fractional": edge.center_fractional,
                    },
                )
            )

        net_plan = NetPlan(
            topology=topology,
            monomer_ids=(node_spec.id, linker_spec.id),
            reaction_ids=(template.id,),
            metadata={
                "planning_mode": "batch-indexed-node-linker",
                "node_monomer_id": node_spec.id,
                "linker_monomer_id": linker_spec.id,
                "connectivities": (len(node_spec.motifs), len(linker_spec.motifs)),
                "topology_metadata": dict(topology.metadata),
                "topology_family": f"indexed-{topology.dimensionality.lower()}",
                "expanded_node_count": len(expanded.node_sites),
                "expanded_edge_count": len(expanded.edge_sites),
            },
        )
        assignment_plan = AssignmentPlan(
            net_plan=net_plan,
            slot_to_monomer=slot_to_monomer,
            slot_to_conformer=slot_to_conformer,
            metadata={"assignment_mode": "batch-indexed-node-linker"},
        )

        events: list[ReactionEvent] = []
        for edge in expanded.edge_sites:
            direction = self._safe_normalize(edge.base_vector)
            linker_layout = self._linker_layout_for_direction(linker_spec, direction, template_id=effective_template_id)
            linker_instance_id = linker_instance_ids[edge.id]
            start_motif_id = str(
                node_placements[edge.start_node_id].motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")]
            )
            end_motif_id = str(
                node_placements[edge.end_node_id].motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")]
            )
            linker_near_start, linker_near_end = linker_layout.motif_ids
            events.append(
                ReactionEvent(
                    id=f"rxn{len(events) + 1}",
                    template_id=template.id,
                    participants=(
                        MotifRef(
                            monomer_instance_id=edge.start_node_id,
                            monomer_id=node_spec.id,
                            motif_id=start_motif_id,
                        ),
                        MotifRef(
                            monomer_instance_id=linker_instance_id,
                            monomer_id=linker_spec.id,
                            motif_id=linker_near_start,
                        ),
                    ),
                    metadata={"topology_edge_id": edge.id, "endpoint": "start"},
                )
            )
            events.append(
                ReactionEvent(
                    id=f"rxn{len(events) + 1}",
                    template_id=template.id,
                    participants=(
                        MotifRef(
                            monomer_instance_id=edge.end_node_id,
                            monomer_id=node_spec.id,
                            motif_id=end_motif_id,
                            periodic_image=edge.end_image,
                        ),
                        MotifRef(
                            monomer_instance_id=linker_instance_id,
                            monomer_id=linker_spec.id,
                            motif_id=linker_near_end,
                        ),
                    ),
                    metadata={"topology_edge_id": edge.id, "endpoint": "end", "node_image": edge.end_image},
                )
            )

        return AssignmentOutcome(
            assignment_plan=assignment_plan,
            monomer_instances=tuple(instances),
            events=tuple(events),
            unreacted_motifs=(),
            consumed_count=sum(len(event.participants) for event in events),
        )

    def _embed_indexed_node_linker(
        self,
        outcome: AssignmentOutcome,
        node_spec: MonomerSpec,
        linker_spec: MonomerSpec,
        template,
        topology,
        expanded: ExpandedIndexedTopology,
    ) -> EmbeddingResult:
        effective_template_id = template.id
        bridge_target = self._template_target_distance(template)
        node_placements = self._expanded_node_site_placements(
            expanded.node_sites,
            {node_site.id: node_spec for node_site in expanded.node_sites},
            dimensionality=topology.dimensionality,
            template_id=effective_template_id,
        )
        node_motifs = {motif.id: motif for motif in node_spec.motifs}
        linker_layouts: dict[str, _LinkerLayout] = {}

        target_edge_vectors: dict[str, Vec3] = {}
        reactive_site_distances: list[float] = []
        for edge in expanded.edge_sites:
            direction = self._safe_normalize(edge.base_vector)
            start_spec = node_spec
            end_spec = node_spec
            start_placement = node_placements[edge.start_node_id]
            end_placement = node_placements[edge.end_node_id]
            start_motif = node_motifs[
                str(start_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")])
            ]
            end_motif = node_motifs[
                str(end_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")])
            ]
            if topology.dimensionality == "3D":
                start_normal = matmul_vec(start_placement.rotation, start_motif.frame.normal)
                end_normal = matmul_vec(end_placement.rotation, end_motif.frame.normal)
                target_normal = self._bridge_plane_normal(direction, start_normal, end_normal)
                linker_layout = self._linker_layout_for_axes(
                    linker_spec,
                    direction,
                    target_normal,
                    template_id=effective_template_id,
                )
            else:
                linker_layout = self._linker_layout_for_direction(linker_spec, direction, template_id=effective_template_id)
            linker_layouts[edge.id] = linker_layout
            reactive_site_span = bridge_target + linker_layout.span + bridge_target
            edge_vector = add(
                scale(direction, reactive_site_span),
                sub(
                    matmul_vec(
                        start_placement.rotation,
                        effective_motif_origin(effective_template_id, start_spec, start_motif),
                    ),
                    matmul_vec(
                        end_placement.rotation,
                        effective_motif_origin(effective_template_id, end_spec, end_motif),
                    ),
                ),
            )
            target_edge_vectors[edge.id] = edge_vector
            reactive_site_distances.append(norm(edge_vector))

        cell = self._fit_expanded_single_node_cell(
            expanded,
            target_edge_vectors,
            dimensionality=topology.dimensionality,
        )
        monomer_poses: dict[str, Pose] = {}
        pose_details: dict[str, object] = {}
        for node_site in expanded.node_sites:
            placement = node_placements[node_site.id]
            translation = self._fractional_position_in_cell(cell, node_site.fractional_position)
            monomer_poses[node_site.id] = Pose(translation=translation, rotation_matrix=placement.rotation)
            pose_details[node_site.id] = {
                "slot_id": f"node_{node_site.id}",
                "translation": translation,
                "rotation_matrix": placement.rotation,
                "motif_count": len(node_spec.motifs),
                "role": "node",
                "topology_node_id": node_site.id,
                "source_node_index": node_site.source_node_index,
                "source_label": node_site.source_label,
                "connectivity": node_site.connectivity,
                "fractional_position": node_site.fractional_position,
                "radial_offsets": tuple(round(value, 6) for value in placement.offsets),
            }

        for edge in expanded.edge_sites:
            linker_layout = linker_layouts[edge.id]
            linker_instance_id = next(
                instance.id
                for instance in outcome.monomer_instances
                if instance.metadata.get("topology_edge_id") == edge.id
            )
            start_spec = node_spec
            end_spec = node_spec
            start_placement = node_placements[edge.start_node_id]
            end_placement = node_placements[edge.end_node_id]
            start_motif = node_motifs[
                str(start_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")])
            ]
            end_motif = node_motifs[
                str(end_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")])
            ]
            start_origin = self._world_motif_origin(
                cell,
                monomer_poses[edge.start_node_id],
                start_spec,
                start_motif,
                (0, 0, 0),
                template_id=effective_template_id,
            )
            end_origin = self._world_motif_origin(
                cell,
                monomer_poses[edge.end_node_id],
                end_spec,
                end_motif,
                edge.end_image,
                template_id=effective_template_id,
            )
            bridge_direction = self._safe_normalize(sub(end_origin, start_origin))
            desired_start = add(start_origin, scale(bridge_direction, bridge_target))
            desired_end = add(end_origin, scale(bridge_direction, -bridge_target))
            linker_local_a, linker_local_b = linker_layout.local_origins
            translation_a = sub(desired_start, matmul_vec(linker_layout.rotation, linker_local_a))
            translation_b = sub(desired_end, matmul_vec(linker_layout.rotation, linker_local_b))
            linker_translation = scale(add(translation_a, translation_b), 0.5)
            monomer_poses[linker_instance_id] = Pose(
                translation=linker_translation,
                rotation_matrix=linker_layout.rotation,
            )
            pose_details[linker_instance_id] = {
                "slot_id": f"linker_{edge.id}",
                "translation": linker_translation,
                "rotation_matrix": linker_layout.rotation,
                "motif_count": len(linker_spec.motifs),
                "role": "linker",
                "topology_edge_id": edge.id,
                "node_image": edge.end_image,
                "edge_direction": tuple(round(value, 6) for value in self._safe_normalize(edge.base_vector)),
                "layout_span": linker_layout.span,
                "layout_offsets": (round(linker_layout.offset_a, 6), round(linker_layout.offset_b, 6)),
                "translation_mismatch": round(norm(sub(translation_b, translation_a)), 6),
            }

        return EmbeddingResult(
            state=AssemblyState(cell=cell, monomer_poses=monomer_poses, stacking_state="disabled"),
            metadata={
                "mode": "topology-guided",
                "topology": topology.id,
                "target_distance": bridge_target,
                "reactive_site_distance": sum(reactive_site_distances) / len(reactive_site_distances),
                "edge_reactive_site_distances": tuple(round(value, 6) for value in reactive_site_distances),
                "cell_kind": self._topology_cell_kind(cell, topology.dimensionality),
                "stacking_enabled": False,
                "placement_mode": (
                    "node-linker-indexed-topology-3d"
                    if topology.dimensionality == "3D"
                    else "node-linker-indexed-topology"
                ),
                "node_instance_count": len(expanded.node_sites),
                "linker_instance_count": len(expanded.edge_sites),
                "topology_family": f"indexed-{topology.dimensionality.lower()}",
                "expanded_node_orbit_size": len(expanded.node_sites),
                "expanded_edge_count": len(expanded.edge_sites),
                "poses": pose_details,
            },
        )

    def _build_indexed_node_node_outcome(
        self,
        first_spec: MonomerSpec,
        second_spec: MonomerSpec,
        template,
        topology,
        expanded: ExpandedIndexedTopology,
        spec_by_node_id: Mapping[str, MonomerSpec],
    ) -> AssignmentOutcome:
        effective_template_id = template.id
        placements = self._expanded_node_site_placements(
            expanded.node_sites,
            spec_by_node_id,
            dimensionality=topology.dimensionality,
            template_id=effective_template_id,
        )

        slot_to_monomer: dict[str, str] = {}
        slot_to_conformer: dict[str, str] = {}
        instances: list[MonomerInstance] = []
        for node_site in expanded.node_sites:
            spec = spec_by_node_id[node_site.id]
            slot_id = f"node_{node_site.id}"
            slot_to_monomer[slot_id] = spec.id
            if spec.conformer_ids:
                slot_to_conformer[slot_id] = spec.conformer_ids[0]
            instances.append(
                MonomerInstance(
                    id=node_site.id,
                    monomer_id=spec.id,
                    conformer_id=spec.conformer_ids[0] if spec.conformer_ids else None,
                    metadata={
                        "slot_id": slot_id,
                        "role": "node",
                        "topology_node_id": node_site.id,
                        "source_node_index": node_site.source_node_index,
                        "source_label": node_site.source_label,
                        "connectivity": node_site.connectivity,
                        "fractional_position": node_site.fractional_position,
                    },
                )
            )

        net_plan = NetPlan(
            topology=topology,
            monomer_ids=(first_spec.id, second_spec.id),
            reaction_ids=(template.id,),
            metadata={
                "planning_mode": "batch-indexed-node-node",
                "connectivities": (len(first_spec.motifs), len(second_spec.motifs)),
                "topology_metadata": dict(topology.metadata),
                "topology_family": f"indexed-{topology.dimensionality.lower()}",
                "expanded_node_count": len(expanded.node_sites),
                "expanded_edge_count": len(expanded.edge_sites),
            },
        )
        assignment_plan = AssignmentPlan(
            net_plan=net_plan,
            slot_to_monomer=slot_to_monomer,
            slot_to_conformer=slot_to_conformer,
            metadata={"assignment_mode": "batch-indexed-node-node"},
        )

        events: list[ReactionEvent] = []
        for edge in expanded.edge_sites:
            start_spec = spec_by_node_id[edge.start_node_id]
            end_spec = spec_by_node_id[edge.end_node_id]
            start_motif_id = str(
                placements[edge.start_node_id].motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")]
            )
            end_motif_id = str(
                placements[edge.end_node_id].motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")]
            )
            events.append(
                ReactionEvent(
                    id=f"rxn{len(events) + 1}",
                    template_id=template.id,
                    participants=(
                        MotifRef(
                            monomer_instance_id=edge.start_node_id,
                            monomer_id=start_spec.id,
                            motif_id=start_motif_id,
                        ),
                        MotifRef(
                            monomer_instance_id=edge.end_node_id,
                            monomer_id=end_spec.id,
                            motif_id=end_motif_id,
                            periodic_image=edge.end_image,
                        ),
                    ),
                    metadata={"topology_edge_id": edge.id, "node_image": edge.end_image},
                )
            )

        return AssignmentOutcome(
            assignment_plan=assignment_plan,
            monomer_instances=tuple(instances),
            events=tuple(events),
            unreacted_motifs=(),
            consumed_count=sum(len(event.participants) for event in events),
        )

    def _embed_indexed_node_node(
        self,
        outcome: AssignmentOutcome,
        first_spec: MonomerSpec,
        second_spec: MonomerSpec,
        template,
        topology,
        expanded: ExpandedIndexedTopology,
        spec_by_node_id: Mapping[str, MonomerSpec],
    ) -> EmbeddingResult:
        effective_template_id = template.id
        bridge_target = self._template_target_distance(template)
        placements = self._expanded_node_site_placements(
            expanded.node_sites,
            spec_by_node_id,
            dimensionality=topology.dimensionality,
            template_id=effective_template_id,
        )
        motifs_by_spec = {
            first_spec.id: {motif.id: motif for motif in first_spec.motifs},
            second_spec.id: {motif.id: motif for motif in second_spec.motifs},
        }

        target_edge_vectors: dict[str, Vec3] = {}
        reactive_site_distances: list[float] = []
        for edge in expanded.edge_sites:
            direction = self._safe_normalize(edge.base_vector)
            start_spec = spec_by_node_id[edge.start_node_id]
            end_spec = spec_by_node_id[edge.end_node_id]
            start_placement = placements[edge.start_node_id]
            end_placement = placements[edge.end_node_id]
            start_motif = motifs_by_spec[start_spec.id][
                str(start_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "start")])
            ]
            end_motif = motifs_by_spec[end_spec.id][
                str(end_placement.motif_id_by_edge[self._edge_endpoint_key(edge.id, "end")])
            ]
            edge_vector = add(
                scale(direction, bridge_target),
                sub(
                    matmul_vec(
                        start_placement.rotation,
                        effective_motif_origin(effective_template_id, start_spec, start_motif),
                    ),
                    matmul_vec(
                        end_placement.rotation,
                        effective_motif_origin(effective_template_id, end_spec, end_motif),
                    ),
                ),
            )
            target_edge_vectors[edge.id] = edge_vector
            reactive_site_distances.append(norm(edge_vector))

        cell = self._fit_expanded_single_node_cell(
            expanded,
            target_edge_vectors,
            dimensionality=topology.dimensionality,
        )
        monomer_poses: dict[str, Pose] = {}
        pose_details: dict[str, object] = {}
        for node_site in expanded.node_sites:
            spec = spec_by_node_id[node_site.id]
            placement = placements[node_site.id]
            translation = self._fractional_position_in_cell(cell, node_site.fractional_position)
            monomer_poses[node_site.id] = Pose(translation=translation, rotation_matrix=placement.rotation)
            pose_details[node_site.id] = {
                "slot_id": f"node_{node_site.id}",
                "translation": translation,
                "rotation_matrix": placement.rotation,
                "motif_count": len(spec.motifs),
                "role": "node",
                "topology_node_id": node_site.id,
                "source_node_index": node_site.source_node_index,
                "source_label": node_site.source_label,
                "connectivity": node_site.connectivity,
                "fractional_position": node_site.fractional_position,
                "radial_offsets": tuple(round(value, 6) for value in placement.offsets),
            }

        return EmbeddingResult(
            state=AssemblyState(cell=cell, monomer_poses=monomer_poses, stacking_state="disabled"),
            metadata={
                "mode": "topology-guided",
                "topology": topology.id,
                "target_distance": bridge_target,
                "reactive_site_distance": sum(reactive_site_distances) / len(reactive_site_distances),
                "edge_reactive_site_distances": tuple(round(value, 6) for value in reactive_site_distances),
                "cell_kind": self._topology_cell_kind(cell, topology.dimensionality),
                "stacking_enabled": False,
                "placement_mode": (
                    "indexed-topology-node-node-3d"
                    if topology.dimensionality == "3D"
                    else "indexed-topology-node-node"
                ),
                "node_instance_count": len(expanded.node_sites),
                "topology_family": f"indexed-{topology.dimensionality.lower()}",
                "expanded_node_orbit_size": len(expanded.node_sites),
                "expanded_edge_count": len(expanded.edge_sites),
                "poses": pose_details,
            },
        )

    def _expanded_node_site_placements(
        self,
        node_sites: tuple[ExpandedSingleNodeNodeSite | ExpandedIndexedNodeSite, ...],
        spec_by_node_id: Mapping[str, MonomerSpec],
        *,
        dimensionality: str,
        template_id: str | None = None,
    ) -> dict[str, _NodeSitePlacement]:
        placements: dict[str, _NodeSitePlacement] = {}
        for node_site in node_sites:
            spec = spec_by_node_id[node_site.id]
            if dimensionality == "3D":
                rotation, motif_ids, offsets = self._rotation_for_spatial_motifs(
                    spec,
                    node_site.directions,
                    template_id=template_id,
                )
            else:
                rotation, motif_ids, offsets = self._rotation_for_planar_motifs(
                    spec,
                    node_site.directions,
                    template_id=template_id,
                )
            motif_id_by_edge = {
                edge_id: motif_id for edge_id, motif_id in zip(node_site.edge_ids, motif_ids)
            }
            placements[node_site.id] = _NodeSitePlacement(
                rotation=rotation,
                motif_ids=motif_ids,
                offsets=offsets,
                motif_id_by_edge=motif_id_by_edge,
            )
        return placements

    def _build_two_node_single_node_node_node_candidate(
        self,
        amine: MonomerSpec,
        aldehyde: MonomerSpec,
        template,
        project: COFProject,
    ) -> Candidate:
        monomer_specs = {amine.id: amine, aldehyde.id: aldehyde}
        templates = (template,)
        template_map = {template.id: template}
        net_plan = self.net_planner.propose(
            monomers=(amine, aldehyde),
            templates=templates,
            target_dimensionality=project.target_dimensionality,
            target_topologies=project.target_topologies,
        )[0]
        assignment_plan = self.assignment_solver.build_assignment_plans(net_plan, monomer_specs)[0]
        monomer_instances = self.assignment_solver.instantiate_monomers(assignment_plan)
        outcome = self.assignment_solver.solve_events(
            assignment_plan=assignment_plan,
            monomer_instances=monomer_instances,
            monomer_specs=monomer_specs,
            templates=templates,
        )
        graph = PeriodicProductGraph()
        for instance in monomer_instances:
            graph.add_monomer(instance)
        for event in outcome.events:
            graph.add_reaction_event(event, monomer_specs, template_map)
        embedding = self.embedder.embed(outcome, monomer_specs, template_map)
        return self._assemble_candidate(
            candidate_id=f"{amine.id}__{aldehyde.id}__{project.target_topologies[0]}__single_node_node_node",
            project=project,
            outcome=outcome,
            graph=graph,
            embedding=embedding,
            monomer_specs=monomer_specs,
            templates=template_map,
            extra_flags=("batch_single_node_node_node",),
        )

    def _fit_expanded_single_node_cell(
        self,
        expanded: ExpandedSingleNodeTopology | ExpandedIndexedTopology,
        target_edge_vectors: Mapping[str, Vec3],
        *,
        dimensionality: str,
    ) -> tuple[Vec3, Vec3, Vec3]:
        if dimensionality == "3D":
            return self._fit_expanded_single_node_cell_3d(expanded, target_edge_vectors)

        s_xx = 0.0
        s_xy = 0.0
        s_yy = 0.0
        rhs_x_0 = 0.0
        rhs_x_1 = 0.0
        rhs_y_0 = 0.0
        rhs_y_1 = 0.0

        node_positions = {node_site.id: node_site.fractional_position for node_site in expanded.node_sites}
        for edge in expanded.edge_sites:
            delta = (
                node_positions[edge.end_node_id][0] + edge.end_image[0] - node_positions[edge.start_node_id][0],
                node_positions[edge.end_node_id][1] + edge.end_image[1] - node_positions[edge.start_node_id][1],
            )
            target = target_edge_vectors[edge.id]
            s_xx += delta[0] * delta[0]
            s_xy += delta[0] * delta[1]
            s_yy += delta[1] * delta[1]
            rhs_x_0 += delta[0] * target[0]
            rhs_x_1 += delta[1] * target[0]
            rhs_y_0 += delta[0] * target[1]
            rhs_y_1 += delta[1] * target[1]

        determinant = s_xx * s_yy - s_xy * s_xy
        if abs(determinant) < 1e-8:
            return self._fallback_expanded_single_node_cell(expanded, target_edge_vectors)

        cell_a = (
            (rhs_x_0 * s_yy - rhs_x_1 * s_xy) / determinant,
            (rhs_y_0 * s_yy - rhs_y_1 * s_xy) / determinant,
            0.0,
        )
        cell_b = (
            (s_xx * rhs_x_1 - s_xy * rhs_x_0) / determinant,
            (s_xx * rhs_y_1 - s_xy * rhs_y_0) / determinant,
            0.0,
        )
        if norm(cell_a) < 1e-8 or norm(cell_b) < 1e-8:
            return self._fallback_expanded_single_node_cell(expanded, target_edge_vectors)
        return (
            cell_a,
            cell_b,
            (0.0, 0.0, self.config.embedding_config.default_layer_spacing),
        )

    def _fit_expanded_single_node_cell_3d(
        self,
        expanded: ExpandedSingleNodeTopology | ExpandedIndexedTopology,
        target_edge_vectors: Mapping[str, Vec3],
    ) -> tuple[Vec3, Vec3, Vec3]:
        node_positions = {node_site.id: node_site.fractional_position for node_site in expanded.node_sites}
        normal_matrix = [[0.0, 0.0, 0.0] for _ in range(3)]
        rhs_x = [0.0, 0.0, 0.0]
        rhs_y = [0.0, 0.0, 0.0]
        rhs_z = [0.0, 0.0, 0.0]

        for edge in expanded.edge_sites:
            delta = (
                node_positions[edge.end_node_id][0] + edge.end_image[0] - node_positions[edge.start_node_id][0],
                node_positions[edge.end_node_id][1] + edge.end_image[1] - node_positions[edge.start_node_id][1],
                node_positions[edge.end_node_id][2] + edge.end_image[2] - node_positions[edge.start_node_id][2],
            )
            target = target_edge_vectors[edge.id]
            for row in range(3):
                rhs_x[row] += delta[row] * target[0]
                rhs_y[row] += delta[row] * target[1]
                rhs_z[row] += delta[row] * target[2]
                for col in range(3):
                    normal_matrix[row][col] += delta[row] * delta[col]

        solution_x = self._solve_symmetric_3x3(normal_matrix, rhs_x)
        solution_y = self._solve_symmetric_3x3(normal_matrix, rhs_y)
        solution_z = self._solve_symmetric_3x3(normal_matrix, rhs_z)
        if solution_x is None or solution_y is None or solution_z is None:
            return self._fallback_expanded_single_node_cell_3d(expanded, target_edge_vectors)

        cell = (
            (solution_x[0], solution_y[0], solution_z[0]),
            (solution_x[1], solution_y[1], solution_z[1]),
            (solution_x[2], solution_y[2], solution_z[2]),
        )
        if any(norm(vector) < 1e-8 for vector in cell):
            return self._fallback_expanded_single_node_cell_3d(expanded, target_edge_vectors)
        return cell

    def _fallback_expanded_single_node_cell(
        self,
        expanded: ExpandedSingleNodeTopology | ExpandedIndexedTopology,
        target_edge_vectors: Mapping[str, Vec3],
    ) -> tuple[Vec3, Vec3, Vec3]:
        base_lengths = [norm(edge.base_vector) for edge in expanded.edge_sites if norm(edge.base_vector) > 1e-8]
        target_lengths = [norm(target_edge_vectors[edge.id]) for edge in expanded.edge_sites if edge.id in target_edge_vectors]
        scale_factor = 1.0
        if base_lengths and target_lengths:
            scale_factor = sum(target_lengths) / sum(base_lengths)
        return (
            scale(expanded.cell[0], scale_factor),
            scale(expanded.cell[1], scale_factor),
            (0.0, 0.0, self.config.embedding_config.default_layer_spacing),
        )

    def _fallback_expanded_single_node_cell_3d(
        self,
        expanded: ExpandedSingleNodeTopology | ExpandedIndexedTopology,
        target_edge_vectors: Mapping[str, Vec3],
    ) -> tuple[Vec3, Vec3, Vec3]:
        base_lengths = [norm(edge.base_vector) for edge in expanded.edge_sites if norm(edge.base_vector) > 1e-8]
        target_lengths = [norm(target_edge_vectors[edge.id]) for edge in expanded.edge_sites if edge.id in target_edge_vectors]
        scale_factor = 1.0
        if base_lengths and target_lengths:
            scale_factor = sum(target_lengths) / sum(base_lengths)
        return tuple(scale(vector, scale_factor) for vector in expanded.cell)  # type: ignore[return-value]

    def _fractional_position_in_cell(
        self,
        cell: tuple[Vec3, Vec3, Vec3],
        fractional_position: tuple[float, float, float],
    ) -> Vec3:
        return add(
            add(scale(cell[0], fractional_position[0]), scale(cell[1], fractional_position[1])),
            scale(cell[2], fractional_position[2]),
        )

    def _node_node_supercell_topology(
        self,
        expanded: ExpandedSingleNodeTopology,
    ) -> ExpandedSingleNodeTopology:
        if not expanded.is_bipartite:
            raise ValueError(
                f"topology {expanded.topology_id!r} is not bipartite and cannot host a single-node node-node alternation"
            )

        colors, center = self._expanded_topology_color_map(expanded)
        scale_x, scale_y = self._expanded_topology_color_periods(expanded, colors, center)
        if scale_x == 1 and scale_y == 1:
            return expanded

        supercell_cell = (
            scale(expanded.cell[0], scale_x),
            scale(expanded.cell[1], scale_y),
            expanded.cell[2],
        )
        base_node_by_id = {node_site.id: node_site for node_site in expanded.node_sites}
        node_sites: list[ExpandedSingleNodeNodeSite] = []
        node_id_map: dict[tuple[str, int, int], str] = {}
        for tile_x in range(scale_x):
            for tile_y in range(scale_y):
                for base_node in expanded.node_sites:
                    node_id = f"{base_node.id}_x{tile_x}_y{tile_y}"
                    node_id_map[(base_node.id, tile_x, tile_y)] = node_id
                    fractional_position = (
                        (base_node.fractional_position[0] + tile_x) / scale_x,
                        (base_node.fractional_position[1] + tile_y) / scale_y,
                        base_node.fractional_position[2],
                    )
                    node_sites.append(
                        ExpandedSingleNodeNodeSite(
                            id=node_id,
                            fractional_position=fractional_position,
                            cartesian_position=self._fractional_position_in_cell(supercell_cell, fractional_position),
                            directions=base_node.directions,
                            edge_ids=(),
                            sublattice=colors[(base_node.id, center + tile_x, center + tile_y)],
                        )
                    )

        edge_sites: list[ExpandedSingleNodeEdge] = []
        for tile_x in range(scale_x):
            for tile_y in range(scale_y):
                for base_edge in expanded.edge_sites:
                    raw_end_x = tile_x + base_edge.end_image[0]
                    raw_end_y = tile_y + base_edge.end_image[1]
                    wrapped_end_x = raw_end_x % scale_x
                    wrapped_end_y = raw_end_y % scale_y
                    end_image = (
                        (raw_end_x - wrapped_end_x) // scale_x,
                        (raw_end_y - wrapped_end_y) // scale_y,
                        0,
                    )
                    start_node_id = node_id_map[(base_edge.start_node_id, tile_x, tile_y)]
                    end_node_id = node_id_map[(base_edge.end_node_id, wrapped_end_x, wrapped_end_y)]
                    start_fractional = next(
                        node_site.fractional_position
                        for node_site in node_sites
                        if node_site.id == start_node_id
                    )
                    end_fractional = next(
                        node_site.fractional_position
                        for node_site in node_sites
                        if node_site.id == end_node_id
                    )
                    center_fractional = (
                        (start_fractional[0] + end_fractional[0] + end_image[0]) * 0.5 % 1.0,
                        (start_fractional[1] + end_fractional[1] + end_image[1]) * 0.5 % 1.0,
                        0.0,
                    )
                    edge_sites.append(
                        ExpandedSingleNodeEdge(
                            id=f"{base_edge.id}_x{tile_x}_y{tile_y}",
                            start_node_id=start_node_id,
                            end_node_id=end_node_id,
                            end_image=end_image,
                            base_vector=base_edge.base_vector,
                            center_fractional=center_fractional,
                        )
                    )

        incident_by_node_id: dict[str, list[tuple[float, Vec3, str]]] = {
            node_site.id: [] for node_site in node_sites
        }
        for edge in edge_sites:
            incident_by_node_id[edge.start_node_id].append(
                (
                    atan2(edge.base_vector[1], edge.base_vector[0]),
                    self._safe_normalize(edge.base_vector),
                    self._edge_endpoint_key(edge.id, "start"),
                )
            )
            reverse_vector = scale(edge.base_vector, -1.0)
            incident_by_node_id[edge.end_node_id].append(
                (
                    atan2(reverse_vector[1], reverse_vector[0]),
                    self._safe_normalize(reverse_vector),
                    self._edge_endpoint_key(edge.id, "end"),
                )
            )

        supercell_nodes = tuple(
            ExpandedSingleNodeNodeSite(
                id=node_site.id,
                fractional_position=node_site.fractional_position,
                cartesian_position=node_site.cartesian_position,
                directions=tuple(direction for _, direction, _ in sorted(incident_by_node_id[node_site.id], key=lambda item: item[0])),
                edge_ids=tuple(edge_id for _, _, edge_id in sorted(incident_by_node_id[node_site.id], key=lambda item: item[0])),
                sublattice=node_site.sublattice,
            )
            for node_site in node_sites
        )

        return ExpandedSingleNodeTopology(
            topology_id=expanded.topology_id,
            connectivity=expanded.connectivity,
            cell=supercell_cell,
            metric_family=expanded.metric_family,
            node_sites=supercell_nodes,
            edge_sites=tuple(edge_sites),
            is_bipartite=True,
        )

    def _expanded_topology_color_map(
        self,
        expanded: ExpandedSingleNodeTopology,
    ) -> tuple[dict[tuple[str, int, int], int], int]:
        max_shift = max(
            max(abs(edge.end_image[0]), abs(edge.end_image[1]))
            for edge in expanded.edge_sites
        )
        reps = 2 * max_shift + 5
        center = max_shift + 2
        colors: dict[tuple[str, int, int], int] = {}
        for node_site in expanded.node_sites:
            seed = (node_site.id, center, center)
            if seed in colors:
                continue
            colors[seed] = 0
            queue = [seed]
            while queue:
                node_id, cell_x, cell_y = queue.pop(0)
                current_color = colors[(node_id, cell_x, cell_y)]
                for edge in expanded.edge_sites:
                    neighbors: list[tuple[str, int, int]] = []
                    if edge.start_node_id == node_id:
                        neighbors.append((edge.end_node_id, cell_x + edge.end_image[0], cell_y + edge.end_image[1]))
                    if edge.end_node_id == node_id:
                        neighbors.append((edge.start_node_id, cell_x - edge.end_image[0], cell_y - edge.end_image[1]))
                    for other in neighbors:
                        if not (0 <= other[1] < reps and 0 <= other[2] < reps):
                            continue
                        if other not in colors:
                            colors[other] = 1 - current_color
                            queue.append(other)
                            continue
                        if colors[other] == current_color:
                            raise ValueError(
                                f"topology {expanded.topology_id!r} is not bipartite and cannot host a single-node node-node alternation"
                            )
        return colors, center

    def _expanded_topology_color_periods(
        self,
        expanded: ExpandedSingleNodeTopology,
        colors: Mapping[tuple[str, int, int], int],
        center: int,
    ) -> tuple[int, int]:
        for scale_x in range(1, 5):
            for scale_y in range(1, 5):
                periodic = True
                for node_site in expanded.node_sites:
                    for delta_x in range(2):
                        for delta_y in range(2):
                            base_key = (node_site.id, center + delta_x, center + delta_y)
                            x_key = (node_site.id, center + delta_x + scale_x, center + delta_y)
                            y_key = (node_site.id, center + delta_x, center + delta_y + scale_y)
                            if colors.get(base_key) != colors.get(x_key) or colors.get(base_key) != colors.get(y_key):
                                periodic = False
                                break
                        if not periodic:
                            break
                    if not periodic:
                        break
                if periodic:
                    return scale_x, scale_y
        raise ValueError(
            f"topology {expanded.topology_id!r} requires a larger color supercell than the current batch builder supports"
        )

    def _assemble_candidate(
        self,
        *,
        candidate_id: str,
        project: COFProject,
        outcome: AssignmentOutcome,
        graph: PeriodicProductGraph,
        embedding: EmbeddingResult,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, object],
        extra_flags: tuple[str, ...] = (),
    ) -> Candidate:
        optimization = self.optimizer.optimize(outcome, embedding.state, monomer_specs, templates)
        scoring = self.scorer.score(outcome, optimization.state, monomer_specs, templates)

        flags: list[str] = list(extra_flags)
        if outcome.unreacted_motifs:
            flags.append(f"unreacted_motifs:{len(outcome.unreacted_motifs)}")
        if outcome.assignment_plan.net_plan.topology is None:
            flags.append("no_topology_hint")
        if project.stacking_mode == "disabled":
            flags.append("stacking_disabled")

        candidate = Candidate(
            id=candidate_id,
            score=scoring.total,
            state=optimization.state,
            events=outcome.events,
            flags=tuple(flags),
            metadata={
                "graph_summary": graph.summary(),
                "target_topologies": project.target_topologies,
                "allowed_reactions": project.allowed_reactions,
                "stacking_mode": project.stacking_mode,
                "net_plan": {
                    "topology": outcome.assignment_plan.net_plan.topology.id
                    if outcome.assignment_plan.net_plan.topology is not None
                    else None,
                    "metadata": dict(outcome.assignment_plan.net_plan.metadata),
                },
                "assignment": dict(outcome.assignment_plan.slot_to_monomer),
                "instance_to_monomer": {
                    instance.id: instance.monomer_id for instance in outcome.monomer_instances
                },
                "instance_to_slot": {
                    instance.id: instance.metadata.get("slot_id") for instance in outcome.monomer_instances
                },
                "embedding": dict(embedding.metadata),
                "optimization": dict(optimization.metrics),
                "score_breakdown": dict(scoring.breakdown),
                "score_metadata": dict(scoring.metadata),
            },
        )
        return annotate_post_build_conversions(
            candidate,
            monomer_specs,
            self.config.post_build_conversions,
        )

    def _resolve_supported_single_node_layout(self, topology_id: str):
        try:
            return resolve_single_node_topology_layout(topology_id)
        except (KeyError, ValueError):
            return resolve_three_d_single_node_topology_layout(topology_id)

    def _expand_supported_single_node_topology(self, topology_id: str) -> ExpandedSingleNodeTopology:
        try:
            return expand_single_node_topology(topology_id)
        except (KeyError, ValueError):
            return expand_three_d_single_node_topology(topology_id)

    def _configured_topology_ids(self) -> tuple[str, ...]:
        return self.config.topology_ids or self.config.single_node_topology_ids

    def _topology_mode_token(self, connectivities: tuple[int, int], pair_mode: str) -> str:
        first_connectivity, second_connectivity = connectivities
        if pair_mode == "node-linker":
            return f"{max(first_connectivity, second_connectivity)}+2"
        low, high = sorted((first_connectivity, second_connectivity))
        return f"{low}+{high}"

    def _topology_mode_metadata_key(self, pair_mode: str) -> str:
        if pair_mode == "node-node":
            return "two_monomer_node_node_modes"
        return "two_monomer_node_linker_modes"

    def _topology_ids_for_pair(
        self,
        *,
        connectivities: tuple[int, int],
        pair_mode: str,
    ) -> tuple[str, ...]:
        key = (tuple(sorted(connectivities)), pair_mode)
        with self._topology_id_cache_lock:
            cached = self._topology_id_cache.get(key)
        if cached is not None:
            return cached

        configured = self._configured_topology_ids()
        if configured:
            candidate_ids = configured
        else:
            legacy_ids = self._legacy_default_topology_ids(connectivities=connectivities, pair_mode=pair_mode)
            curated_ids = self._curated_default_topology_ids(connectivities=connectivities, pair_mode=pair_mode)
            candidate_ids = tuple(dict.fromkeys(legacy_ids + curated_ids))
            if self.config.use_indexed_topology_defaults or not legacy_ids:
                indexed_ids = self._indexed_topology_ids(connectivities=connectivities, pair_mode=pair_mode)
                candidate_ids = tuple(dict.fromkeys(candidate_ids + indexed_ids))

        mode_token = self._topology_mode_token(connectivities, pair_mode)
        mode_key = self._topology_mode_metadata_key(pair_mode)
        selected: list[str] = []
        for topology_id in candidate_ids:
            try:
                hint = self._topology_repository.get_hint(topology_id)
            except KeyError:
                continue
            topology_modes = tuple(str(value) for value in hint.metadata.get(mode_key, ()))
            if mode_token not in topology_modes:
                continue
            builder_error = self._topology_builder_compatibility_error(
                topology_id=topology_id,
                connectivities=connectivities,
                pair_mode=pair_mode,
            )
            if builder_error is not None:
                continue
            if topology_id not in selected:
                selected.append(topology_id)

        result = tuple(selected)
        with self._topology_id_cache_lock:
            existing = self._topology_id_cache.get(key)
            if existing is not None:
                return existing
            self._topology_id_cache[key] = result
        return result

    def _curated_default_topology_ids(
        self,
        *,
        connectivities: tuple[int, int],
        pair_mode: str,
    ) -> tuple[str, ...]:
        mode_token = self._topology_mode_token(connectivities, pair_mode)
        mode_key = self._topology_mode_metadata_key(pair_mode)
        selected: list[str] = []
        for topology_id in _CURATED_DEFAULT_TOPOLOGY_IDS:
            try:
                hint = self._topology_repository.get_hint(topology_id)
            except KeyError:
                continue
            topology_modes = tuple(str(value) for value in hint.metadata.get(mode_key, ()))
            if mode_token not in topology_modes:
                continue
            builder_error = self._topology_builder_compatibility_error(
                topology_id=topology_id,
                connectivities=connectivities,
                pair_mode=pair_mode,
            )
            if builder_error is not None:
                continue
            selected.append(topology_id)
        return tuple(selected)

    def _indexed_topology_ids(
        self,
        *,
        connectivities: tuple[int, int],
        pair_mode: str,
    ) -> tuple[str, ...]:
        mode_token = self._topology_mode_token(connectivities, pair_mode)
        mode_key = self._topology_mode_metadata_key(pair_mode)
        return tuple(
            entry.id
            for entry in self._topology_repository.list_index()
            if mode_token in tuple(str(value) for value in entry.metadata.get(mode_key, ()))
        )

    def _single_node_topology_ids(self, *, connectivity: int, pair_mode: str) -> tuple[str, ...]:
        if pair_mode == "node-linker":
            connectivities = (connectivity, 2)
        else:
            connectivities = (connectivity, connectivity)
        return tuple(
            topology_id
            for topology_id in self._topology_ids_for_pair(connectivities=connectivities, pair_mode=pair_mode)
            if topology_id in self._legacy_default_topology_ids(connectivities=connectivities, pair_mode=pair_mode)
        )

    def _topology_unavailable_errors(
        self,
        *,
        connectivities: tuple[int, int],
        pair_mode: str,
        generic_message: str,
    ) -> dict[str, str]:
        configured = self._configured_topology_ids()
        if not configured:
            return {"__pair__": generic_message}

        mode_token = self._topology_mode_token(connectivities, pair_mode)
        mode_key = self._topology_mode_metadata_key(pair_mode)
        errors: dict[str, str] = {}
        for topology_id in configured:
            try:
                hint = self._topology_repository.get_hint(topology_id)
            except KeyError as exc:
                errors[topology_id] = f"{type(exc).__name__}: {exc}"
                continue
            topology_modes = tuple(str(value) for value in hint.metadata.get(mode_key, ()))
            if mode_token not in topology_modes:
                reason_key = (
                    "two_monomer_node_node_reason"
                    if pair_mode == "node-node"
                    else "two_monomer_node_linker_reason"
                )
                reason = hint.metadata.get(reason_key) or hint.metadata.get("two_monomer_reason")
                errors[topology_id] = str(reason or f"topology {topology_id!r} is not chemically compatible with mode {mode_token!r}")
                continue
            builder_error = self._topology_builder_compatibility_error(
                topology_id=topology_id,
                connectivities=connectivities,
                pair_mode=pair_mode,
            )
            if builder_error is not None:
                errors[topology_id] = builder_error
                continue
        if not errors:
            errors["__pair__"] = generic_message
        elif "__pair__" not in errors:
            errors["__pair__"] = next(iter(errors.values()))
        return errors

    def _topology_builder_compatibility_error(
        self,
        *,
        topology_id: str,
        connectivities: tuple[int, int],
        pair_mode: str,
    ) -> str | None:
        try:
            layout = self._resolve_supported_single_node_layout(topology_id)
        except (KeyError, ValueError):
            return None

        if not layout.supports_current_builder:
            return layout.unsupported_reason or (
                f"topology {topology_id!r} is not supported by the current single-node builder"
            )

        if pair_mode == "node-linker":
            if 2 not in connectivities:
                return f"topology {topology_id!r} requires a ditopic linker in node-linker mode"
            node_connectivity = max(connectivities)
            if layout.connectivity != node_connectivity:
                return (
                    f"topology {topology_id!r} has connectivity {layout.connectivity}, expected {node_connectivity}"
                )
            if not layout.supports_node_linker:
                return f"topology {topology_id!r} does not support node-linker generation"
            return None

        first_connectivity, second_connectivity = connectivities
        if first_connectivity != second_connectivity:
            return (
                f"topology {topology_id!r} uses a single-node builder and cannot realize mixed-connectivity "
                f"node-node pairing {min(connectivities)}+{max(connectivities)}"
            )
        if layout.connectivity != first_connectivity:
            return f"topology {topology_id!r} has connectivity {layout.connectivity}, expected {first_connectivity}"
        if not layout.supports_node_node:
            return f"topology {topology_id!r} does not support node-node generation"
        return None

    def _legacy_default_topology_ids(
        self,
        *,
        connectivities: tuple[int, int],
        pair_mode: str,
    ) -> tuple[str, ...]:
        first_connectivity, second_connectivity = connectivities
        if pair_mode == "node-node" and first_connectivity == second_connectivity:
            connectivity = first_connectivity
            if connectivity == 3:
                return (self.config.hcb_topology_id,) + tuple(
                    topology_id
                    for topology_id in list_supported_single_node_topology_ids(connectivity, pair_mode=pair_mode)
                    if topology_id != self.config.hcb_topology_id
                )
            if connectivity == 4:
                return ("dia",)
        if pair_mode == "node-linker":
            connectivity = max(first_connectivity, second_connectivity)
            if connectivity == 3:
                return (self.config.hcb_topology_id,) + tuple(
                    topology_id
                    for topology_id in list_supported_single_node_topology_ids(connectivity, pair_mode=pair_mode)
                    if topology_id != self.config.hcb_topology_id
                )
            if connectivity == 4:
                return ("dia",)
            if connectivity == 6:
                return ("pcu",)
        return ()

    def _pair_mode_from_connectivities(self, amine_connectivity: int, aldehyde_connectivity: int) -> str:
        if amine_connectivity >= 3 and aldehyde_connectivity >= 3:
            low, high = sorted((amine_connectivity, aldehyde_connectivity))
            return f"{low}+{high}-node-node"
        if min(amine_connectivity, aldehyde_connectivity) == 2 and max(amine_connectivity, aldehyde_connectivity) >= 3:
            return f"{max(amine_connectivity, aldehyde_connectivity)}+2-node-linker"
        return "unsupported"

    def _batch_pair_iterables(
        self,
        valid: Mapping[str, tuple[BatchMonomerRecord, ...]],
        *,
        template: ReactionTemplate | None = None,
    ) -> tuple[Iterable[tuple[BatchMonomerRecord, BatchMonomerRecord]], ...]:
        selected_template = template or self._selected_binary_bridge_template()
        first_role, second_role = self._binary_bridge_roles(selected_template)
        first_prefix = f"{(first_role.library_prefix or f'{first_role.motif_kind}s')}_count_"
        second_prefix = f"{(second_role.library_prefix or f'{second_role.motif_kind}s')}_count_"
        first_by_connectivity = self._records_by_connectivity(valid, prefix=first_prefix)
        second_by_connectivity = self._records_by_connectivity(valid, prefix=second_prefix)
        pair_iterables: list[Iterable[tuple[BatchMonomerRecord, BatchMonomerRecord]]] = []

        for first_connectivity in sorted(first_by_connectivity):
            if first_connectivity < 3:
                continue
            for second_connectivity in sorted(second_by_connectivity):
                if second_connectivity < 3:
                    continue
                if self._topology_ids_for_pair(
                    connectivities=(first_connectivity, second_connectivity),
                    pair_mode="node-node",
                ):
                    pair_iterables.append(
                        product(first_by_connectivity[first_connectivity], second_by_connectivity[second_connectivity])
                    )

        if 2 in second_by_connectivity:
            for connectivity in sorted(set(first_by_connectivity) - {2}):
                if connectivity < 3:
                    continue
                if self._topology_ids_for_pair(connectivities=(connectivity, 2), pair_mode="node-linker"):
                    pair_iterables.append(product(first_by_connectivity[connectivity], second_by_connectivity[2]))

        if 2 in first_by_connectivity:
            for connectivity in sorted(set(second_by_connectivity) - {2}):
                if connectivity < 3:
                    continue
                if self._topology_ids_for_pair(connectivities=(2, connectivity), pair_mode="node-linker"):
                    pair_iterables.append(product(first_by_connectivity[2], second_by_connectivity[connectivity]))

        return tuple(pair_iterables)

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

    def _single_motif_kind(self, monomer: MonomerSpec) -> str:
        kinds = {motif.kind for motif in monomer.motifs}
        if len(kinds) != 1:
            raise ValueError(
                f"monomer {monomer.id!r} must expose exactly one motif kind for binary bridge generation, got {tuple(sorted(kinds))!r}"
            )
        return next(iter(kinds))

    def _candidate_to_summary(
        self,
        *,
        pair_id: str,
        structure_id: str,
        pair_mode: str,
        reactant_a_id: str,
        reactant_b_id: str,
        reactant_a_connectivity: int,
        reactant_b_connectivity: int,
        candidate: Candidate,
        cif_path: str | None,
        cif_export_index: int,
        topology_rank: int,
        topology_count: int,
        reactant_record_ids: Mapping[str, str],
        reactant_connectivities: Mapping[str, int],
        reactant_roles: tuple[str, str],
        template_id: str,
        hard_hard_invalid_reasons: tuple[str, ...] = (),
        hard_hard_invalid_metrics: Mapping[str, object] | None = None,
    ) -> BatchPairSummary:
        export_blocked = bool(hard_hard_invalid_reasons)
        return BatchPairSummary(
            structure_id=structure_id,
            pair_id=pair_id,
            pair_mode=pair_mode,
            status="ok",
            reactant_a_record_id=reactant_a_id,
            reactant_b_record_id=reactant_b_id,
            reactant_a_connectivity=reactant_a_connectivity,
            reactant_b_connectivity=reactant_b_connectivity,
            topology_id=candidate.metadata["net_plan"]["topology"],
            score=candidate.score,
            flags=tuple(candidate.flags) + tuple(
                f"hard_hard_invalid:{reason}" for reason in hard_hard_invalid_reasons
            ),
            cif_path=cif_path,
            metadata={
                "cif_export_index": cif_export_index,
                "cif_export_blocked": export_blocked,
                "hard_hard_invalid_reasons": hard_hard_invalid_reasons,
                "hard_hard_invalid_metrics": dict(hard_hard_invalid_metrics or {}),
                "graph_summary": dict(candidate.metadata["graph_summary"]),
                "embedding": dict(candidate.metadata["embedding"]),
                "optimization": dict(candidate.metadata["optimization"]),
                "score_breakdown": dict(candidate.metadata["score_breakdown"]),
                "score_metadata": dict(candidate.metadata["score_metadata"]),
                "pair_mode": pair_mode,
                "topology_rank": topology_rank,
                "topology_count": topology_count,
                "topology_selection": dict(candidate.metadata.get("topology_selection", {})),
                "reactant_record_ids": dict(reactant_record_ids),
                "reactant_connectivities": dict(reactant_connectivities),
                "reactant_roles": reactant_roles,
                "template_id": template_id,
                **(
                    {
                        "post_build_conversions": dict(
                            candidate.metadata["post_build_conversions"]
                        )
                    }
                    if "post_build_conversions" in candidate.metadata
                    else {}
                ),
            },
        )

    def _validation_metadata(self, report: CoarseValidationReport, *, cif_path: str | None) -> dict[str, object]:
        metrics = dict(report.metrics)
        metrics["cif_path"] = cif_path
        return {
            "classification": report.classification,
            "is_valid": report.is_valid,
            "passes_hard_validation": report.passes_hard_validation,
            "blocks_cif_export": report.blocks_cif_export,
            "reasons": list(report.reasons),
            "hard_hard_invalid_reasons": list(report.hard_hard_invalid_reasons),
            "warning_reasons": list(report.warning_reasons),
            "hard_invalid_reasons": list(report.hard_invalid_reasons),
            "metrics": self._json_safe(metrics),
        }

    def _hard_hard_validation_metadata(
        self,
        *,
        reasons: tuple[str, ...],
        metrics: Mapping[str, object],
    ) -> dict[str, object]:
        return {
            "classification": "hard_hard_invalid",
            "is_valid": False,
            "passes_hard_validation": False,
            "blocks_cif_export": True,
            "reasons": list(reasons),
            "hard_hard_invalid_reasons": list(reasons),
            "warning_reasons": [],
            "hard_invalid_reasons": [],
            "metrics": self._json_safe({"cif_path": None, **dict(metrics)}),
        }

    def _classified_cif_destination(
        self,
        out_dir: str | Path,
        *,
        structure_id: str,
        classification: str,
    ) -> Path:
        output_root = Path(out_dir)
        if not self.config.separate_cif_outputs_by_validation:
            return output_root / f"{structure_id}.cif"
        bucket = "valid"
        if classification == "warning":
            bucket = "warning"
        elif classification in {"hard_invalid", "hard_hard_invalid"}:
            bucket = "invalid"
        return output_root / bucket / f"{structure_id}.cif"

    def _write_classified_candidate_cif(
        self,
        *,
        out_dir: str | Path,
        structure_id: str,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec],
        provisional_summary: BatchPairSummary,
    ) -> tuple[str | None, dict[str, object] | None]:
        if not self.config.separate_cif_outputs_by_validation:
            out_path = Path(out_dir) / f"{structure_id}.cif"
            self.cif_writer.write_candidate(
                out_path,
                candidate,
                monomer_specs,
                data_name=structure_id,
            )
            return str(out_path), None

        staging_dir = Path(out_dir) / ".staging"
        staging_dir.mkdir(parents=True, exist_ok=True)
        staging_path = staging_dir / f"{structure_id}.cif"
        self.cif_writer.write_candidate(
            staging_path,
            candidate,
            monomer_specs,
            data_name=structure_id,
        )
        validation_record = self._summary_to_dict(replace(provisional_summary, cif_path=str(staging_path)))
        report = self.structure_validator.validate_manifest_record(validation_record)
        if report.classification == "hard_hard_invalid":
            staging_path.unlink(missing_ok=True)
            return None, self._validation_metadata(report, cif_path=None)
        final_path = self._classified_cif_destination(
            out_dir,
            structure_id=structure_id,
            classification=report.classification,
        )
        final_path.parent.mkdir(parents=True, exist_ok=True)
        staging_path.replace(final_path)
        return str(final_path), self._validation_metadata(report, cif_path=str(final_path))

    def _generate_pair_candidates_from_monomers(
        self,
        *,
        pair: _BinaryBridgePairContext,
        out_dir: str | Path | None,
        write_cif: bool | None,
        cif_export_start_index: int,
    ) -> tuple[tuple[BatchPairSummary, ...], tuple[Candidate, ...], int]:
        first = pair.first
        second = pair.second
        first_id = pair.first_id
        second_id = pair.second_id
        pair_id = pair.pair_id
        first_role, second_role = pair.role_ids
        role_record_ids = {
            first_role: first_id,
            second_role: second_id,
        }
        role_connectivities = {
            first_role: len(first.motifs),
            second_role: len(second.motifs),
        }

        try:
            evaluation = self._evaluate_pair_topologies(pair)
        except ValueError as exc:
            pair_mode = self._pair_mode_from_connectivities(len(first.motifs), len(second.motifs))
            return (
                (
                    BatchPairSummary(
                        structure_id=pair_id,
                        pair_id=pair_id,
                        pair_mode=pair_mode,
                        status="unsupported-connectivity",
                        reactant_a_record_id=first_id,
                        reactant_b_record_id=second_id,
                        reactant_a_connectivity=len(first.motifs),
                        reactant_b_connectivity=len(second.motifs),
                        metadata={
                            "error": f"{type(exc).__name__}: {exc}",
                            "reactant_record_ids": role_record_ids,
                            "reactant_connectivities": role_connectivities,
                            "reactant_roles": pair.role_ids,
                            "template_id": pair.template.id,
                        },
                    ),
                ),
                (),
                0,
            )

        attempted_structures = len(evaluation.available_topologies)
        if not evaluation.candidates:
            error_message = (
                str(evaluation.failed_topologies.get("__pair__"))
                if "__pair__" in evaluation.failed_topologies
                else "all topology-specific generation attempts failed"
            )
            return (
                (
                    BatchPairSummary(
                        structure_id=pair_id,
                        pair_id=pair_id,
                        pair_mode=evaluation.pair_mode,
                        status="generation-failed",
                        reactant_a_record_id=first_id,
                        reactant_b_record_id=second_id,
                        reactant_a_connectivity=len(first.motifs),
                        reactant_b_connectivity=len(second.motifs),
                        metadata={
                            "available_topologies": evaluation.available_topologies,
                            "failed_topologies": dict(evaluation.failed_topologies),
                            "error": error_message,
                            "reactant_record_ids": role_record_ids,
                            "reactant_connectivities": role_connectivities,
                            "reactant_roles": pair.role_ids,
                            "template_id": pair.template.id,
                        },
                    ),
                ),
                (),
                attempted_structures,
            )

        effective_write_cif = self.config.write_cif if write_cif is None else write_cif
        summaries: list[BatchPairSummary] = []
        local_cifs_written = 0
        ordered_candidates = tuple(sorted(evaluation.candidates, key=lambda candidate: candidate.score, reverse=True))
        for topology_rank, candidate in enumerate(ordered_candidates, start=1):
            topology_id = candidate.metadata["net_plan"]["topology"]
            structure_id = f"{pair_id}__{topology_id}" if topology_id is not None else pair_id
            hard_hard_invalid_reasons, hard_hard_invalid_metrics = self._hard_hard_invalid_reasons(candidate)
            export_index = 0
            provisional_summary = self._candidate_to_summary(
                pair_id=pair_id,
                structure_id=structure_id,
                pair_mode=evaluation.pair_mode,
                reactant_a_id=first_id,
                reactant_b_id=second_id,
                reactant_a_connectivity=len(first.motifs),
                reactant_b_connectivity=len(second.motifs),
                candidate=candidate,
                cif_path=None,
                cif_export_index=export_index,
                topology_rank=topology_rank,
                topology_count=len(ordered_candidates),
                reactant_record_ids=role_record_ids,
                reactant_connectivities=role_connectivities,
                reactant_roles=pair.role_ids,
                template_id=pair.template.id,
                hard_hard_invalid_reasons=hard_hard_invalid_reasons,
                hard_hard_invalid_metrics=hard_hard_invalid_metrics,
            )
            cif_path = None
            validation_metadata = None
            if (
                effective_write_cif
                and out_dir is not None
                and not hard_hard_invalid_reasons
                and self._should_write_cif(write_cif, cif_export_start_index + local_cifs_written)
            ):
                export_index = cif_export_start_index + local_cifs_written + 1
                provisional_summary = replace(
                    provisional_summary,
                    metadata={
                        **dict(provisional_summary.metadata),
                        "cif_export_index": export_index,
                    },
                )
                cif_path, validation_metadata = self._write_classified_candidate_cif(
                    out_dir=out_dir,
                    structure_id=structure_id,
                    candidate=candidate,
                    monomer_specs=(first, second),
                    provisional_summary=provisional_summary,
                )
                if cif_path is not None:
                    local_cifs_written += 1
            elif hard_hard_invalid_reasons:
                validation_metadata = self._hard_hard_validation_metadata(
                    reasons=hard_hard_invalid_reasons,
                    metrics=hard_hard_invalid_metrics,
                )

            summary = provisional_summary
            if cif_path is not None:
                summary = replace(
                    summary,
                    cif_path=cif_path,
                    metadata={
                        **dict(summary.metadata),
                        **({"validation": validation_metadata} if validation_metadata is not None else {}),
                    },
                )
            elif validation_metadata is not None:
                summary = replace(
                    summary,
                    metadata={
                        **dict(summary.metadata),
                        "validation": validation_metadata,
                    },
                )

            summaries.append(summary)
        return tuple(summaries), ordered_candidates, attempted_structures

    def _hard_hard_invalid_reasons(self, candidate: Candidate) -> tuple[tuple[str, ...], dict[str, object]]:
        bridge_metrics = tuple(self._bridge_event_metrics(candidate))
        actual_distances = tuple(
            float(item["actual_distance"])
            for item in bridge_metrics
            if item.get("actual_distance") is not None
        )
        max_actual_distance = max(actual_distances) if actual_distances else None
        reasons: list[str] = []
        if (
            max_actual_distance is not None
            and max_actual_distance >= self.config.hard_hard_max_bridge_distance
        ):
            reasons.append("bridge_distance_exceeds_cif_export_limit")
        return tuple(reasons), {
            "max_actual_bridge_distance": max_actual_distance,
            "hard_hard_max_bridge_distance": self.config.hard_hard_max_bridge_distance,
        }

    def _bridge_event_metrics(self, candidate: Candidate) -> Iterable[Mapping[str, object]]:
        score_metadata = candidate.metadata.get("score_metadata", {})
        if not isinstance(score_metadata, Mapping):
            return ()
        metrics = score_metadata.get("bridge_event_metrics", ())
        return (
            item
            for item in metrics
            if isinstance(item, Mapping)
        )

    def _annotate_topology_selection(
        self,
        candidate: Candidate,
        *,
        pair_mode: str,
        available_topologies: tuple[str, ...],
        failed_topologies: Mapping[str, str],
        selected_topology: str | None,
        evaluation_mode: str,
    ) -> Candidate:
        metadata = dict(candidate.metadata)
        metadata["topology_selection"] = {
            "pair_mode": pair_mode,
            "evaluation_mode": evaluation_mode,
            "available_topologies": tuple(available_topologies),
            "evaluated_topologies": tuple(
                topology_id for topology_id in available_topologies if topology_id not in failed_topologies
            ),
            "failed_topologies": dict(failed_topologies),
            "selected_topology": selected_topology,
        }
        return replace(candidate, metadata=metadata)

    def _single_node_cell_kind(self, cell: tuple[Vec3, Vec3, Vec3]) -> str:
        first, second, _ = cell
        first_norm = norm(first)
        second_norm = norm(second)
        if first_norm < 1e-8 or second_norm < 1e-8:
            return "oblique"
        cosine = dot(first, second) / (first_norm * second_norm)
        if abs(first_norm - second_norm) < 1e-3 and abs(cosine) < 1e-3:
            return "square"
        if abs(first_norm - second_norm) < 1e-3 and abs(cosine - 0.5) < 1e-3:
            return "hexagonal"
        return "oblique"

    def _topology_cell_kind(self, cell: tuple[Vec3, Vec3, Vec3], dimensionality: str) -> str:
        if dimensionality == "3D":
            return self._spatial_cell_kind(cell)
        return self._single_node_cell_kind(cell)

    def _spatial_cell_kind(self, cell: tuple[Vec3, Vec3, Vec3]) -> str:
        lengths = tuple(norm(vector) for vector in cell)
        if any(length < 1e-8 for length in lengths):
            return "triclinic"
        cosines = (
            dot(cell[0], cell[1]) / (lengths[0] * lengths[1]),
            dot(cell[0], cell[2]) / (lengths[0] * lengths[2]),
            dot(cell[1], cell[2]) / (lengths[1] * lengths[2]),
        )
        if all(abs(cosine) < 1e-3 for cosine in cosines):
            if max(lengths) - min(lengths) < 1e-3:
                return "cubic"
            return "orthorhombic"
        return "triclinic"

    def _ordered_planar_motif_ids(self, spec: MonomerSpec) -> tuple[str, ...]:
        motifs = sorted(
            spec.motifs,
            key=lambda motif: (atan2(motif.frame.origin[1], motif.frame.origin[0]), motif.id),
        )
        return tuple(motif.id for motif in motifs)

    def _rotation_for_planar_motifs(
        self,
        spec: MonomerSpec,
        target_directions: tuple[Vec3, ...],
        *,
        template_id: str | None = None,
    ) -> tuple[Mat3, tuple[str, ...], tuple[float, ...]]:
        motifs = sorted(
            spec.motifs,
            key=lambda motif: (atan2(motif.frame.origin[1], motif.frame.origin[0]), motif.id),
        )
        motif_vectors = [
            (
                effective_motif_origin(template_id, spec, motif)
                if norm(effective_motif_origin(template_id, spec, motif)) > 1e-6
                else motif.frame.primary
            )
            for motif in motifs
        ]
        if len(motif_vectors) != len(target_directions):
            rotation = (
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            )
            return rotation, tuple(motif.id for motif in motifs), tuple(0.0 for _ in target_directions)

        local_angles = [atan2(vector[1], vector[0]) for vector in motif_vectors]
        target_angles = [atan2(direction[1], direction[0]) for direction in target_directions]
        diffs = [
            self._wrap_angle(target_angle - local_angle)
            for local_angle, target_angle in zip(local_angles, target_angles)
        ]
        best_angle = sum(diffs) / len(diffs)
        rotation = (
            (cos(best_angle), -sin(best_angle), 0.0),
            (sin(best_angle), cos(best_angle), 0.0),
            (0.0, 0.0, 1.0),
        )
        projected = tuple(
            dot(matmul_vec(rotation, vector), direction)
            for vector, direction in zip(motif_vectors, target_directions)
        )
        return rotation, tuple(motif.id for motif in motifs), projected

    def _rotation_for_spatial_motifs(
        self,
        spec: MonomerSpec,
        target_directions: tuple[Vec3, ...],
        *,
        template_id: str | None = None,
    ) -> tuple[Mat3, tuple[str, ...], tuple[float, ...]]:
        cache_key = (
            spec.id,
            template_id,
            tuple(tuple(round(value, 6) for value in direction) for direction in target_directions),
        )
        with self._spatial_rotation_cache_lock:
            cached = self._spatial_rotation_cache.get(cache_key)
        if cached is not None:
            return cached

        motifs = tuple(spec.motifs)
        if len(motifs) != len(target_directions):
            result = (
                (
                    (1.0, 0.0, 0.0),
                    (0.0, 1.0, 0.0),
                    (0.0, 0.0, 1.0),
                ),
                tuple(motif.id for motif in motifs),
                tuple(0.0 for _ in target_directions),
            )
            with self._spatial_rotation_cache_lock:
                existing = self._spatial_rotation_cache.get(cache_key)
                if existing is not None:
                    return existing
                self._spatial_rotation_cache[cache_key] = result
            return result

        local_directions = tuple(self._motif_direction_seed(motif) for motif in motifs)
        normalized_targets = tuple(self._safe_normalize(direction) for direction in target_directions)
        best_score = None
        best_result = None
        for motif_indices in permutations(range(len(motifs))):
            rotation = self._spatial_rotation_from_assignment(local_directions, normalized_targets, motif_indices)
            if rotation is None:
                continue
            alignment_score = 0.0
            offsets: list[float] = []
            for target_index, motif_index in enumerate(motif_indices):
                rotated_direction = self._safe_normalize(matmul_vec(rotation, local_directions[motif_index]))
                alignment_score += dot(rotated_direction, normalized_targets[target_index])
                offsets.append(
                    dot(
                        matmul_vec(
                            rotation,
                            effective_motif_origin(template_id, spec, motifs[motif_index]),
                        ),
                        normalized_targets[target_index],
                    )
                )
            if best_score is None or alignment_score > best_score:
                best_score = alignment_score
                best_result = (
                    rotation,
                    tuple(motifs[motif_index].id for motif_index in motif_indices),
                    tuple(offsets),
                )

        if best_result is None:
            best_result = (
                (
                    (1.0, 0.0, 0.0),
                    (0.0, 1.0, 0.0),
                    (0.0, 0.0, 1.0),
                ),
                tuple(motif.id for motif in motifs),
                tuple(0.0 for _ in target_directions),
            )
        with self._spatial_rotation_cache_lock:
            existing = self._spatial_rotation_cache.get(cache_key)
            if existing is not None:
                return existing
            self._spatial_rotation_cache[cache_key] = best_result
        return best_result

    def _spatial_rotation_from_assignment(
        self,
        local_directions: tuple[Vec3, ...],
        target_directions: tuple[Vec3, ...],
        motif_indices: tuple[int, ...],
    ) -> Mat3 | None:
        for second_target_index in range(1, len(target_directions)):
            first_local = local_directions[motif_indices[0]]
            second_local = local_directions[motif_indices[second_target_index]]
            local_normal = cross(first_local, second_local)
            if norm(local_normal) < 1e-6:
                continue

            first_target = target_directions[0]
            second_target = target_directions[second_target_index]
            target_normal = cross(first_target, second_target)
            if norm(target_normal) < 1e-6:
                continue

            return rotation_from_frame_to_axes(
                Frame(origin=(0.0, 0.0, 0.0), primary=first_local, normal=local_normal),
                first_target,
                target_normal,
            )
        return None

    def _motif_direction_seed(self, motif) -> Vec3:
        if norm(motif.frame.origin) > 1e-6:
            return self._safe_normalize(motif.frame.origin)
        return self._safe_normalize(motif.frame.primary)

    def _linker_layout_for_axes(
        self,
        spec: MonomerSpec,
        target_direction: Vec3,
        target_normal: Vec3,
        *,
        template_id: str | None = None,
    ) -> _LinkerLayout:
        if len(spec.motifs) != 2:
            raise ValueError("single-node node-linker placement requires a ditopic linker")
        first, second = spec.motifs
        first_origin = effective_motif_origin(template_id, spec, first)
        second_origin = effective_motif_origin(template_id, spec, second)
        local_primary = sub(second_origin, first_origin)
        if norm(local_primary) < 1e-8:
            local_primary = first.frame.primary
        normal_seed = add(first.frame.normal, second.frame.normal)
        if norm(normal_seed) < 1e-8:
            normal_seed = target_normal
        synthetic_frame = Frame(origin=(0.0, 0.0, 0.0), primary=local_primary, normal=normal_seed)
        rotation = rotation_from_frame_to_axes(synthetic_frame, target_direction, target_normal)
        transformed_a = matmul_vec(rotation, first_origin)
        transformed_b = matmul_vec(rotation, second_origin)
        projected_a = dot(transformed_a, target_direction)
        projected_b = dot(transformed_b, target_direction)
        if projected_a <= projected_b:
            first_motif, second_motif = first, second
            first_local_origin, second_local_origin = first_origin, second_origin
            first_projected, second_projected = projected_a, projected_b
        else:
            first_motif, second_motif = second, first
            first_local_origin, second_local_origin = second_origin, first_origin
            first_projected, second_projected = projected_b, projected_a
        midpoint = 0.5 * (first_projected + second_projected)
        return _LinkerLayout(
            rotation=rotation,
            motif_ids=(first_motif.id, second_motif.id),
            local_origins=(first_local_origin, second_local_origin),
            offset_a=midpoint - first_projected,
            offset_b=second_projected - midpoint,
            span=norm(sub(second_origin, first_origin)),
        )

    def _linker_layout_for_direction(
        self,
        spec: MonomerSpec,
        target_direction: Vec3,
        *,
        template_id: str | None = None,
    ) -> _LinkerLayout:
        return self._linker_layout_for_axes(spec, target_direction, (0.0, 0.0, 1.0), template_id=template_id)

    def _bridge_plane_normal(self, bridge_direction: Vec3, start_normal: Vec3, end_normal: Vec3) -> Vec3:
        combined = add(start_normal, end_normal)
        projected = sub(combined, scale(bridge_direction, dot(combined, bridge_direction)))
        if norm(projected) < 1e-8:
            fallback = (0.0, 0.0, 1.0) if abs(bridge_direction[2]) < 0.9 else (1.0, 0.0, 0.0)
            projected = sub(fallback, scale(bridge_direction, dot(fallback, bridge_direction)))
        return self._safe_normalize(projected)

    def _node_linker_edge_direction(
        self,
        *,
        motif,
        monomer: MonomerSpec,
        node_rotation_a: Mat3,
        node_rotation_b: Mat3,
        template_id: str | None = None,
    ) -> Vec3:
        local_origin = effective_motif_origin(template_id, monomer, motif)
        local_direction = local_origin if norm(local_origin) > 1e-6 else motif.frame.primary
        direction_a = self._safe_normalize(matmul_vec(node_rotation_a, local_direction))
        direction_b = scale(self._safe_normalize(matmul_vec(node_rotation_b, local_direction)), -1.0)
        combined = add(direction_a, direction_b)
        if norm(combined) < 1e-8:
            combined = direction_a
        return self._safe_normalize(combined)

    def _generalized_single_node_edge_vector(
        self,
        motif,
        monomer: MonomerSpec,
        *,
        node_rotation_a: Mat3,
        node_rotation_b: Mat3,
        direction: Vec3,
        reactive_site_span: float,
        template_id: str | None = None,
    ) -> Vec3:
        local_origin = effective_motif_origin(template_id, monomer, motif)
        rotated_a = matmul_vec(node_rotation_a, local_origin)
        rotated_b = matmul_vec(node_rotation_b, local_origin)
        return add(scale(direction, reactive_site_span), sub(rotated_a, rotated_b))

    def _generalized_single_node_cell(self, edge_vectors: tuple[Vec3, Vec3, Vec3]) -> tuple[Vec3, Vec3, Vec3]:
        edge_one, edge_two, edge_three = edge_vectors
        return (
            sub(edge_one, edge_two),
            sub(edge_one, edge_three),
            (0.0, 0.0, self.config.embedding_config.default_layer_spacing),
        )

    def _world_motif_origin(
        self,
        cell: tuple[Vec3, Vec3, Vec3],
        pose: Pose,
        monomer: MonomerSpec,
        motif,
        image: tuple[int, int, int],
        *,
        template_id: str | None = None,
    ) -> Vec3:
        base = add(
            pose.translation,
            matmul_vec(
                pose.rotation_matrix,
                effective_motif_origin(template_id, monomer, motif),
            ),
        )
        return add(base, self._periodic_offset(cell, image))

    def _periodic_offset(self, cell: tuple[Vec3, Vec3, Vec3], image: tuple[int, int, int]) -> Vec3:
        return add(add(scale(cell[0], image[0]), scale(cell[1], image[1])), scale(cell[2], image[2]))

    def _safe_normalize(self, vector: Vec3) -> Vec3:
        if norm(vector) < 1e-8:
            return (1.0, 0.0, 0.0)
        return normalize(vector)

    def _solve_symmetric_3x3(self, matrix: list[list[float]], rhs: list[float]) -> tuple[float, float, float] | None:
        a = [row[:] for row in matrix]
        b = rhs[:]
        size = 3
        for pivot_index in range(size):
            pivot_row = max(range(pivot_index, size), key=lambda row_index: abs(a[row_index][pivot_index]))
            if abs(a[pivot_row][pivot_index]) < 1e-8:
                return None
            if pivot_row != pivot_index:
                a[pivot_index], a[pivot_row] = a[pivot_row], a[pivot_index]
                b[pivot_index], b[pivot_row] = b[pivot_row], b[pivot_index]
            pivot = a[pivot_index][pivot_index]
            for column in range(pivot_index, size):
                a[pivot_index][column] /= pivot
            b[pivot_index] /= pivot
            for row_index in range(size):
                if row_index == pivot_index:
                    continue
                factor = a[row_index][pivot_index]
                if abs(factor) < 1e-12:
                    continue
                for column in range(pivot_index, size):
                    a[row_index][column] -= factor * a[pivot_index][column]
                b[row_index] -= factor * b[pivot_index]
        return (b[0], b[1], b[2])

    def _edge_endpoint_key(self, edge_id: str, endpoint: str) -> str:
        return f"{edge_id}:{endpoint}"

    def _template_target_distance(self, template) -> float:
        return bridge_target_distance(template, default_bridge_distance=self.config.embedding_config.bridge_target_distance)

    def _trigonal_directions(self, start_angle: float) -> tuple[Vec3, Vec3, Vec3]:
        return tuple(
            (cos(start_angle + offset), sin(start_angle + offset), 0.0)
            for offset in (0.0, 2.0 * pi / 3.0, -2.0 * pi / 3.0)
        )  # type: ignore[return-value]

    def _wrap_angle(self, value: float) -> float:
        while value <= -pi:
            value += 2.0 * pi
        while value > pi:
            value -= 2.0 * pi
        return value

    def _should_write_cif(self, write_cif: bool | None, cifs_written: int) -> bool:
        effective = self.config.write_cif if write_cif is None else write_cif
        if not effective:
            return False
        if self.config.max_cif_exports is None:
            return True
        return cifs_written < self.config.max_cif_exports

    def _summary_to_dict(self, summary: BatchPairSummary) -> dict[str, object]:
        return {
            "structure_id": summary.structure_id,
            "pair_id": summary.pair_id,
            "pair_mode": summary.pair_mode,
            "status": summary.status,
            "reactant_a_record_id": summary.reactant_a_record_id,
            "reactant_b_record_id": summary.reactant_b_record_id,
            "reactant_a_connectivity": summary.reactant_a_connectivity,
            "reactant_b_connectivity": summary.reactant_b_connectivity,
            "amine_record_id": summary.amine_record_id,
            "aldehyde_record_id": summary.aldehyde_record_id,
            "amine_connectivity": summary.amine_connectivity,
            "aldehyde_connectivity": summary.aldehyde_connectivity,
            "reactant_record_ids": dict(summary.reactant_record_ids),
            "reactant_connectivities": dict(summary.reactant_connectivities),
            "reactant_roles": tuple(summary.reactant_roles),
            "topology_id": summary.topology_id,
            "score": summary.score,
            "flags": summary.flags,
            "cif_path": summary.cif_path,
            "metadata": summary.metadata,
        }

    def _write_summary_report(self, path: Path, summary: BatchRunSummary) -> None:
        lines = [
            "# Batch binary-bridge generation summary",
            "",
            f"- Input directory: `{summary.input_dir}`",
            f"- Output directory: `{summary.output_dir}`",
            f"- Manifest: `{summary.manifest_path}`",
            f"- Built monomers: {summary.built_monomers}",
            f"- Failed monomer builds: {summary.failed_monomers}",
            f"- Attempted pairs: {summary.attempted_pairs}",
            f"- Successful pairs: {summary.successful_pairs}",
            f"- Attempted structures: {summary.attempted_structures}",
            f"- Successful structures: {summary.successful_structures}",
            f"- CIFs written: {summary.cifs_written}",
            f"- Mode counts (structures): {dict(summary.mode_counts)}",
            f"- Topology counts: {dict(summary.topology_counts)}",
            "",
            "## Top results",
            "",
        ]
        if not summary.top_results:
            lines.append("- No successful pairs were generated.")
        else:
            for result in summary.top_results[:10]:
                lines.append(
                    f"- `{result.structure_id}` ({result.pair_mode}) score `{result.score:.6f}` topology `{result.topology_id}`"
                )
        if summary.build_failures:
            lines.extend(
                [
                    "",
                    "## Build failures",
                    "",
                ]
            )
            for record_id, error in sorted(summary.build_failures.items())[:20]:
                lines.append(f"- `{record_id}`: {error}")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def _json_safe(self, value):
        if isinstance(value, (str, int, float, bool)) or value is None:
            return value
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, Mapping):
            return {str(key): self._json_safe(item) for key, item in value.items()}
        if isinstance(value, tuple):
            return [self._json_safe(item) for item in value]
        if isinstance(value, list):
            return [self._json_safe(item) for item in value]
        return str(value)


__all__ = [
    "BatchGenerationConfig",
    "BatchMonomerRecord",
    "BatchPairSummary",
    "BatchRunSummary",
    "BatchStructureGenerator",
    "BuiltBatchMonomer",
]
