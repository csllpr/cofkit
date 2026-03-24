from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from .engine import COFProject
from .indexed_topology_layouts import ExpandedIndexedTopology
from .model import Candidate, MonomerSpec, ReactionTemplate
from .planner import TopologyHint
from .single_node_topologies import ExpandedSingleNodeTopology


@dataclass(frozen=True)
class PairTopologyBuildRequest:
    project: COFProject
    monomers: tuple[MonomerSpec, MonomerSpec]
    template: ReactionTemplate
    topology: TopologyHint
    expanded: ExpandedSingleNodeTopology | ExpandedIndexedTopology
    builder_family: str
    pair_mode: str
    connectivities: tuple[int, int]


@dataclass(frozen=True)
class PairTopologyBuilder:
    id: str
    description: str
    supports: Callable[[PairTopologyBuildRequest], bool]
    build: Callable[[PairTopologyBuildRequest], Candidate]


class PairTopologyBuilderRegistry:
    def __init__(self, builders: tuple[PairTopologyBuilder, ...]):
        self._builders = builders

    def matching_builders(self, request: PairTopologyBuildRequest) -> tuple[PairTopologyBuilder, ...]:
        return tuple(builder for builder in self._builders if builder.supports(request))

    def can_build(self, request: PairTopologyBuildRequest) -> bool:
        return any(builder.supports(request) for builder in self._builders)

    def build(self, request: PairTopologyBuildRequest) -> Candidate:
        for builder in self._builders:
            if builder.supports(request):
                return builder.build(request)
        raise ValueError(
            "no registered topology builder can handle "
            f"topology={request.topology.id!r}, pair_mode={request.pair_mode!r}, "
            f"connectivities={request.connectivities!r}, dimensionality={request.topology.dimensionality!r}, "
            f"builder_family={request.builder_family!r}"
        )
