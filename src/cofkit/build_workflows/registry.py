from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from ..reactions import ReactionLibrary


@dataclass(frozen=True)
class BuildWorkflowDefinition:
    id: str
    description: str
    implementation_status: str
    builder_module: str
    isolation_boundary: str
    template_workflow_family: str | None = None
    explicit_selection_required: bool = False
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class BuildWorkflowRegistry:
    workflows: Mapping[str, BuildWorkflowDefinition]

    def get(self, workflow_id: str) -> BuildWorkflowDefinition:
        try:
            return self.workflows[workflow_id]
        except KeyError as exc:
            raise KeyError(f"unknown build workflow {workflow_id!r}") from exc

    def list(self) -> tuple[BuildWorkflowDefinition, ...]:
        return tuple(self.workflows.values())

    def resolve_template_workflow(
        self,
        reaction_library: ReactionLibrary,
        template_ids: tuple[str, ...],
    ) -> BuildWorkflowDefinition:
        if not template_ids:
            raise ValueError("at least one template id is required to resolve a build workflow")

        families = tuple(
            dict.fromkeys(reaction_library.workflow_family(template_id) for template_id in template_ids)
        )
        if len(families) != 1:
            raise ValueError(
                "build workflows must be isolated by workflow family; "
                f"got mixed families {families!r} for templates {template_ids!r}"
            )

        family = families[0]
        for workflow in self.workflows.values():
            if workflow.template_workflow_family == family:
                return workflow
        raise KeyError(f"no build workflow is registered for template workflow family {family!r}")

    @classmethod
    def builtin(cls) -> "BuildWorkflowRegistry":
        workflows = (
            BuildWorkflowDefinition(
                id="binary_bridge",
                description=(
                    "Current production workflow for two-role topology-guided bridge-forming COF assembly."
                ),
                implementation_status="available",
                builder_module="cofkit.build_workflows.binary_bridge",
                isolation_boundary=(
                    "Owns the existing BatchStructureGenerator path and remains isolated from ring-forming "
                    "and composite-topology builders."
                ),
                template_workflow_family="binary_bridge",
                metadata={
                    "entrypoints": ("BatchStructureGenerator", "cofkit build single-pair", "cofkit build batch-binary-bridge"),
                    "topology_assignment_mode": "topology_guided",
                },
            ),
            BuildWorkflowDefinition(
                id="ring_forming",
                description=(
                    "Reserved workflow family for ring-forming linkages that bypass the current binary-bridge "
                    "pair-generation assumptions."
                ),
                implementation_status="planned",
                builder_module="cofkit.build_workflows.ring_forming",
                isolation_boundary=(
                    "Requires separate event generation, realization, geometry priors, and validation so current "
                    "binary-bridge behavior stays unchanged."
                ),
                template_workflow_family="ring_forming",
                metadata={
                    "topology_assignment_mode": "topology_bypass",
                },
            ),
            BuildWorkflowDefinition(
                id="composite_topology_unit",
                description=(
                    "Reserved workflow family for monomers that map to multiple topology roles or encode "
                    "composite node-edge units inside one precursor."
                ),
                implementation_status="planned",
                builder_module="cofkit.build_workflows.composite_topology_unit",
                isolation_boundary=(
                    "Requires a dedicated topology-assignment/build path instead of extending the current "
                    "one-monomer-one-role binary-bridge model."
                ),
                explicit_selection_required=True,
                metadata={
                    "selection_mode": "explicit_builder",
                },
            ),
        )
        return cls(workflows={workflow.id: workflow for workflow in workflows})


def builtin_build_workflow_registry() -> BuildWorkflowRegistry:
    return BuildWorkflowRegistry.builtin()
