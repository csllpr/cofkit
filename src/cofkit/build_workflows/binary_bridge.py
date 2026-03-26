from __future__ import annotations

from ..batch import BatchGenerationConfig, BatchStructureGenerator

WORKFLOW_ID = "binary_bridge"

# This adapter module is intentionally thin. The existing production builder stays in
# cofkit.batch, while future build families grow in sibling modules instead of editing
# binary-bridge control flow in place.
DEFAULT_GENERATOR_CLASS = BatchStructureGenerator
DEFAULT_CONFIG_CLASS = BatchGenerationConfig

__all__ = [
    "DEFAULT_CONFIG_CLASS",
    "DEFAULT_GENERATOR_CLASS",
    "WORKFLOW_ID",
]
