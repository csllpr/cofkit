"""Wiring assertions for the single-owner cell-default constants.

Guards W1.1/W1.7 of agent-docs/MAGIC_NUMBER_FIX_PLAN.md: every cell-default
owner and the ring-forming CLI argparse default must reference the constants
in ``cofkit.constants``, never retype the literal.
"""

from cofkit import constants, topology_analysis, topology_symmetry
from cofkit.build_workflows.ring_forming import RingFormationConfig
from cofkit.cli import build_parser
from cofkit.embedding import EmbeddingConfig
from cofkit.engine import COFEngineConfig
from cofkit.model import AssemblyState


def test_embedding_config_references_cell_constants():
    config = EmbeddingConfig()
    assert config.default_layer_spacing == constants.DEFAULT_MONOLAYER_C_ANGSTROM
    assert config.default_lateral_span == constants.DEFAULT_LATERAL_SPAN_ANGSTROM


def test_engine_config_references_cell_constants():
    config = COFEngineConfig()
    assert config.default_layer_spacing == constants.DEFAULT_MONOLAYER_C_ANGSTROM
    assert config.default_ring_layer_spacing == constants.DEFAULT_MONOLAYER_C_ANGSTROM
    assert config.default_lateral_span == constants.DEFAULT_LATERAL_SPAN_ANGSTROM


def test_assembly_state_default_cell_references_cell_constants():
    cell = AssemblyState().cell
    assert cell == (
        (constants.DEFAULT_LATERAL_SPAN_ANGSTROM, 0.0, 0.0),
        (0.0, constants.DEFAULT_LATERAL_SPAN_ANGSTROM, 0.0),
        (0.0, 0.0, constants.DEFAULT_MONOLAYER_C_ANGSTROM),
    )


def test_ring_formation_config_references_cell_constants():
    config = RingFormationConfig()
    assert config.layer_spacing == constants.DEFAULT_MONOLAYER_C_ANGSTROM


def test_ring_forming_cli_layer_spacing_default_references_constant():
    args = build_parser().parse_args(["build", "ring-forming"])
    assert args.layer_spacing == constants.DEFAULT_MONOLAYER_C_ANGSTROM


def test_fractional_wrap_tolerance_has_single_owner():
    assert topology_analysis._FRACTIONAL_WRAP_TOLERANCE is topology_symmetry._FRACTIONAL_WRAP_TOLERANCE
