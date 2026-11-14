"""Prepared topology caches preserve graph semantics and source isolation."""

from dataclasses import replace
import math
import random
from unittest.mock import patch

import pytest

from cofkit import decompose, topology_data, topology_symmetry as symmetry
from cofkit.topology_index import TopologyDefinition, TopologyNodeDefinition, TopologyEdgeDefinition
from cofkit.topology_repository import TopologyRepository
from topology_reference import exhaustive_topology_expansion


def p1_definition(gain=1):
    return TopologyDefinition(
        id="cache-test", name="cache-test", dimensionality="3D",
        node_definitions=(TopologyNodeDefinition("node", (0, 0, 0), 6),),
        edge_definitions=tuple(TopologyEdgeDefinition((0, 0, 0), endpoint) for endpoint in (
            (1, 0, 0), (0, 1, 0), (0, 0, gain),
        )), metadata={"space_group": "P 1"},
    )


def test_workspace_repository_reuses_parsed_index():
    with patch.object(topology_data, "_read_json", wraps=topology_data._read_json) as read:
        repository = TopologyRepository.from_workspace_data()
        for identifier in ("hcb", "sql", "pcu", "dia"):
            assert repository.load(identifier).id == identifier
        assert read.call_count == 1


def test_same_identifier_in_different_repositories_cannot_share_wrong_graph():
    first, second = p1_definition(), p1_definition(2)
    a = TopologyRepository(definitions={first.id: first})
    b = TopologyRepository(definitions={second.id: second})
    graph_a = decompose._topology_definition_graph(a.load(first.id))
    graph_b = decompose._topology_definition_graph(b.load(second.id))
    assert graph_a != graph_b
    assert graph_a.gain_edges[-1][2] == (0, 0, -1)
    assert graph_b.gain_edges[-1][2] == (0, 0, -2)


def test_cache_keys_include_structure_and_mutable_space_group_metadata():
    definition = p1_definition()
    symmetry._cached_topology_expansion.cache_clear()
    with patch.object(symmetry, "_gemmi", wraps=symmetry._gemmi) as gemmi_loader:
        first = symmetry.expand_topology_definition(definition)
        definition.metadata["annotation"] = "irrelevant to graph"
        assert symmetry.expand_topology_definition(definition) is first
        assert gemmi_loader.call_count == 1
        changed = replace(definition, node_definitions=(replace(definition.node_definitions[0], connectivity=4),))
        assert symmetry.expand_topology_definition(changed).node_sites[0].connectivity == 4
        assert gemmi_loader.call_count == 2
        definition.metadata["space_group"] = "nonexistent-space-group"
        with pytest.raises(ValueError, match="not recognized"):
            symmetry.expand_topology_definition(definition)


def test_derived_graph_reuses_only_immutable_expansion_contents():
    decompose._expanded_topology_graph.cache_clear()
    first = decompose._topology_definition_graph(p1_definition())
    equivalent = decompose._topology_definition_graph(p1_definition())
    assert equivalent is first
    assert decompose._topology_definition_graph(p1_definition(2)) != first


def test_malformed_definition_keeps_linear_oracle_error():
    definition = replace(p1_definition(), node_definitions=(TopologyNodeDefinition("empty", (), 6),))
    errors = []
    for method in (symmetry.expand_topology_definition, exhaustive_topology_expansion):
        with pytest.raises(ValueError) as error:
            method(definition)
        errors.append(str(error.value))
    assert errors[0] == errors[1]


@pytest.mark.parametrize("identifier", [
    "hcb", "sql", "kgm", "hca", "dia", "pcu", "bex", "car", "bew", "cor", "cpa",
])
def test_indexed_symmetry_expansion_matches_frozen_linear_oracle(identifier):
    definition = TopologyRepository.from_workspace_data().load(identifier)
    assert symmetry.expand_topology_definition(definition) == exhaustive_topology_expansion(definition)


@pytest.mark.parametrize("dimension", [2, 3])
def test_coordinate_index_preserves_boundary_ties_and_first_match(dimension):
    rng = random.Random(76123)
    tolerance = symmetry._FRACTIONAL_WRAP_TOLERANCE
    positions = [tuple(rng.random() for _ in range(dimension)) for _ in range(100)]
    for boundary in (0, tolerance, 0.25, .5, 1 - tolerance, 1):
        for first in (math.nextafter(boundary, -math.inf), boundary, math.nextafter(boundary, math.inf)):
            positions.append((first,) + (0.25,) * (dimension - 1))
    index = symmetry._PositionIndex()
    for position in positions:
        expected = next((i for i, existing in enumerate(index.positions) if symmetry._same_fractional_position(position, existing)), None)
        assert index.match(position) == expected
        index.add(position)
    for position in positions:
        for offset in (-tolerance, 0, tolerance):
            query = (position[0] + offset,) + position[1:]
            expected = next((i for i, existing in enumerate(index.positions) if symmetry._same_fractional_position(query, existing)), None)
            assert index.match(query) == expected


def test_coordinate_index_prunes_linear_exact_comparisons():
    index = symmetry._PositionIndex()
    for i in range(1000):
        index.add((i / 1000, .25, .25))
    with patch.object(symmetry, "_same_fractional_position", wraps=symmetry._same_fractional_position) as same:
        assert index.match((.999, .25, .25)) == 999
    assert same.call_count < 5
