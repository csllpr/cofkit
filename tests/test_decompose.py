import contextlib
import io
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit import BatchGenerationConfig, BatchMonomerRecord, BatchStructureGenerator, CoarseValidationThresholds
from cofkit.cli import main as cli_main
from cofkit.decompose import decompose_cif_to_cofid
from cofkit.reactions import ReactionLibrary


TAPB = "C1=CC(=CC=C1C2=CC(=CC(=C2)C3=CC=C(C=C3)N)C4=CC=C(C=C4)N)N"
TFB = "C1=C(C=C(C=C1C=O)C=O)C=O"
PDA = "O=Cc1ccc(C=O)cc1"
COF42_HYDRAZIDE = "CCOc1cc(C(=O)NN)cc(C(=O)NN)c1OCC"
HYDRAZINE = "NN"
TP = "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O"
PPD = "Nc1ccc(N)cc1"
BDBA = "OB(O)c1ccc(B(O)O)cc1"
HHTP = "OC1=C(O)C=C2C(=C1)C1=CC(O)=C(O)C=C1C1=CC(O)=C(O)C=C21"
TMT = "Cc1nc(C)nc(C)n1"
BEX_D2H_ALDEHYDE = (
    "C1C=C(N(C2C=CC(C3C4C(=NON=4)C(C4C=CC(N(C5C=CC(C=O)=CC=5)C5C=CC(C=O)=CC=5)=CC=4)=CC=3)=CC=2)"
    "C2C=CC(C=O)=CC=2)C=CC=1C=O"
)
BEX_D2H_AMINE = (
    "Nc1ccc(-c2ccc3c(c2)C(=C2c4cc(-c5ccc(N)cc5)ccc4-c4ccc(-c5ccc(N)cc5)cc42)c2cc(-c4ccc(N)cc4)ccc2-3)cc1"
)
TEREPHTHALALDEHYDE = "O=Cc1ccc(C=O)cc1"
TETRA_AMINE = "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1"
TETRA_ALDEHYDE = "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1"
HEXA_AMINE = "Nc1cc2c3cc(N)c(N)cc3c3cc(N)c(N)cc3c2cc1N"


def _generator(template_id: str = "imine_bridge") -> BatchStructureGenerator:
    return BatchStructureGenerator(
        BatchGenerationConfig(
            allowed_reactions=(template_id,),
            single_node_topology_ids=("hcb",),
            topology_ids=("hcb",),
            enumerate_all_topologies=False,
            write_cif=True,
            rdkit_num_conformers=2,
        )
    )


def _record(
    record_id: str,
    smiles: str,
    motif_kind: str,
    connectivity: int,
) -> BatchMonomerRecord:
    return BatchMonomerRecord(
        id=record_id,
        name=record_id,
        smiles=smiles,
        motif_kind=motif_kind,
        expected_connectivity=connectivity,
    )


ALL_BINARY_LINKAGE_ROUND_TRIP_CASES = (
    (
        "hydrazone_bridge",
        "hydrazone",
        _record("hdz", COF42_HYDRAZIDE, "hydrazide", 2),
        _record("tfb", TFB, "aldehyde", 3),
        6,
    ),
    (
        "azine_bridge",
        "azine",
        _record("hyd", HYDRAZINE, "hydrazine", 2),
        _record("tfb", TFB, "aldehyde", 3),
        6,
    ),
    (
        "keto_enamine_bridge",
        "bken",
        _record("tp", TP, "keto_aldehyde", 3),
        _record("ppd", PPD, "amine", 2),
        6,
    ),
    (
        "boronate_ester_bridge",
        "boest",
        _record("bdba", BDBA, "boronic_acid", 2),
        _record("hhtp", HHTP, "catechol", 3),
        12,
    ),
    (
        "vinylene_bridge",
        "vinylene",
        _record("tmt", TMT, "activated_methylene", 3),
        _record("pda", PDA, "aldehyde", 2),
        6,
    ),
)
ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES = (
    (
        "imine_bridge",
        "imine",
        _record("tapb", TAPB, "amine", 3),
        _record("tfb", TFB, "aldehyde", 3),
        3,
    ),
) + ALL_BINARY_LINKAGE_ROUND_TRIP_CASES

TOPOLOGY_ENABLED_ROUND_TRIP_POOL = (
    (
        "hcb_3_3_imine",
        "imine_bridge",
        "hcb",
        _record("tapb", TAPB, "amine", 3),
        _record("tfb", TFB, "aldehyde", 3),
        False,
    ),
    (
        "hcb_3_2_imine",
        "imine_bridge",
        "hcb",
        _record("tapb", TAPB, "amine", 3),
        _record("pda", PDA, "aldehyde", 2),
        False,
    ),
    (
        "hcb_2_3_hydrazone",
        "hydrazone_bridge",
        "hcb",
        _record("hdz", COF42_HYDRAZIDE, "hydrazide", 2),
        _record("tfb", TFB, "aldehyde", 3),
        False,
    ),
    (
        "hcb_3_2_bken",
        "keto_enamine_bridge",
        "hcb",
        _record("tp", TP, "keto_aldehyde", 3),
        _record("ppd", PPD, "amine", 2),
        False,
    ),
    (
        "sql_4_2_imine",
        "imine_bridge",
        "sql",
        _record("tetra_amine", TETRA_AMINE, "amine", 4),
        _record("tpal", TEREPHTHALALDEHYDE, "aldehyde", 2),
        True,
    ),
    (
        "sql_4_4_imine",
        "imine_bridge",
        "sql",
        _record("tetra_amine", TETRA_AMINE, "amine", 4),
        _record("tetra_aldehyde", TETRA_ALDEHYDE, "aldehyde", 4),
        True,
    ),
    (
        "hxl_6_2_imine",
        "imine_bridge",
        "hxl",
        _record("hexa_amine", HEXA_AMINE, "amine", 6),
        _record("tpal", TEREPHTHALALDEHYDE, "aldehyde", 2),
        False,
    ),
    (
        "bex_4_4_imine",
        "imine_bridge",
        "bex",
        _record("bex_d2h_amine", BEX_D2H_AMINE, "amine", 4),
        _record("bex_d2h_aldehyde", BEX_D2H_ALDEHYDE, "aldehyde", 4),
        True,
    ),
)


def _cif_lines_without_cofid(source_path: str | Path) -> list[str]:
    lines = Path(source_path).read_text(encoding="utf-8").splitlines()
    if lines and lines[0].startswith("# COFid: "):
        lines = lines[1:]
    return lines


def _write_cif_lines(lines: list[str], target_path: Path) -> Path:
    target_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target_path


def _without_cofid_comment(source_path: str | Path, target_path: Path) -> Path:
    return _write_cif_lines(_cif_lines_without_cofid(source_path), target_path)


def _without_cif_bond_type_column(source_path: str | Path, target_path: Path) -> Path:
    lines = _cif_lines_without_cofid(source_path)
    output: list[str] = []
    index = 0
    while index < len(lines):
        if lines[index].strip() == "loop_":
            next_index = index + 1
            headers: list[str] = []
            while next_index < len(lines) and lines[next_index].lstrip().startswith("_"):
                headers.append(lines[next_index].strip())
                next_index += 1
            if "_geom_bond_atom_site_label_1" in headers and "_geom_bond_atom_site_label_2" in headers:
                remove_indices = {
                    header_index
                    for header_index, header in enumerate(headers)
                    if header in {"_ccdc_geom_bond_type", "_geom_bond_type"}
                }
                output.append(lines[index])
                output.extend(
                    header for header_index, header in enumerate(headers) if header_index not in remove_indices
                )
                while (
                    next_index < len(lines)
                    and lines[next_index].strip()
                    and lines[next_index].strip() != "loop_"
                    and not lines[next_index].lstrip().startswith("_")
                ):
                    row = lines[next_index].split()
                    output.append(
                        " ".join(value for value_index, value in enumerate(row) if value_index not in remove_indices)
                    )
                    next_index += 1
                index = next_index
                continue
        output.append(lines[index])
        index += 1
    return _write_cif_lines(output, target_path)


def _without_cif_bond_loop(source_path: str | Path, target_path: Path) -> Path:
    lines = _cif_lines_without_cofid(source_path)
    output: list[str] = []
    index = 0
    while index < len(lines):
        if lines[index].strip() == "loop_":
            next_index = index + 1
            headers: list[str] = []
            while next_index < len(lines) and lines[next_index].lstrip().startswith("_"):
                headers.append(lines[next_index].strip())
                next_index += 1
            if "_geom_bond_atom_site_label_1" in headers and "_geom_bond_atom_site_label_2" in headers:
                while (
                    next_index < len(lines)
                    and lines[next_index].strip()
                    and lines[next_index].strip() != "loop_"
                    and not lines[next_index].lstrip().startswith("_")
                ):
                    next_index += 1
                index = next_index
                continue
        output.append(lines[index])
        index += 1
    return _write_cif_lines(output, target_path)


def _with_generic_atom_labels(source_path: str | Path, target_path: Path) -> Path:
    lines = Path(source_path).read_text(encoding="utf-8").splitlines()
    output: list[str] = []
    index = 0
    while index < len(lines):
        if lines[index].strip() == "loop_":
            next_index = index + 1
            headers: list[str] = []
            while next_index < len(lines) and lines[next_index].lstrip().startswith("_"):
                headers.append(lines[next_index].strip())
                next_index += 1
            if "_atom_site_label" in headers:
                label_index = headers.index("_atom_site_label")
                type_index = headers.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in headers else None
                output.append(lines[index])
                output.extend(headers)
                atom_index = 1
                while (
                    next_index < len(lines)
                    and lines[next_index].strip()
                    and lines[next_index].strip() != "loop_"
                    and not lines[next_index].lstrip().startswith("_")
                ):
                    row = lines[next_index].split()
                    symbol = row[type_index] if type_index is not None else "X"
                    row[label_index] = f"{symbol}{atom_index}"
                    output.append(" ".join(row))
                    atom_index += 1
                    next_index += 1
                index = next_index
                continue
        output.append(lines[index])
        index += 1
    return _write_cif_lines(output, target_path)


class DecomposeRoundTripTests(unittest.TestCase):
    def test_decompose_generated_hcb_three_three_imine_returns_original_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 3)
        self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_decompose_auto_detects_topology_from_cofid_comment(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )

            result = decompose_cif_to_cofid(summary.cif_path)

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.topology, "hcb")
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        detection = result.metadata["topology_detection"]
        self.assertEqual(detection["confidence"], "exact")
        self.assertEqual(detection["metadata"]["source"], "cofid_comment")

    def test_decompose_generated_hcb_three_two_imine_returns_original_inputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("pda", PDA, "aldehyde", 2),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 6)
        self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_decompose_generated_decorated_bex_imine_recovers_monomers(self):
        generator = BatchStructureGenerator(
            BatchGenerationConfig(
                topology_ids=("bex",),
                stacking_ids=("AA",),
                enumerate_all_topologies=False,
                write_cif=True,
                rdkit_num_conformers=1,
                hard_hard_max_bridge_distance=10.0,
                validation_thresholds=CoarseValidationThresholds(hard_hard_max_bridge_distance=10.0),
            )
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = generator.generate_pair_candidate(
                _record("bex_d2h_amine", BEX_D2H_AMINE, "amine", 4),
                _record("bex_d2h_aldehyde", BEX_D2H_ALDEHYDE, "aldehyde", 4),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped_bex.cif")

            result = decompose_cif_to_cofid(input_cif, topology="bex")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 8)
        self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_analyze_decompose_cli_prints_recovered_cofid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_main(["analyze", "decompose", str(input_cif), "--topology", "hcb"])

        self.assertEqual(buffer.getvalue().strip(), summary.metadata["cofid"])

    def test_analyze_decompose_cli_auto_detects_topology_when_omitted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                cli_main(["analyze", "decompose", str(input_cif)])

        self.assertEqual(buffer.getvalue().strip(), summary.metadata["cofid"])

    def test_decompose_generated_hcb_cofs_for_all_build_binary_linkages(self):
        for template_id, linkage, first, second, expected_linkage_bonds in ALL_BINARY_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

                result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage=linkage)

                self.assertTrue(result.ok, result.reason)
                self.assertEqual(result.cofid, summary.metadata["cofid"])
                self.assertEqual(result.metadata[f"n_{linkage}_linkage_bonds"], expected_linkage_bonds)
                self.assertEqual(result.metadata["n_unique_monomers"], 2)

    def test_decompose_infers_bond_orders_when_cif_bond_type_column_is_absent(self):
        for template_id, linkage, first, second, expected_linkage_bonds in ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cif_bond_type_column(summary.cif_path, temp_path / "without_bond_type.cif")

                result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage=linkage)

                self.assertTrue(result.ok, result.reason)
                self.assertEqual(result.cofid, summary.metadata["cofid"])
                self.assertEqual(result.metadata["bond_source"], "explicit_cif")
                self.assertGreater(result.metadata["n_missing_cif_bond_orders"], 0)
                self.assertGreater(result.metadata["n_bond_orders_inferred"], 0)
                self.assertEqual(result.metadata[f"n_{linkage}_linkage_bonds"], expected_linkage_bonds)

    def test_decompose_infers_bonds_when_cif_bond_loop_is_absent(self):
        for template_id, linkage, first, second, expected_linkage_bonds in ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES:
            with self.subTest(template_id=template_id), tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                summary, _candidate = _generator(template_id).generate_pair_candidate(
                    first,
                    second,
                    out_dir=temp_path,
                    write_cif=True,
                )
                input_cif = _without_cif_bond_loop(summary.cif_path, temp_path / "without_bond_loop.cif")

                result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage=linkage)

                self.assertTrue(result.ok, result.reason)
                self.assertEqual(result.cofid, summary.metadata["cofid"])
                self.assertEqual(result.metadata["bond_source"], "distance_inferred")
                self.assertGreater(result.metadata["n_distance_inferred_bonds"], 0)
                self.assertGreater(result.metadata["n_bond_orders_inferred"], 0)
                self.assertEqual(result.metadata[f"n_{linkage}_linkage_bonds"], expected_linkage_bonds)

    def test_decompose_auto_detects_topology_without_cofid_bonds_or_instance_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator().generate_pair_candidate(
                _record("tapb", TAPB, "amine", 3),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            no_bonds = _without_cif_bond_loop(summary.cif_path, temp_path / "without_bond_loop.cif")
            generic_labels = _with_generic_atom_labels(no_bonds, temp_path / "generic_labels.cif")

            result = decompose_cif_to_cofid(generic_labels)

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.topology, "hcb")
        self.assertEqual(result.cofid, summary.metadata["cofid"])
        self.assertEqual(result.metadata["bond_source"], "distance_inferred")
        self.assertEqual(result.metadata["n_imine_linkage_bonds"], 3)
        detection = result.metadata["topology_detection"]
        self.assertEqual(detection["selected_topology"], "hcb")
        self.assertIn(detection["confidence"], {"exact", "high"})

    def test_generated_decompose_cases_cover_buildable_binary_linkages(self):
        library = ReactionLibrary.builtin()
        buildable = {
            template_id
            for template_id, profile in library.linkage_profiles.items()
            if profile.supports_binary_bridge_pair_generation and profile.supports_atomistic_realization
        }
        covered = {
            template_id
            for template_id, _linkage, _first, _second, _expected in ALL_BUILD_BINARY_LINKAGE_ROUND_TRIP_CASES
        }

        self.assertEqual(covered, buildable)

    def test_decompose_accepts_template_id_linkage_alias(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            summary, _candidate = _generator("hydrazone_bridge").generate_pair_candidate(
                _record("hdz", COF42_HYDRAZIDE, "hydrazide", 2),
                _record("tfb", TFB, "aldehyde", 3),
                out_dir=temp_path,
                write_cif=True,
            )
            input_cif = _without_cofid_comment(summary.cif_path, temp_path / "stripped.cif")

            result = decompose_cif_to_cofid(input_cif, topology="hcb", linkage="hydrazone_bridge")

        self.assertTrue(result.ok, result.reason)
        self.assertEqual(result.linkage, "hydrazone")
        self.assertEqual(result.cofid, summary.metadata["cofid"])

    def test_random_topology_enabled_round_trips_cover_connectivity_combinations(self):
        selected_cases = random.Random(20260625).sample(
            TOPOLOGY_ENABLED_ROUND_TRIP_POOL,
            k=min(8, len(TOPOLOGY_ENABLED_ROUND_TRIP_POOL)),
        )
        observed_connectivity_pairs: set[tuple[int, int]] = set()
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            for case_name, template_id, topology, first, second, allow_repairable_geometry in selected_cases:
                with self.subTest(case=case_name):
                    config_kwargs = {
                        "allowed_reactions": (template_id,),
                        "topology_ids": (topology,),
                        "single_node_topology_ids": (topology,),
                        "enumerate_all_topologies": False,
                        "write_cif": True,
                        "rdkit_num_conformers": 1,
                        "retain_top_results": 1,
                    }
                    if allow_repairable_geometry:
                        config_kwargs.update(
                            {
                                "hard_hard_max_bridge_distance": 10.0,
                                "validation_thresholds": CoarseValidationThresholds(
                                    hard_hard_max_bridge_distance=10.0
                                ),
                            }
                        )
                    generator = BatchStructureGenerator(BatchGenerationConfig(**config_kwargs))
                    summary, _candidate = generator.generate_pair_candidate(
                        first,
                        second,
                        out_dir=temp_path / case_name,
                        write_cif=True,
                    )
                    self.assertEqual(summary.status, "ok")
                    self.assertEqual(summary.topology_id, topology)
                    self.assertIsNotNone(summary.cif_path)
                    input_cif = _without_cofid_comment(
                        summary.cif_path,
                        temp_path / case_name / "stripped.cif",
                    )

                    result = decompose_cif_to_cofid(input_cif, topology=topology, linkage=template_id)

                    self.assertTrue(result.ok, result.reason)
                    self.assertEqual(result.cofid, summary.metadata["cofid"])
                    observed_connectivity_pairs.add(
                        tuple(sorted((first.expected_connectivity, second.expected_connectivity), reverse=True))
                    )

        self.assertGreaterEqual(
            observed_connectivity_pairs,
            {(3, 3), (3, 2), (4, 2), (4, 4), (6, 2)},
        )


if __name__ == "__main__":
    unittest.main()
