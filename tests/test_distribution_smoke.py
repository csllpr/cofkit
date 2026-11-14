"""Portable installation checks: CLI, bundled assets, and a chemistry round trip."""

import subprocess
import sys

from cofkit import (
    CIFWriter, RingFormingStructureGenerator, build_rdkit_monomer,
    decompose_cif_to_cofid, default_topology_repository, packaged_forcefield_catalog,
)


def test_packaged_runtime_assets():
    assert default_topology_repository().load("hcb").node_definitions
    assert packaged_forcefield_catalog()


def test_cli_entrypoint():
    result = subprocess.run([sys.executable, "-m", "cofkit", "--version"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("cofkit ")
    result = subprocess.run([sys.executable, "-m", "cofkit", "--help"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert all(group in result.stdout for group in ("build", "analyze", "calculate", "validate"))


def test_installed_build_export_decomposition_round_trip(tmp_path):
    precursor = build_rdkit_monomer("ctf", "terephthalonitrile", "N#Cc1ccc(C#N)cc1", "nitrile", num_conformers=1)
    candidate = RingFormingStructureGenerator().generate(precursor, "triazine_trimerization", topology_id="hcb")
    cif = tmp_path / "COF with spaces.cif"
    CIFWriter().write_candidate(cif, candidate, {precursor.id: precursor})
    result = decompose_cif_to_cofid(cif, topology="hcb", linkage="triazine")
    assert result.ok, result.to_dict()
    assert result.cofid
