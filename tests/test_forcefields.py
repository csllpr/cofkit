import copy
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.forcefields import (
    ForceFieldMetadataError,
    load_forcefield_metadata,
    load_packaged_forcefield_metadata,
    packaged_forcefield_metadata_path,
    resolve_forcefield_metadata,
    supported_forcefield_families,
    supported_forcefield_selectors,
    verify_forcefield_artifacts,
)
from cofkit.guest_forcefields import packaged_guest_forcefield_catalog


class ForceFieldMetadataTests(unittest.TestCase):
    def test_packaged_registry_exposes_stable_default_parameter_set_ids(self):
        entries = load_packaged_forcefield_metadata()

        self.assertEqual(supported_forcefield_families(), ("uff", "dreiding"))
        self.assertEqual(
            supported_forcefield_selectors(),
            (
                "uff",
                "uff-openbabel-3.1.0-cofkit-1.0",
                "dreiding",
                "dreiding-standard-1990-cofkit-1.0",
            ),
        )
        self.assertEqual(
            {entry.id for entry in entries},
            {
                "uff-openbabel-3.1.0-cofkit-1.0",
                "dreiding-standard-1990-cofkit-1.0",
            },
        )
        for family in supported_forcefield_families():
            metadata = resolve_forcefield_metadata(family)
            self.assertTrue(metadata.default_for_family)
            self.assertIs(resolve_forcefield_metadata(metadata.id), metadata)
            self.assertEqual(metadata.to_dict()["schema_version"], 1)

    def test_packaged_registry_records_scientific_and_implementation_boundaries(self):
        uff = resolve_forcefield_metadata("uff")
        dreiding = resolve_forcefield_metadata("dreiding")

        for metadata in (uff, dreiding):
            self.assertEqual(metadata.electrostatics.native_charge_model, "none")
            self.assertEqual(
                metadata.electrostatics.supported_charge_sources,
                ("none", "input_cif", "eqeq"),
            )
            self.assertEqual(metadata.nonbonded.epsilon_mixing_rule, "geometric")
            self.assertEqual(metadata.nonbonded.sigma_mixing_rule, "arithmetic")
            self.assertEqual(
                metadata.nonbonded.pair_coefficient_mode_by_backend,
                {"lammps": "explicit_pairij", "raspa": "mixing_rules"},
            )
            self.assertEqual(
                {typing.backends for typing in metadata.atom_typing},
                {("lammps",), ("raspa",)},
            )
            self.assertEqual(metadata.intramolecular_scaling.backends, ("lammps",))
            self.assertEqual(metadata.intramolecular_scaling.lj_14, 1.0)
            self.assertEqual(metadata.intramolecular_scaling.coulomb_14, 1.0)
            self.assertEqual(metadata.coverage.validated_linkages, ())
            self.assertEqual(metadata.coverage.linkage_validation_status, "not_validated")
            self.assertEqual(metadata.license.status, "unreviewed")
            self.assertTrue(metadata.parameter_source.citations[0].doi)

        self.assertIn("Xe", uff.coverage.parameterized_elements)
        self.assertNotIn("Xe", dreiding.coverage.parameterized_elements)
        self.assertEqual(uff.coverage.backend_elements["raspa"], dreiding.coverage.backend_elements["raspa"])

    def test_packaged_parameter_artifacts_are_sha256_verified(self):
        package_root = packaged_forcefield_metadata_path().parents[2]

        for metadata in load_packaged_forcefield_metadata():
            verified = verify_forcefield_artifacts(metadata, artifact_root=package_root)
            self.assertEqual(verified, {"parameter_table": metadata.parameter_source.artifacts[0].sha256})
            artifact = metadata.parameter_source.artifacts[0]
            observed = hashlib.sha256((package_root / artifact.path).read_bytes()).hexdigest()
            self.assertEqual(observed, artifact.sha256)

    def test_registry_rejects_stale_parameter_artifact(self):
        raw = json.loads(packaged_forcefield_metadata_path().read_text(encoding="utf-8"))
        raw["forcefields"] = [raw["forcefields"][0]]
        raw["forcefields"][0]["parameter_source"]["artifacts"][0].update(
            {"path": "artifact.prm", "sha256": "0" * 64}
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "artifact.prm").write_text("changed parameters\n", encoding="utf-8")
            registry_path = root / "registry.json"
            registry_path.write_text(json.dumps(raw), encoding="utf-8")

            with self.assertRaisesRegex(ForceFieldMetadataError, "failed SHA-256 verification"):
                load_forcefield_metadata(registry_path)

    def test_registry_rejects_multiple_defaults_for_one_family(self):
        raw = json.loads(packaged_forcefield_metadata_path().read_text(encoding="utf-8"))
        raw["forcefields"] = [raw["forcefields"][0]]
        duplicate = copy.deepcopy(raw["forcefields"][0])
        duplicate["id"] = "uff-future-variant-2.0"
        duplicate["variant"] = "future-variant"
        raw["forcefields"].append(duplicate)

        with tempfile.TemporaryDirectory() as temp_dir:
            registry_path = Path(temp_dir) / "registry.json"
            registry_path.write_text(json.dumps(raw), encoding="utf-8")

            with self.assertRaisesRegex(ForceFieldMetadataError, "exactly one default variant"):
                load_forcefield_metadata(registry_path, verify_artifacts=False)

    def test_guest_compatibility_accepts_registered_family_or_stable_id(self):
        xenon = packaged_guest_forcefield_catalog()["Xe"]
        carbon_dioxide = packaged_guest_forcefield_catalog()["CO2"]

        self.assertTrue(xenon.supports("uff-openbabel-3.1.0-cofkit-1.0"))
        self.assertTrue(xenon.supports("dreiding-standard-1990-cofkit-1.0"))
        self.assertFalse(carbon_dioxide.supports("uff-openbabel-3.1.0-cofkit-1.0"))


if __name__ == "__main__":
    unittest.main()
