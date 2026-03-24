from pathlib import Path

from cofkit import COFEngine, write_candidate_cif

from mock_imine import build_project


if __name__ == "__main__":
    project = build_project()
    best = COFEngine().run(project).top(1)[0]
    out_path = Path("out/mock_imine.cif")
    result = write_candidate_cif(out_path, best, project.monomers, data_name="mock_imine")
    print("wrote:", out_path)
    print("mode:", result.mode)
    print("sites:", result.n_sites)
