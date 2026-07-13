# Native Windows Jupyter tutorials

Language: **English** | [中文](../cn/README.md)

Platform: [Bash (Linux/macOS/WSL)](../../en/README.md) | **Native Windows PowerShell**

Run Jupyter from the repository root and follow the notebooks in order:

1. [`01_environment_and_first_build.ipynb`](01_environment_and_first_build.ipynb) creates or updates the `cofkit` Conda environment, registers its kernel, and builds a first TAPB–TFB imine COF.
2. [`02_topologies_connectivities_and_linkages.ipynb`](02_topologies_connectivities_and_linkages.ipynb) compares selected 2D/3D topologies, connection combinations, and linkage chemistries.
3. [`03_zeopp_pore_analysis.ipynb`](03_zeopp_pore_analysis.ipynb) discovers Notebook 1's CIF and provides a simple in-notebook field where students paste the path to a Windows-runnable Zeo++ executable.

The executable cells use `%%script powershell -NoProfile` and `conda run -n cofkit`. Corresponding Python API examples remain fully commented. Zeo++ 0.3 documents Windows compilation through Cygwin and GNU Make; in restricted environments, an instructor-provided executable is the practical option. Students do not need to configure environment variables before launching Jupyter.
