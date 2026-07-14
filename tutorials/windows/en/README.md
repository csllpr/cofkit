# Native Windows Jupyter tutorials

Language: **English** | [中文](../cn/README.md)

Platform: [Bash (Linux/macOS/WSL)](../../en/README.md) | **Native Windows with Python cells**

Run Jupyter from the repository root and follow the notebooks in order:

1. [`01_environment_setup.ipynb`](01_environment_setup.ipynb) creates or updates the `cofkit` Conda environment and registers its Jupyter kernel.
2. [`02_first_cof_build.ipynb`](02_first_cof_build.ipynb) builds the first TAPB–TFB imine COF.
3. [`03_topologies_connectivities_and_linkages.ipynb`](03_topologies_connectivities_and_linkages.ipynb) compares selected 2D/3D topologies, connection combinations, and linkage chemistries.
4. [`04_zeopp_pore_analysis.ipynb`](04_zeopp_pore_analysis.ipynb) discovers Notebook 2's CIF and provides a simple in-notebook field where students paste the path to a Windows-runnable Zeo++ executable.

The executable cells are regular Python. They pass `conda run -n cofkit` and cofkit CLI arguments safely through `subprocess`, with no PowerShell cell magic or shell quoting. Corresponding direct Python API examples remain fully commented. Zeo++ 0.3 documents Windows compilation through Cygwin and GNU Make; in restricted environments, an instructor-provided executable is the practical option. Students do not need to configure environment variables before launching Jupyter.
