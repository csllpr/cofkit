# 原生 Windows Jupyter 教程

语言：[English](../en/README.md) | **中文**

平台：[Bash（Linux/macOS/WSL）](../../cn/README.md) | **使用 Python 单元格的原生 Windows**

请从仓库根目录启动 Jupyter，并按以下顺序学习：

1. [`01_environment_setup.ipynb`](01_environment_setup.ipynb) 创建或更新 `cofkit` Conda 环境，并注册其 Jupyter 内核。
2. [`02_first_cof_build.ipynb`](02_first_cof_build.ipynb) 首次构建 TAPB–TFB 亚胺 COF。
3. [`03_topologies_connectivities_and_linkages.ipynb`](03_topologies_connectivities_and_linkages.ipynb) 比较具有代表性的二维/三维拓扑、连接数组合和连接键化学。
4. [`04_zeopp_pore_analysis.ipynb`](04_zeopp_pore_analysis.ipynb) 自动查找 Notebook 2 生成的 CIF，并提供一个简单的 Notebook 内路径输入位置，让学生粘贴可在 Windows 上运行的 Zeo++ 可执行文件路径。

可执行单元格均为普通 Python。它们通过 `subprocess` 安全传递 `conda run -n cofkit` 和 cofkit CLI 参数，不使用 PowerShell 单元格魔法或 Shell 引号。相应的直接 Python API 示例仍保持完全注释状态。Zeo++ 0.3 文档说明可通过 Cygwin 和 GNU Make 在 Windows 上编译；在受限环境中，由教师提供可执行文件是更实际的方案。学生无需在启动 Jupyter 前配置环境变量。
