# 原生 Windows Jupyter 教程

语言：[English](../en/README.md) | **中文**

平台：[Bash（Linux/macOS/WSL）](../../cn/README.md) | **原生 Windows PowerShell**

请从仓库根目录启动 Jupyter，并按以下顺序学习：

1. [`01_environment_and_first_build.ipynb`](01_environment_and_first_build.ipynb) 创建或更新 `cofkit` Conda 环境、注册其 Jupyter 内核，并首次构建 TAPB–TFB 亚胺 COF。
2. [`02_topologies_connectivities_and_linkages.ipynb`](02_topologies_connectivities_and_linkages.ipynb) 比较具有代表性的二维/三维拓扑、连接数组合和连接键化学。
3. [`03_zeopp_pore_analysis.ipynb`](03_zeopp_pore_analysis.ipynb) 自动查找 Notebook 1 生成的 CIF，并提供一个简单的 Notebook 内路径输入位置，让学生粘贴可在 Windows 上运行的 Zeo++ 可执行文件路径。

可执行单元格使用 `%%script powershell -NoProfile` 和 `conda run -n cofkit`。相应的 Python API 示例仍保持完全注释状态。Zeo++ 0.3 文档说明可通过 Cygwin 和 GNU Make 在 Windows 上编译；在受限环境中，由教师提供可执行文件是更实际的方案。学生无需在启动 Jupyter 前配置环境变量。
