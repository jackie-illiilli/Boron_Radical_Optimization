# 仓库说明（中文）

本仓库是论文 **Programmable multistep radical dechlorinative functionalization of polychloroarenes through hierarchical active learning** 的配套代码和数据仓库。

项目目标是用分层主动学习寻找多氯芳烃逐步、位点选择性脱氯官能化的反应条件。条件空间由以下组合构成：

```text
24 个硼烷-路易斯碱复合物 x 14 个硫醇 x 17 个溶剂体系 x 9 个引发剂 x 4 个硼烷当量 = 205,632 个候选条件
```

核心思路是把 C-Cl 键的内在反应性层级转化为主动学习的搜索顺序：先寻找较强的 meta 位活化条件，再排除强反应窗口，继续寻找 ortho 和 para 位所需的更弱、更选择性的条件。

## 主要文件

```text
Bayesian Optimization.ipynb
    主动学习与条件推荐 notebook。负责构建完整条件空间、训练 CatBoost 模型、计算不确定性，并输出下一轮推荐实验。

Modelling&Validation.ipynb
    模型验证与解释性分析 notebook。包括交叉验证、分类阈值分析、特征重要性和 PySR 分析。

Calculate PhysOrg.ipynb
    物理有机描述符生成流程记录。默认复现不需要重跑 DFT，只需要使用已有 pkl 缓存。

DFTStructureGenerator/
    作者常用的本地 DFT/结构/描述符工具箱，不是本项目专用包。这里主要用于反应物预处理、Gaussian/xTB 输入生成、log 解析和描述符构建。

Data/
    反应物、描述符缓存、主动学习迭代数据和迁移学习数据。

Figure/
    论文相关图表输出。
```

## 默认复现范围

默认复现从已有描述符和实验结果开始，不要求重新运行 xTB/Gaussian。

需要使用的核心缓存包括：

- `Data/PhysOrgdes_new.pkl`
- `Data/Cldes_new.pkl`
- `Data/Fingerprint.pkl`
- `Data/Iteration/*.xlsx`
- `Data/Iteration2/*.xlsx`

`Calculate PhysOrg.ipynb` 中的 DFT 部分主要用于记录描述符来源。如果需要完全重建描述符，需要额外准备 Gaussian、xTB/CREST、morfeus 等依赖，并修改 notebook 中的本地路径。

## 环境

本项目使用的环境是 `main_py3_12`，Python 版本为 3.12。

如果用 conda 重建环境：

```bash
conda env create -f environment.yml
conda activate main_py3_12
```

如果已有 Python 3.12 环境：

```bash
pip install -r requirements.txt
```

DFT 描述符重建的可选依赖放在：

```bash
requirements-dft-optional.txt
```

`PySR` 只用于 `Modelling&Validation.ipynb` 中的符号回归部分，通常还需要 Julia 环境：

```bash
pip install -r requirements-analysis-optional.txt
```

## 推荐运行顺序

请从仓库根目录启动 Jupyter，保证相对路径正常。

1. 运行 `Modelling&Validation.ipynb`

   用 `Data/PhysOrgdes_new.pkl`、`Data/Cldes_new.pkl` 和 `Data/Iteration/Result_sum_00022.xlsx` 复现最终模型验证、分类表现和特征重要性分析。

2. 运行 `Bayesian Optimization.ipynb`

   用已有描述符构建 205,632 个候选条件空间，并根据设定的 `Times`、`Cl_Index` 和 `Cl_Atomid` 复现主优化或迁移学习的下一轮推荐。

3. 查看 `Calculate PhysOrg.ipynb`

   用于理解描述符缓存如何由结构优化、DFT 结果和物理有机参数构建。默认不建议从头重跑。

## 数据目录说明

- `Data/Reactants.csv`：候选 B-N 复合物、硫醇、溶剂和引发剂清单。
- `Data/Processed_Reactants.csv`：规范化后的反应物、索引、关键原子编号和 DFT 能量。
- `Data/Iteration/first_Data.csv`：初始 50 个随机实验。
- `Data/Iteration/data_*.csv`：主分层主动学习每轮推荐实验。
- `Data/Iteration/Result_sum_*.xlsx`：每轮累计实验结果，用于训练和验证模型。
- `Data/Iteration/drop_n.pkl` 与 `drop_p.pkl`：排除强反应窗口的条件掩码。
- `Data/Iteration2/`：迁移学习相关的推荐批次和结果。

## 版本控制说明

仓库应保留源码、notebook、论文复现所需数据、描述符缓存和图表。运行时自动生成的 `catboost_info/`、`outputs/`、`__pycache__/`、Office 锁文件和临时测试输出不再纳入版本控制。
