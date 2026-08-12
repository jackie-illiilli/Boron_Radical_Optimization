# Reactivity-guided hierarchical active learning enables programmable radical dechlorination of polychloroarenes

---

[![DOI](https://zenodo.org/badge/969837101.svg)](https://doi.org/10.5281/zenodo.21897320)

This repository contains the data, notebooks, and helper code used for the manuscript:

**Programmable multistep radical dechlorinative functionalization of polychloroarenes through hierarchical active learning**

The project combines physical-organic reaction descriptors with hierarchical active learning to identify stage-specific conditions for selective C-Cl activation in polychloroarenes.

## Search-space definition

The initially enumerated library contained 222,768 combinations constructed from 26 borane-Lewis base candidates, 14 thiols, 17 solvent systems, 9 initiators, and 4 borane loadings. Before active learning, candidates that were experimentally infeasible because a reagent was unavailable or expected to decompose under the reaction conditions were removed. The curated, experimentally accessible search space therefore contains 205,632 conditions:

```text
24 borane-Lewis base complexes x 14 thiols x 17 solvent systems x 9 initiators x 4 borane loadings
= 205,632 candidate conditions
```

Across the main workflow, CatBoost models and uncertainty-guided acquisition were used to navigate this curated space with 270 experiments and identify complementary condition sets for meta-, ortho-, and para-selective dechlorinative functionalization.

## Repository layout

```text
.
|-- demo.py                          # Fast end-to-end modelling demo
|-- Bayesian Optimization.ipynb      # Active-learning selection and optimization figures
|-- Calculate PhysOrg.ipynb          # Descriptor provenance; DFT reruns are optional
|-- Modelling&Validation.ipynb       # Validation, classification, importance, and PySR
|-- DFTStructureGenerator/           # Local DFT/structure/descriptor utilities
|-- Data/
|   |-- Reactants.csv                # Original reagent list
|   |-- Processed_Reactants.csv      # Curated reagents and descriptor metadata
|   |-- PhysOrgdes_new.pkl           # Physical-organic reagent descriptors
|   |-- Cldes_new.pkl                # Chloroarene/site descriptors
|   |-- Fingerprint.pkl              # Fingerprint baseline descriptors
|   |-- Iteration/                   # Hierarchical active-learning history
|   `-- Iteration2/                  # Transfer-learning experiments
|-- Figure/                          # Generated manuscript-facing figures
|-- docs/
|   |-- ALGORITHM.md                 # Detailed algorithm and pseudocode
|   |-- CODE_CHECKLIST_STATUS.md     # Nature checklist status and remaining actions
|   |-- Manuscript dechlorination.docx
|   `-- README_CN.md                 # Chinese repository guide
|-- environment.yml                  # Reproducible conda environment
|-- requirements.txt                 # Core pip dependencies
|-- requirements-analysis-optional.txt # Optional PySR dependency
`-- requirements-dft-optional.txt    # Optional descriptor-regeneration dependencies
```

## System requirements

### Core workflow

- Python 3.12; the environment is named `main_py3_12`.
- Dependency versions and compatibility ranges are specified in `environment.yml` and `requirements.txt`.
- No GPU or other non-standard hardware is required. The demo and model-validation workflow run on a conventional 64-bit desktop CPU; 8 GB RAM is recommended.
- The code has been tested on Windows 10 Pro 64-bit, build 19045. A fresh installation from `environment.yml` was verified with Python 3.12.13, NumPy 1.26.4, pandas 2.3.3, openpyxl 3.1.5, scikit-learn 1.6.1, CatBoost 1.2.7, RDKit 2024.09.6, Matplotlib 3.10.9, seaborn 0.13.2, tqdm 4.70.0, JupyterLab 4.6.3, Notebook 7.6.2, and ipykernel 6.31.0.
- Linux and macOS have not yet been independently tested.

### Optional workflows

- Symbolic regression requires `requirements-analysis-optional.txt` and a working Julia installation managed by PySR. This is not required for the demo or the core CatBoost results.
- Descriptor regeneration requires `requirements-dft-optional.txt`, xTB/CREST, Gaussian, and access to the original quantum-chemistry files and paths. External quantum-chemistry executables are not bundled or licensed by this repository.
- The default manuscript-level reproduction starts from the included descriptor caches and does not require xTB, CREST, Gaussian, Julia, a GPU, or a computing cluster.

## Installation guide

Clone or unpack the repository, open a terminal in its root directory, and create the conda environment:

```bash
conda env create -f environment.yml
conda activate main_py3_12
```

If an environment with that name already exists, update it instead:

```bash
conda env update -n main_py3_12 -f environment.yml --prune
conda activate main_py3_12
```

Alternatively, with an existing Python 3.12 installation:

```bash
python -m pip install -r requirements.txt
```

On a normal desktop with a broadband connection, the core conda installation is expected to take approximately 10-20 minutes. A clean Windows test took about 6 minutes; package-download speed and cache state are the main sources of variation. Installation does not compile project-specific native code.

Start notebooks from the repository root so that relative paths resolve correctly:

```bash
jupyter lab
```

## Demo

The demo loads the final experimental results and cached descriptors, constructs the same 82-column physical-organic representation used by the modelling workflow, and performs deterministic five-fold CatBoost cross-validation with a reduced iteration count.

Run:

```bash
python demo.py
```

Input data:

- `Data/Iteration/Result_sum_00022.xlsx`
- `Data/PhysOrgdes_new.pkl`
- `Data/Cldes_new.pkl`

Expected terminal output and `demo_output/metrics.json` contents, allowing for small platform-level floating-point differences:

```text
samples: 257
features: 82
mean_r2: approximately 0.692
mean_mae_percentage_points: approximately 9.88
classification_accuracy at the 20% yield threshold: approximately 0.864
classification_f1: approximately 0.900
```

The source workbook contains 270 experimental rows; 257 complete cases remain after the same `dropna()` preprocessing used by the notebooks. On the tested desktop (Intel Core i5-13600K), the demo required about 6 seconds including conda process startup and about 1.5 seconds inside Python.

To reproduce the manuscript-scale 10,000-iteration validation setting with the same demo entry point:

```bash
python demo.py --iterations 10000 --output demo_output/metrics_10000.json
```

The tested result was mean R2 = 0.7009989 and mean MAE = 9.6753 percentage points; the end-to-end runtime was about 29 seconds on the test machine.

## Running the workflow on compatible data

For data that use the included reagent and chloroarene descriptor caches, prepare an `.xlsx` file with one experimental condition per row and the following columns:

| Column        | Meaning                                                 | Requirement                            |
| ------------- | ------------------------------------------------------- | -------------------------------------- |
| `B_Index`   | Borane-Lewis base index from`Processed_Reactants.csv` | Required                               |
| `eqs`       | Borane loading in equivalents                           | Required                               |
| `S_Index`   | Thiol index                                             | Required                               |
| `ini_Index` | Initiator index                                         | Required                               |
| `sol_Index` | Solvent-system index                                    | Required                               |
| `Cl_Index`  | Chloroarene descriptor index                            | Required for manuscript substrates     |
| `Cl_Atomid` | Target chlorine-site atom index                         | Required for site-specific entries     |
| `yield`     | Experimental yield in percent                           | Required for model training/validation |
| `sol_name`  | Human-readable solvent label                            | Optional; ignored by the model         |

All categorical indices must exist in `Data/Processed_Reactants.csv`, `Data/PhysOrgdes_new.pkl`, or `Data/Cldes_new.pkl`. Then construct the feature matrix with:

```python
import pickle
from DFTStructureGenerator import B_N_Cl

with open("Data/PhysOrgdes_new.pkl", "rb") as handle:
    reagent_descriptors = pickle.load(handle)
with open("Data/Cldes_new.pkl", "rb") as handle:
    chloroarene_descriptors = pickle.load(handle)

X, y = B_N_Cl.descriptor_generator(
    "path/to/your_results.xlsx",
    reagent_descriptors,
    chloroarene_descriptors,
)
```

New reagents or substrates that are absent from the caches require descriptor generation. That advanced path depends on external quantum-chemistry software and must be configured for the user's local installation; it is documented for provenance in `Calculate PhysOrg.ipynb` but is not part of the standard reproduction workflow.

## Manuscript-result reproduction

### 1. Model validation and interpretation

Open `Modelling&Validation.ipynb` and run the core CatBoost sections using `Data/Iteration/Result_sum_00022.xlsx`, `Data/PhysOrgdes_new.pkl`, and `Data/Cldes_new.pkl`. The fixed settings are five shuffled folds (`random_state=1`) and CatBoost models with 10,000 iterations, depth 2, and `random_seed=0`.

Expected core regression results are mean R2 approximately 0.701 and mean MAE approximately 9.68 percentage points. Subsequent cells generate:

- `Data/Correlation_matrix.csv`
- `Data/Feature_Importance.csv`
- `Figure/Feature_Importance.png`

The PySR cells are optional and require the separate analysis dependency plus Julia.

### 2. Hierarchical active-learning recommendation

Open `Bayesian Optimization.ipynb`. The notebook reconstructs all 205,632 feasible candidate conditions, loads the cached descriptors, trains the uncertainty-aware CatBoost model, and selects a diverse batch of ten conditions. Historical main-search results are stored in `Data/Iteration/`; transfer-learning results are stored in `Data/Iteration2/`.

Before generating a recommendation, set `Times`, `Cl`, and `Cl_atomid` to the intended historical or new stage. Recommendation cells write new picker arrays and CSV files, so use a copy of the repository when replaying earlier iterations. The active-learning pseudocode and parameter definitions are in `docs/ALGORITHM.md`.

### 3. Descriptor provenance

`Calculate PhysOrg.ipynb` documents how reactant descriptors and fingerprint caches were generated. For standard reproduction, run only the cache-based analysis paths. Do not run the xTB/Gaussian sections unless the external programs, source calculation files, and local paths have been configured.

## Data notes

- `Data/Iteration/first_Data.csv` contains the initial 50 randomly selected conditions.
- `Data/Iteration/data_*.csv` contain recommended experimental batches from the main hierarchical search.
- `Data/Iteration/Result_sum_*.xlsx` contain cumulative experimental results used for model fitting.
- `Data/Iteration/drop_n.pkl` and `drop_p.pkl` store masks for excluding stronger reactivity windows.
- `Data/Iteration2/` stores transfer-learning batches for a shifted substrate space.

## Code and data availability

An immutable release of the source code, notebooks, cached descriptors, and supporting data is archived on Zenodo at [https://doi.org/10.5281/zenodo.21897320](https://doi.org/10.5281/zenodo.21897320). The GitHub repository remains private during manuscript preparation; its public URL will be added here and to the manuscript before publication.

The source code is released under the OSI-approved BSD 3-Clause License. See `LICENSE` for the complete terms.

## Version-control hygiene

Generated runtime folders such as `demo_output/`, `catboost_info/`, `outputs/`, Python bytecode caches, office lock files, and scratch outputs are ignored by git. Source code, notebooks, manuscript-facing data, descriptor caches, and figures required for reproduction should remain version controlled.
