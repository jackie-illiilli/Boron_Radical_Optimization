# Programmable Multistep Radical Dechlorinative Functionalization

This repository contains the data, notebooks, and helper code used for the manuscript:

**Programmable multistep radical dechlorinative functionalization of polychloroarenes through hierarchical active learning**

The project combines physical-organic reaction descriptors with hierarchical active learning to identify stage-specific conditions for selective C-Cl activation in polychloroarenes. The search space contains 205,632 candidate conditions:

```text
24 borane-Lewis base complexes x 14 thiols x 17 solvent systems x 9 initiators x 4 borane loadings
```

Across the main workflow, CatBoost models and uncertainty-guided acquisition were used to navigate this space with approximately 270 experiments and identify orthogonal condition sets for meta-, ortho-, and para-selective dechlorinative functionalization.

## Repository Layout

```text
.
|-- Bayesian Optimization.ipynb      # Active-learning condition selection and optimization figures
|-- Calculate PhysOrg.ipynb          # Descriptor construction workflow; DFT reruns are optional
|-- Modelling&Validation.ipynb       # Cross-validation, classification, feature importance, PySR analysis
|-- DFTStructureGenerator/           # Local DFT/structure/descriptor utility toolkit
|-- Data/
|   |-- Reactants.csv                # Candidate boranes, thiols, solvents, and initiators
|   |-- Processed_Reactants.csv      # Canonicalized reactants with atom IDs and descriptor metadata
|   |-- PhysOrgdes_new.pkl           # Physical-organic descriptors used for modelling
|   |-- Cldes_new.pkl                # Chloroarene/site descriptors used for staged searches
|   |-- Fingerprint.pkl              # Fingerprint baseline descriptors
|   |-- Iteration/                   # Main hierarchical active-learning rounds
|   `-- Iteration2/                  # Transfer-learning rounds for shifted substrate space
|-- Figure/                          # Generated figure panels used in the manuscript
|-- docs/
|   |-- Manuscript dechlorination-0511.docx
|   `-- README_CN.md                 # Chinese repository guide
|-- environment.yml                  # Conda environment for the reproducible workflow
|-- requirements.txt                 # Core Python requirements
|-- requirements-analysis-optional.txt # Optional PySR dependency for symbolic regression
`-- requirements-dft-optional.txt    # Optional dependencies for descriptor regeneration
```

## Reproducibility Scope

The standard reproducible workflow starts from the included descriptor caches and experimental result files:

- `Data/PhysOrgdes_new.pkl`
- `Data/Cldes_new.pkl`
- `Data/Fingerprint.pkl`
- `Data/Iteration/*.xlsx`
- `Data/Iteration2/*.xlsx`

The raw xTB/Gaussian calculations do **not** need to be rerun for manuscript-level model and figure reproduction. `DFTStructureGenerator/` and the early sections of `Calculate PhysOrg.ipynb` document the internal workflow used to build the descriptor caches, but those sections require external quantum-chemistry software and local paths that are not part of the default reproduction path.

## Environment

The working environment used for this project is `main_py3_12` with Python 3.12.

To recreate the core analysis environment with conda:

```bash
conda env create -f environment.yml
conda activate main_py3_12
```

If you already have Python 3.12 available:

```bash
pip install -r requirements.txt
```

Optional descriptor-regeneration dependencies are listed separately:

```bash
pip install -r requirements-dft-optional.txt
```

`PySR` is only required for the symbolic-regression section of `Modelling&Validation.ipynb` and may require a working Julia installation:

```bash
pip install -r requirements-analysis-optional.txt
```

## Recommended Workflow

Run notebooks from the repository root so relative paths resolve correctly.

1. **Model validation and interpretation**

   Open `Modelling&Validation.ipynb` and run the sections using `Data/PhysOrgdes_new.pkl`, `Data/Cldes_new.pkl`, and `Data/Iteration/Result_sum_00022.xlsx`. This reproduces the CatBoost regression/classification validation and feature-importance analysis from the final main optimization dataset.

2. **Hierarchical active-learning recommendation**

   Open `Bayesian Optimization.ipynb`. The notebook builds the full 205,632-condition space, loads descriptor caches, trains CatBoost uncertainty models, and generates the next recommended batch for a selected iteration. The main search history is stored in `Data/Iteration/`; transfer-learning runs are stored in `Data/Iteration2/`.

3. **Descriptor cache provenance**

   `Calculate PhysOrg.ipynb` records how reactant descriptors were generated from structure optimization and DFT parsing. For normal reproduction, use the included `.pkl` descriptor caches instead of rerunning DFT.

## Data Notes

- `Data/Iteration/first_Data.csv` contains the initial 50 randomly selected experiments.
- `Data/Iteration/data_*.csv` contain recommended experimental batches from the main hierarchical search.
- `Data/Iteration/Result_sum_*.xlsx` contain accumulated experimental results used for model fitting.
- `Data/Iteration/drop_n.pkl` and `Data/Iteration/drop_p.pkl` store filtered condition masks used to exclude stronger reactivity windows.
- `Data/Iteration2/` stores transfer-learning batches used for reactivity-shifted substrate space.

## Local Utility Code

`DFTStructureGenerator/` is a local utility toolkit used by the author across DFT projects. In this repository it supports:

- reactant preprocessing and atom-index extraction;
- xTB/Gaussian input generation;
- Gaussian log parsing and energy collection;
- construction of physical-organic descriptors for modelling.

Only a subset of this toolkit is required for the default notebook workflow because descriptor caches are already included.

## Version-Control Hygiene

Generated runtime folders such as `catboost_info/`, `outputs/`, Python bytecode caches, office lock files, and scratch outputs are ignored by git. The tracked repository should contain source code, notebooks, manuscript-facing data, descriptor caches, and figure assets needed for reproduction.
