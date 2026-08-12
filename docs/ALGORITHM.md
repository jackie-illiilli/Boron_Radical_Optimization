# Hierarchical Active-Learning Algorithm

This document describes the computational workflow implemented in `Bayesian Optimization.ipynb`, `Modelling&Validation.ipynb`, and `DFTStructureGenerator/B_N_Cl.py`.

## Inputs and feasible search space

The initial combinatorial enumeration comprised 222,768 conditions from 26 borane-Lewis base candidates, 14 thiols, 17 solvent systems, 9 initiators, and 4 borane loadings. Candidates involving unavailable reagents or reagents expected to decompose under the reaction conditions were removed before model-guided experimentation. The curated space used by the code contains:

```text
24 borane-Lewis base complexes
14 thiols
17 solvent systems
9 initiators
4 borane loadings (1.5, 2.0, 2.5, and 3.0 equivalents)
Total: 24 x 14 x 17 x 9 x 4 = 205,632 conditions
```

Each condition is identified by `B_Index`, `S_Index`, `ini_Index`, `sol_Index`, and `eqs`. A substrate and target C-Cl site are identified by `Cl_Index` and `Cl_Atomid`.

## Reaction representation

`B_N_Cl.descriptor_generator` concatenates the cached descriptor vectors in the following order:

1. target chloroarene/C-Cl-site descriptors;
2. borane-Lewis base descriptors;
3. initiator descriptors;
4. solvent descriptors;
5. thiol descriptors;
6. borane equivalents.

For the manuscript dataset this produces an 82-column feature matrix. Rows containing missing experimental fields are removed using `dropna()` before descriptor construction.

## Model and validation settings

- Regression: `CatBoostRegressor(iterations=10000, depth=2, random_seed=0)`.
- Uncertainty prediction: the CatBoost `RMSEWithUncertainty` loss returns a predicted mean and uncertainty for every candidate.
- Reactivity classification: yields above 20% are labelled active; CatBoost classification defines the exclusion mask for the next, milder stage.
- Validation: five-fold shuffled `KFold` with `random_state=1`.
- Reported regression metrics: R2 and mean absolute error (MAE).
- Threshold metrics: accuracy, F1 score, and the confusion matrix at the 20% yield cutoff.

## Batch acquisition

Candidate priority combines predicted yield and uncertainty, while conditions outside the current admissible subspace receive an availability penalty. Batch diversity is enforced with RDKit's `MaxMinPicker`, using similarities calculated across the borane, thiol, initiator, solvent, and loading components. The initial dataset contains 50 randomly selected conditions; subsequent batches contain 10 conditions.

## Pseudocode

```text
enumerate all 222,768 nominal combinations
remove experimentally infeasible combinations
feasible_space <- 205,632 conditions

experimental_data <- 50 initial random conditions
admissible_space <- feasible_space

for stage in [strongest activation, intermediate activation, mildest activation]:
    repeat until the stage stopping criterion is met:
        X_train, y_train <- descriptors and yields from all available experiments
        fit uncertainty-aware CatBoost regressor with random_seed = 0
        mean_yield, uncertainty <- predict every condition in admissible_space
        acquisition <- mean_yield + uncertainty
        select 10 high-acquisition, compositionally diverse conditions
        perform experiments and append measured yields to experimental_data

    active <- measured_yield > 20%
    fit CatBoost classifier for the stage activation boundary
    predicted_active <- classifier(feasible_space)
    admissible_space <- admissible_space minus predicted_active

return stage-matched condition sets and all experimental/model histories
```

The computational search proceeds from the strongest activation requirement toward milder windows so that overly reactive conditions identified at one stage can be excluded before the next stage. The experimentally executed synthetic sequence can then apply the matched conditions in the required order for programmable functionalization.

## Outputs

The main workflow records:

- cumulative experimental workbooks: `Data/Iteration/Result_sum_*.xlsx`;
- recommended batches: `Data/Iteration/data_*.csv`;
- selected-condition indices: `Data/Iteration/picker_*.npy`;
- stage-exclusion masks: `Data/Iteration/drop_n.pkl` and `drop_p.pkl`;
- transfer-learning histories: `Data/Iteration2/`;
- validation and interpretation outputs in `Data/` and `Figure/`.

All random seeds used by the manuscript-facing CatBoost and cross-validation code are stated above. The notebooks retain the complete historical iteration files so that reported learning trajectories can be inspected without rerunning experiments.
