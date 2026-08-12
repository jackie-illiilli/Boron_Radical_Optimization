"""Run a lightweight, deterministic demo of the manuscript modelling workflow.

The demo uses the included final experimental dataset and cached physical-organic
descriptors. It performs five-fold cross-validation with a reduced CatBoost
iteration count so that installation and data flow can be checked quickly on a
standard desktop without rerunning quantum-chemical calculations.
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
from pathlib import Path

import numpy as np
from catboost import CatBoostRegressor
from sklearn.metrics import accuracy_score, f1_score, mean_absolute_error, r2_score
from sklearn.model_selection import KFold

from DFTStructureGenerator import B_N_Cl


ROOT = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--iterations",
        type=int,
        default=500,
        help="CatBoost iterations per fold (default: 500; manuscript analyses use 10,000).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "demo_output" / "metrics.json",
        help="Path for the JSON summary.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.perf_counter()

    with (ROOT / "Data" / "PhysOrgdes_new.pkl").open("rb") as handle:
        descriptor_map = pickle.load(handle)
    with (ROOT / "Data" / "Cldes_new.pkl").open("rb") as handle:
        chloroarene_descriptor_map = pickle.load(handle)

    features, yields = B_N_Cl.descriptor_generator(
        str(ROOT / "Data" / "Iteration" / "Result_sum_00022.xlsx"),
        descriptor_map,
        chloroarene_descriptor_map,
    )

    splitter = KFold(n_splits=5, shuffle=True, random_state=1)
    predictions = np.zeros(len(yields), dtype=float)
    fold_r2 = []
    fold_mae = []

    for train_indices, test_indices in splitter.split(features):
        model = CatBoostRegressor(
            iterations=args.iterations,
            depth=2,
            random_seed=0,
            verbose=False,
            allow_writing_files=False,
            thread_count=1,
        )
        model.fit(features[train_indices], yields[train_indices])
        fold_predictions = model.predict(features[test_indices])
        predictions[test_indices] = fold_predictions
        fold_r2.append(r2_score(yields[test_indices], fold_predictions))
        fold_mae.append(mean_absolute_error(yields[test_indices], fold_predictions))

    observed_active = yields > 20
    predicted_active = predictions > 20
    summary = {
        "dataset": "Data/Iteration/Result_sum_00022.xlsx",
        "samples": int(features.shape[0]),
        "features": int(features.shape[1]),
        "cross_validation": "5-fold shuffled KFold (random_state=1)",
        "model": {
            "name": "CatBoostRegressor",
            "iterations": args.iterations,
            "depth": 2,
            "random_seed": 0,
        },
        "mean_r2": float(np.mean(fold_r2)),
        "mean_mae_percentage_points": float(np.mean(fold_mae)),
        "classification_threshold_percent": 20.0,
        "classification_accuracy": float(accuracy_score(observed_active, predicted_active)),
        "classification_f1": float(f1_score(observed_active, predicted_active)),
        "runtime_seconds": float(time.perf_counter() - started),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
