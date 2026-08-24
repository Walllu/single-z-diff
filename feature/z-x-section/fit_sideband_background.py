#!/usr/bin/env python3
"""Fit a provisional smooth nonprompt mass model in prompt-subtracted sidebands."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from scipy.optimize import minimize
import uproot

import plot_abcd_diagnostics as diagnostics


HERE = Path(__file__).resolve().parent
EDGES = np.arange(60.0, 121.0, 1.0)
CENTERS = 0.5 * (EDGES[:-1] + EDGES[1:])
LOW = (CENTERS >= 60.0) & (CENTERS < 75.0)
SIGNAL = (CENTERS >= 80.0) & (CENTERS < 100.0)
HIGH = (CENTERS >= 105.0) & (CENTERS < 120.0)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--label", default="local_functional_background")
    parser.add_argument("--prompt-scale", type=float)
    parser.add_argument("--isolated-max", type=float, default=0.15)
    parser.add_argument("--anti-isolated-min", type=float, default=0.25)
    return parser.parse_args()


def load_corrected_selected(directory: Path, sample: str) -> dict[str, np.ndarray]:
    filename = directory / f"{sample}_selected.root"
    weight_branch = "event_weight_nominal" if sample == "mc" else None
    branches = ["dimuon_mass"] + ([weight_branch] if weight_branch else [])
    with uproot.open(f"{filename}:Events") as tree:
        arrays = tree.arrays(branches, library="np")
    result = {"dimuon_mass": np.asarray(arrays["dimuon_mass"])}
    result["weight"] = (
        np.ones_like(result["dimuon_mass"], dtype=float)
        if weight_branch is None else np.asarray(arrays[weight_branch], dtype=float)
    )
    return result


def hist(sample: dict[str, np.ndarray], mask: np.ndarray,
         weight: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    selected = None if weight is None else weight[mask]
    value, _ = np.histogram(sample["dimuon_mass"][mask], bins=EDGES, weights=selected)
    variance, _ = np.histogram(
        sample["dimuon_mass"][mask], bins=EDGES,
        weights=None if selected is None else np.square(selected),
    )
    if selected is None:
        variance = value.astype(float)
    return value.astype(float), variance.astype(float)


def background(params: np.ndarray) -> np.ndarray:
    log_events_per_gev, slope = params
    return np.exp(np.clip(log_events_per_gev + slope * (CENTERS - 90.0), -50.0, 50.0))


def fit(data: np.ndarray, prompt: np.ndarray, mask: np.ndarray) -> dict[str, Any]:
    residual = max(1.0, float(np.maximum(data[mask] - prompt[mask], 0.0).sum()))
    initial = np.asarray([math.log(residual / max(1, int(mask.sum()))), 0.0])

    def objective(params: np.ndarray) -> float:
        expected = np.maximum(prompt + background(params), 1.0e-12)
        observed = data[mask]
        mu = expected[mask]
        return float(2.0 * np.sum(mu - observed + np.where(observed > 0.0, observed * np.log(observed / mu), 0.0)))

    result = minimize(objective, initial, method="L-BFGS-B",
                      bounds=((-20.0, 20.0), (-0.5, 0.5)))
    if not result.success:
        raise RuntimeError(f"Sideband fit failed: {result.message}")
    model = background(result.x)
    return {
        "parameters": {"log_events_per_gev_at_90": float(result.x[0]),
                       "slope_per_gev": float(result.x[1])},
        "deviance": float(result.fun),
        "bins": int(mask.sum()),
        "background": model,
        "background_signal_yield": float(model[SIGNAL].sum()),
    }


def region_deviance(data: np.ndarray, expected: np.ndarray, mask: np.ndarray) -> float:
    observed, mu = data[mask], np.maximum(expected[mask], 1.0e-12)
    return float(2.0 * np.sum(mu - observed + np.where(observed > 0.0, observed * np.log(observed / mu), 0.0)))


def main() -> int:
    args = arguments()
    corrected_input = ((args.selection_dir / "data_selected.root").is_file()
                       and (args.selection_dir / "mc_selected.root").is_file())
    if corrected_input:
        selected_data = load_corrected_selected(args.selection_dir, "data")
        selected_mc = load_corrected_selected(args.selection_dir, "mc")
        data_a, _ = np.histogram(selected_data["dimuon_mass"], bins=EDGES)
        data_a = data_a.astype(float)
        data_var = data_a.copy()
        mc_a, _ = np.histogram(selected_mc["dimuon_mass"], bins=EDGES,
                               weights=selected_mc["weight"])
        mc_var, _ = np.histogram(selected_mc["dimuon_mass"], bins=EDGES,
                                 weights=np.square(selected_mc["weight"]))
        data_b = data_c = data_d = np.zeros_like(data_a)
    else:
        data = diagnostics.load_skim(args.selection_dir / "data_abcd_skim.root")
        mc = diagnostics.load_skim(args.selection_dir / "mc_abcd_skim.root")
        dc = diagnostics.categories(data, "both", args.isolated_max, args.anti_isolated_min)
        pc = diagnostics.categories(mc, "both", args.isolated_max, args.anti_isolated_min)
        data_a, data_var = hist(data, dc["A"])
        mc_a, mc_var = hist(mc, pc["A"], mc["event_weight"])
        data_b, _ = hist(data, dc["B"])
        data_c, _ = hist(data, dc["C"])
        data_d, _ = hist(data, dc["D"])
    if args.prompt_scale is None:
        alpha = float(data_a[SIGNAL].sum() / mc_a[SIGNAL].sum())
        normalization_method = "provisional OS-isolated 80-100 GeV shape normalization"
    else:
        alpha = float(args.prompt_scale)
        normalization_method = "fixed externally supplied prompt scale"
    prompt = alpha * mc_a
    combined = fit(data_a, prompt, LOW | HIGH)
    low_fit = fit(data_a, prompt, LOW)
    high_fit = fit(data_a, prompt, HIGH)
    expected = prompt + combined["background"]
    low_expected = prompt + low_fit["background"]
    high_expected = prompt + high_fit["background"]

    b, c, d = data_b[SIGNAL].sum(), data_c[SIGNAL].sum(), data_d[SIGNAL].sum()
    abcd = None if d <= 0.0 else float(b * c / d)
    output = args.output_dir / args.label
    output.mkdir(parents=True, exist_ok=True)
    hep.style.use("CMS")
    fig, (upper, lower) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06},
    )
    upper.errorbar(CENTERS, data_a, yerr=np.sqrt(data_var), fmt="o", color="black",
                   markersize=3.5, label="OS isolated data")
    hep.histplot(prompt, EDGES, ax=upper, histtype="step", color="#3b82f6",
                 linestyle="--", label="Scaled prompt DY")
    hep.histplot(combined["background"], EDGES, ax=upper, histtype="step",
                 color="#f59e0b", label="Exponential background")
    hep.histplot(expected, EDGES, ax=upper, histtype="step", color="#ef4444",
                 linewidth=1.8, label="Prompt + background")
    pull = (data_a - expected) / np.sqrt(np.maximum(data_var + alpha * alpha * mc_var, 1.0))
    lower.axhspan(-1, 1, color="#10b981", alpha=0.12)
    lower.axhline(0.0, color="black", linewidth=1.0)
    lower.plot(CENTERS, pull, "o", color="#374151", markersize=3.2)
    for boundary in (75.0, 80.0, 100.0, 105.0):
        upper.axvline(boundary, color="0.5", linestyle="--", linewidth=0.9)
        lower.axvline(boundary, color="0.5", linestyle="--", linewidth=0.9)
    upper.set_yscale("log")
    upper.set_ylim(bottom=0.5)
    upper.set_ylabel("Events / GeV")
    upper.legend(fontsize=10, ncol=2)
    upper.text(0.03, 0.70, "Fit regions: 60-75 and 105-120 GeV\n"
               f"Signal interpolation: {combined['background_signal_yield']:.1f} events",
               transform=upper.transAxes, fontsize=10)
    lower.set_ylabel("Pull")
    lower.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    lower.set_ylim(-5.0, 5.0)
    hep.cms.label("Open Data", data=True, com=13, ax=upper)
    fig.savefig(output / "sideband_exponential_interpolation.svg", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "title": "Provisional prompt-subtracted exponential sideband interpolation",
        "selection_directory": str(args.selection_dir.resolve()),
        "corrected_selected_input": corrected_input,
        "prompt_scale": alpha,
        "prompt_normalization_method": normalization_method,
        "combined_sideband_fit": {k: v for k, v in combined.items() if k != "background"},
        "holdout_validation": {
            "low_fit_predicts_high_deviance": region_deviance(data_a, low_expected, HIGH),
            "high_fit_predicts_low_deviance": region_deviance(data_a, high_expected, LOW),
            "low_fit": {k: v for k, v in low_fit.items() if k != "background"},
            "high_fit": {k: v for k, v in high_fit.items() if k != "background"},
        },
        "signal_window": {
            "bounds_gev": [80.0, 100.0],
            "data": float(data_a[SIGNAL].sum()),
            "prompt": float(prompt[SIGNAL].sum()),
            "functional_background": combined["background_signal_yield"],
            "raw_abcd_background": abcd,
        },
        "caveats": [
            "Without an external prompt scale, normalizing DY in the signal window is circular.",
            "The model describes a smooth residual only and is not a substitute for missing prompt-background samples.",
            "Low/high holdout disagreement must be treated as a model-failure diagnostic, not tuned away.",
        ],
        "plot": "sideband_exponential_interpolation.svg",
    }
    (output / "background_fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["signal_window"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
