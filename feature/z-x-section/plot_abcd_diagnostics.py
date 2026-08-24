#!/usr/bin/env python3
"""Prompt-subtracted ABCD closure calculation and SVG diagnostics."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot


HERE = Path(__file__).resolve().parent

MASS_CONFIGURATION = {
    "plot_range_gev": [60.0, 120.0],
    "histogram_bin_width_gev": 1.0,
    "ratio_edges_gev": [60.0, 65.0, 70.0, 75.0, 80.0, 90.0,
                        100.0, 105.0, 110.0, 115.0, 120.0],
    "regions": {
        "low_sideband": [60.0, 75.0],
        "signal": [80.0, 100.0],
        "high_sideband": [105.0, 120.0],
    },
    "vertical_boundaries_gev": [75.0, 80.0, 100.0, 105.0],
}

ISOLATION_CONFIGURATION = {
    "isolated_max": 0.15,
    "anti_isolated_min": 0.25,
    "modes": ("both", "at_least_one"),
}

REGION_LABELS = {
    "A": "A: OS, isolated",
    "B": "B: SS, isolated",
    "C": "C: OS, anti-isolated",
    "D": "D: SS, anti-isolated",
}

COLORS = {"A": "black", "B": "#3b82f6", "C": "#ef4444", "D": "#10b981"}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "selection_dir", type=Path,
        help="Output directory from a 60-120 GeV run_z_selection.py skim production",
    )
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--label", default="local_sideband_closure")
    parser.add_argument("--isolated-max", type=float,
                        default=ISOLATION_CONFIGURATION["isolated_max"])
    parser.add_argument("--anti-isolated-min", type=float,
                        default=ISOLATION_CONFIGURATION["anti_isolated_min"])
    parser.add_argument("--maximum-iterations", type=int, default=20)
    parser.add_argument(
        "--prompt-scale", type=float,
        help=("Fix the prompt-MC normalization instead of deriving it from the "
              "OS-isolated peak; use this for a luminosity-normalized production"),
    )
    parser.add_argument(
        "--prompt-scale-uncertainty", type=float, default=0.0,
        help="Absolute uncertainty on --prompt-scale",
    )
    return parser.parse_args()


def load_skim(filename: Path) -> dict[str, np.ndarray]:
    if not filename.is_file():
        raise FileNotFoundError(f"Missing diagnostic skim: {filename}")
    branches = [
        "event_weight", "lead_pt", "sublead_pt", "lead_eta", "sublead_eta",
        "lead_charge", "sublead_charge", "lead_rel_iso", "sublead_rel_iso",
        "dimuon_mass",
    ]
    with uproot.open(f"{filename}:Events") as tree:
        arrays = tree.arrays(branches, library="ak")
    return {name: ak.to_numpy(arrays[name]) for name in branches}


def mass_mask(values: np.ndarray, bounds: list[float]) -> np.ndarray:
    return (values >= bounds[0]) & (values < bounds[1])


def categories(sample: dict[str, np.ndarray], mode: str,
               isolated_max: float, anti_isolated_min: float) -> dict[str, np.ndarray]:
    opposite = sample["lead_charge"] * sample["sublead_charge"] < 0
    isolated = ((sample["lead_rel_iso"] < isolated_max)
                & (sample["sublead_rel_iso"] < isolated_max))
    if mode == "both":
        anti = ((sample["lead_rel_iso"] > anti_isolated_min)
                & (sample["sublead_rel_iso"] > anti_isolated_min))
    elif mode == "at_least_one":
        anti = ((sample["lead_rel_iso"] > anti_isolated_min)
                | (sample["sublead_rel_iso"] > anti_isolated_min))
    else:
        raise ValueError(f"Unknown anti-isolation mode {mode}")
    return {
        "A": opposite & isolated,
        "B": ~opposite & isolated,
        "C": opposite & anti,
        "D": ~opposite & anti,
    }


def weighted_yield(sample: dict[str, np.ndarray], mask: np.ndarray) -> tuple[float, float]:
    weights = sample["event_weight"][mask].astype(float)
    return float(weights.sum()), float(np.square(weights).sum())


def component(data: dict[str, np.ndarray], prompt: dict[str, np.ndarray],
              data_mask: np.ndarray, prompt_mask: np.ndarray,
              alpha: float, alpha_variance: float) -> dict[str, float]:
    data_yield, data_variance = weighted_yield(data, data_mask)
    mc_yield, mc_variance = weighted_yield(prompt, prompt_mask)
    prompt_yield = alpha * mc_yield
    prompt_variance = alpha * alpha * mc_variance + mc_yield * mc_yield * alpha_variance
    return {
        "data": data_yield,
        "data_variance": data_variance,
        "prompt": prompt_yield,
        "prompt_variance": prompt_variance,
        "nonprompt": data_yield - prompt_yield,
        "nonprompt_variance": data_variance + prompt_variance,
        "raw_prompt_mc": mc_yield,
        "raw_prompt_mc_variance": mc_variance,
    }


def base_prediction(parts: dict[str, dict[str, float]]) -> tuple[float | None, float | None]:
    b, c, d = (parts[key]["nonprompt"] for key in ("B", "C", "D"))
    vb, vc, vd = (parts[key]["nonprompt_variance"] for key in ("B", "C", "D"))
    if b <= 0.0 or c <= 0.0 or d <= 0.0:
        return None, None
    prediction = b * c / d
    variance = ((c / d) ** 2 * vb + (b / d) ** 2 * vc
                + (b * c / (d * d)) ** 2 * vd)
    return prediction, variance


def window_components(data: dict[str, np.ndarray], prompt: dict[str, np.ndarray],
                      data_categories: dict[str, np.ndarray],
                      prompt_categories: dict[str, np.ndarray], bounds: list[float],
                      alpha: float, alpha_variance: float) -> dict[str, dict[str, float]]:
    dm = mass_mask(data["dimuon_mass"], bounds)
    pm = mass_mask(prompt["dimuon_mass"], bounds)
    return {
        region: component(data, prompt, dm & data_categories[region],
                          pm & prompt_categories[region], alpha, alpha_variance)
        for region in REGION_LABELS
    }


def combined_sideband_mask(sample: dict[str, np.ndarray]) -> np.ndarray:
    regions = MASS_CONFIGURATION["regions"]
    return (mass_mask(sample["dimuon_mass"], regions["low_sideband"])
            | mass_mask(sample["dimuon_mass"], regions["high_sideband"]))


def fit_closure_factor(data: dict[str, np.ndarray], prompt: dict[str, np.ndarray],
                       data_categories: dict[str, np.ndarray],
                       prompt_categories: dict[str, np.ndarray], alpha: float,
                       alpha_variance: float) -> dict[str, float | None]:
    dm, pm = combined_sideband_mask(data), combined_sideband_mask(prompt)
    parts = {
        region: component(data, prompt, dm & data_categories[region],
                          pm & prompt_categories[region], alpha, alpha_variance)
        for region in REGION_LABELS
    }
    base, base_variance = base_prediction(parts)
    observed = parts["A"]["nonprompt"]
    observed_variance = parts["A"]["nonprompt_variance"]
    if base is None or base <= 0.0 or observed <= 0.0:
        return {
            "value": None, "variance": None, "stat_uncertainty": None,
            "observed_sideband_nonprompt": observed,
            "base_sideband_prediction": base,
        }
    value = observed / base
    variance = value * value * (
        observed_variance / (observed * observed)
        + (base_variance or 0.0) / (base * base)
    )
    return {
        "value": value,
        "variance": variance,
        "stat_uncertainty": math.sqrt(max(0.0, variance)),
        "observed_sideband_nonprompt": observed,
        "base_sideband_prediction": base,
    }


def solve_normalization_and_closure(
    data: dict[str, np.ndarray], prompt: dict[str, np.ndarray], mode: str,
    isolated_max: float, anti_isolated_min: float, maximum_iterations: int,
    fixed_prompt_scale: float | None = None,
    fixed_prompt_scale_uncertainty: float = 0.0,
) -> dict[str, Any]:
    dc = categories(data, mode, isolated_max, anti_isolated_min)
    pc = categories(prompt, mode, isolated_max, anti_isolated_min)
    signal = MASS_CONFIGURATION["regions"]["signal"]
    dm, pm = mass_mask(data["dimuon_mass"], signal), mass_mask(prompt["dimuon_mass"], signal)
    data_a, data_a_variance = weighted_yield(data, dm & dc["A"])
    mc_a, mc_a_variance = weighted_yield(prompt, pm & pc["A"])
    if mc_a <= 0.0:
        raise RuntimeError("The prompt MC signal-region normalization yield is non-positive")
    closure: dict[str, Any] = {}
    if fixed_prompt_scale is not None:
        if fixed_prompt_scale < 0.0 or fixed_prompt_scale_uncertainty < 0.0:
            raise ValueError("Fixed prompt scale and its uncertainty must be non-negative")
        alpha = fixed_prompt_scale
        alpha_variance = fixed_prompt_scale_uncertainty ** 2
        converged = True
        iterations = 0
        normalization_method = "fixed prompt-MC scale supplied by user"
    else:
        alpha = data_a / mc_a
        alpha_variance = data_a_variance / (mc_a * mc_a)
        alpha_variance += data_a * data_a * mc_a_variance / (mc_a ** 4)
        converged = False
        iterations = maximum_iterations
        normalization_method = "iterative OS-isolated peak normalization after ABCD subtraction"
        for iteration in range(maximum_iterations):
            closure = fit_closure_factor(data, prompt, dc, pc, alpha, alpha_variance)
            kappa = closure["value"] if closure["value"] is not None else 1.0
            signal_parts = window_components(
                data, prompt, dc, pc, signal, alpha, alpha_variance
            )
            base, base_variance = base_prediction(signal_parts)
            predicted_nonprompt = kappa * base if base is not None else 0.0
            prediction_variance = 0.0
            if base is not None and base_variance is not None:
                prediction_variance = kappa * kappa * base_variance
                if closure["variance"] is not None:
                    prediction_variance += base * base * closure["variance"]
            new_alpha = (data_a - predicted_nonprompt) / mc_a
            new_alpha_variance = (data_a_variance + prediction_variance) / (mc_a * mc_a)
            new_alpha_variance += ((data_a - predicted_nonprompt) ** 2
                                   * mc_a_variance / (mc_a ** 4))
            if abs(new_alpha - alpha) < 1e-10 * max(1.0, abs(alpha)):
                alpha, alpha_variance, converged = new_alpha, new_alpha_variance, True
                iterations = iteration + 1
                break
            alpha, alpha_variance = new_alpha, new_alpha_variance

    closure = fit_closure_factor(data, prompt, dc, pc, alpha, alpha_variance)
    windows: dict[str, Any] = {}
    for name, bounds in MASS_CONFIGURATION["regions"].items():
        parts = window_components(data, prompt, dc, pc, bounds, alpha, alpha_variance)
        base, base_variance = base_prediction(parts)
        kappa = closure["value"]
        prediction = kappa * base if kappa is not None and base is not None else None
        prediction_variance = None
        if prediction is not None and base_variance is not None:
            prediction_variance = kappa * kappa * base_variance
            if closure["variance"] is not None:
                prediction_variance += base * base * closure["variance"]
        windows[name] = {
            "bounds_gev": bounds,
            "components": parts,
            "base_bc_over_d": base,
            "base_variance": base_variance,
            "kappa_corrected_prediction": prediction,
            "prediction_variance": prediction_variance,
            "observed_A_nonprompt": parts["A"]["nonprompt"],
            "observed_A_variance": parts["A"]["nonprompt_variance"],
        }
    return {
        "mode": mode,
        "prompt_scale": alpha,
        "prompt_scale_variance": alpha_variance,
        "prompt_scale_stat_uncertainty": math.sqrt(max(0.0, alpha_variance)),
        "normalization_method": normalization_method,
        "converged": converged,
        "iterations": iterations,
        "closure_factor": closure,
        "windows": windows,
    }


def histogram(values: np.ndarray, mask: np.ndarray, edges: np.ndarray,
              weights: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    selected_weights = None if weights is None else weights[mask]
    content, _ = np.histogram(values[mask], bins=edges, weights=selected_weights)
    variance_weights = None if selected_weights is None else np.square(selected_weights)
    variance, _ = np.histogram(values[mask], bins=edges, weights=variance_weights)
    if selected_weights is None:
        variance = content.astype(float)
    return content.astype(float), variance.astype(float)


def cms_label(ax: Any) -> None:
    hep.cms.label("Open Data", data=True, com=13, ax=ax)


def draw_mass_regions(data: dict[str, np.ndarray], prompt: dict[str, np.ndarray],
                      result: dict[str, Any], isolated_max: float,
                      anti_isolated_min: float, output: Path) -> None:
    mode = result["mode"]
    dc = categories(data, mode, isolated_max, anti_isolated_min)
    pc = categories(prompt, mode, isolated_max, anti_isolated_min)
    low, high = MASS_CONFIGURATION["plot_range_gev"]
    width = MASS_CONFIGURATION["histogram_bin_width_gev"]
    edges = np.arange(low, high + 0.5 * width, width)
    fig, ax = plt.subplots(figsize=(10, 7))
    for region in REGION_LABELS:
        values, _ = histogram(data["dimuon_mass"], dc[region], edges)
        hep.histplot(values, edges, ax=ax, histtype="step", linewidth=1.8,
                     color=COLORS[region], label=REGION_LABELS[region])
        prompt_values, _ = histogram(
            prompt["dimuon_mass"], pc[region], edges,
            prompt["event_weight"] * result["prompt_scale"],
        )
        hep.histplot(prompt_values, edges, ax=ax, histtype="step", linewidth=1.0,
                     linestyle="--", color=COLORS[region], alpha=0.65)
    for boundary in MASS_CONFIGURATION["vertical_boundaries_gev"]:
        ax.axvline(boundary, color="0.45", linestyle="--", linewidth=1.0)
    ax.set_yscale("log")
    ax.set_xlim(low, high)
    ax.set_ylim(bottom=0.5)
    ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    ax.set_ylabel("Events / GeV")
    ax.legend(ncol=2, fontsize=11)
    ax.text(0.03, 0.75, "Solid: data\nDashed: scaled prompt DY",
            transform=ax.transAxes, fontsize=10)
    cms_label(ax)
    fig.tight_layout()
    fig.savefig(output / f"abcd_mass_regions_{mode}.svg")
    plt.close(fig)


def draw_isolation(data: dict[str, np.ndarray], prompt: dict[str, np.ndarray],
                   prompt_scale: float, isolated_max: float,
                   anti_isolated_min: float, output: Path) -> None:
    opposite_data = data["lead_charge"] * data["sublead_charge"] < 0
    opposite_prompt = prompt["lead_charge"] * prompt["sublead_charge"] < 0
    edges = np.linspace(0.0, 1.0, 51)
    for prefix, branch, label in (
        ("leading", "lead_rel_iso", "leading muon"),
        ("subleading", "sublead_rel_iso", "subleading muon"),
    ):
        fig, ax = plt.subplots(figsize=(9, 7))
        clipped_data = np.clip(data[branch], edges[0], np.nextafter(edges[-1], edges[0]))
        clipped_prompt = np.clip(prompt[branch], edges[0], np.nextafter(edges[-1], edges[0]))
        for mask, color, name in (
            (opposite_data, "black", "OS data"),
            (~opposite_data, "#3b82f6", "SS data"),
        ):
            content, _ = histogram(clipped_data, mask, edges)
            hep.histplot(content, edges, ax=ax, histtype="step", linewidth=1.8,
                         color=color, label=name)
        prompt_content, _ = histogram(
            clipped_prompt, opposite_prompt, edges,
            prompt["event_weight"] * prompt_scale,
        )
        hep.histplot(prompt_content, edges, ax=ax, histtype="step", linewidth=1.4,
                     linestyle="--", color="#ef4444", label="OS scaled prompt DY")
        ax.axvline(isolated_max, color="0.4", linestyle="--")
        ax.axvline(anti_isolated_min, color="0.4", linestyle="--")
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.5)
        ax.set_xlabel(f"{label} relative isolation (overflow in last bin)")
        ax.set_ylabel("Events")
        ax.legend(fontsize=11)
        cms_label(ax)
        fig.tight_layout()
        fig.savefig(output / f"{prefix}_relative_isolation.svg")
        plt.close(fig)


def ratio_and_variance(numerator: np.ndarray, numerator_variance: np.ndarray,
                       denominator: np.ndarray,
                       denominator_variance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ratio = np.full_like(numerator, np.nan, dtype=float)
    variance = np.full_like(numerator, np.nan, dtype=float)
    valid = (numerator > 0.0) & (denominator > 0.0)
    ratio[valid] = numerator[valid] / denominator[valid]
    variance[valid] = np.square(ratio[valid]) * (
        numerator_variance[valid] / np.square(numerator[valid])
        + denominator_variance[valid] / np.square(denominator[valid])
    )
    return ratio, variance


def draw_os_ss_ratio(data: dict[str, np.ndarray], prompt: dict[str, np.ndarray],
                     results: dict[str, dict[str, Any]], isolated_max: float,
                     anti_isolated_min: float, output: Path) -> None:
    edges = np.asarray(MASS_CONFIGURATION["ratio_edges_gev"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    for mode, color, marker in (("both", "#111827", "o"),
                                ("at_least_one", "#ef4444", "s")):
        dc = categories(data, mode, isolated_max, anti_isolated_min)
        pc = categories(prompt, mode, isolated_max, anti_isolated_min)
        os_raw, os_raw_var = histogram(data["dimuon_mass"], dc["C"], edges)
        ss_raw, ss_raw_var = histogram(data["dimuon_mass"], dc["D"], edges)
        raw_ratio, raw_variance = ratio_and_variance(os_raw, os_raw_var, ss_raw, ss_raw_var)
        axes[0].errorbar(centers, raw_ratio, yerr=np.sqrt(raw_variance), fmt=marker,
                         color=color, label=mode.replace("_", " "))

        alpha = results[mode]["prompt_scale"]
        alpha_var = results[mode]["prompt_scale_variance"]
        os_mc, os_mc_var = histogram(
            prompt["dimuon_mass"], pc["C"], edges, prompt["event_weight"]
        )
        ss_mc, ss_mc_var = histogram(
            prompt["dimuon_mass"], pc["D"], edges, prompt["event_weight"]
        )
        os_sub = os_raw - alpha * os_mc
        ss_sub = ss_raw - alpha * ss_mc
        os_sub_var = os_raw_var + alpha * alpha * os_mc_var + np.square(os_mc) * alpha_var
        ss_sub_var = ss_raw_var + alpha * alpha * ss_mc_var + np.square(ss_mc) * alpha_var
        sub_ratio, sub_variance = ratio_and_variance(os_sub, os_sub_var, ss_sub, ss_sub_var)
        axes[1].errorbar(centers, sub_ratio, yerr=np.sqrt(sub_variance), fmt=marker,
                         color=color, label=mode.replace("_", " "))
    for ax, title in zip(axes, ("Raw data", "After scaled prompt-DY subtraction")):
        for boundary in MASS_CONFIGURATION["vertical_boundaries_gev"]:
            ax.axvline(boundary, color="0.6", linestyle="--", linewidth=0.9)
        ax.set_ylabel("OS / SS")
        ax.set_title(title, loc="left", fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(axis="y", alpha=0.2)
    axes[-1].set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    cms_label(axes[0])
    fig.tight_layout()
    fig.savefig(output / "anti_isolated_os_ss_ratio_vs_mass.svg")
    plt.close(fig)


def draw_closure(results: dict[str, dict[str, Any]], output: Path) -> None:
    for mode, result in results.items():
        names = list(MASS_CONFIGURATION["regions"])
        centers = np.asarray([
            sum(MASS_CONFIGURATION["regions"][name]) / 2.0 for name in names
        ])
        xerr = np.asarray([
            (MASS_CONFIGURATION["regions"][name][1]
             - MASS_CONFIGURATION["regions"][name][0]) / 2.0 for name in names
        ])
        observed = np.asarray([result["windows"][name]["observed_A_nonprompt"] for name in names])
        observed_error = np.sqrt(np.asarray([
            result["windows"][name]["observed_A_variance"] for name in names
        ]))
        predicted = np.asarray([
            np.nan if result["windows"][name]["kappa_corrected_prediction"] is None
            else result["windows"][name]["kappa_corrected_prediction"] for name in names
        ])
        predicted_error = np.asarray([
            np.nan if result["windows"][name]["prediction_variance"] is None
            else math.sqrt(max(0.0, result["windows"][name]["prediction_variance"]))
            for name in names
        ])
        ratio, ratio_var = ratio_and_variance(
            observed, np.square(observed_error), predicted, np.square(predicted_error)
        )
        fig, (upper, lower) = plt.subplots(
            2, 1, figsize=(10, 8), sharex=True,
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06},
        )
        upper.errorbar(centers - 0.4, observed, xerr=xerr, yerr=observed_error,
                       fmt="o", color="black", label=r"Observed $A_{data}-A_{prompt}$")
        upper.errorbar(centers + 0.4, predicted, xerr=xerr, yerr=predicted_error,
                       fmt="s", color="#ef4444", label=r"$\kappa(B-B_p)(C-C_p)/(D-D_p)$")
        kappa = result["closure_factor"]["value"]
        kappa_unc = result["closure_factor"]["stat_uncertainty"]
        text = "Sideband closure undefined" if kappa is None else f"sideband kappa = {kappa:.3g} +/- {kappa_unc:.2g}"
        upper.text(0.03, 0.76, text, transform=upper.transAxes, fontsize=11)
        upper.set_ylabel("Nonprompt events")
        upper.legend(fontsize=10)
        upper.grid(axis="y", alpha=0.2)
        lower.errorbar(centers, ratio, xerr=xerr, yerr=np.sqrt(ratio_var),
                       fmt="o", color="#374151")
        lower.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
        lower.set_ylabel("Obs./pred.")
        lower.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
        lower.set_ylim(0.0, 2.5)
        lower.grid(axis="y", alpha=0.2)
        cms_label(upper)
        fig.savefig(output / f"bc_over_d_closure_vs_mass_{mode}.svg", bbox_inches="tight")
        plt.close(fig)


def main() -> int:
    args = arguments()
    if not 0.0 < args.isolated_max < args.anti_isolated_min:
        raise SystemExit("Require 0 < isolated maximum < anti-isolated minimum")
    data = load_skim(args.selection_dir / "data_abcd_skim.root")
    prompt = load_skim(args.selection_dir / "mc_abcd_skim.root")
    output = args.output_dir / args.label
    output.mkdir(parents=True, exist_ok=True)
    hep.style.use("CMS")

    results = {
        mode: solve_normalization_and_closure(
            data, prompt, mode, args.isolated_max, args.anti_isolated_min,
            args.maximum_iterations, args.prompt_scale,
            args.prompt_scale_uncertainty,
        )
        for mode in ISOLATION_CONFIGURATION["modes"]
    }
    for mode, result in results.items():
        draw_mass_regions(data, prompt, result, args.isolated_max,
                          args.anti_isolated_min, output)
    draw_isolation(data, prompt, results["both"]["prompt_scale"],
                   args.isolated_max, args.anti_isolated_min, output)
    draw_os_ss_ratio(data, prompt, results, args.isolated_max,
                     args.anti_isolated_min, output)
    draw_closure(results, output)

    summary = {
        "title": "Prompt-subtracted charge/isolation ABCD sideband closure",
        "selection_directory": str(args.selection_dir.resolve()),
        "configuration": {
            "mass": MASS_CONFIGURATION,
            "isolation": {
                "isolated_max": args.isolated_max,
                "anti_isolated_min": args.anti_isolated_min,
                "nominal_mode": "both",
                "systematic_mode": "at_least_one",
            },
            "equation": "A_nonprompt = kappa*(B_data-B_prompt)*(C_data-C_prompt)/(D_data-D_prompt)",
            "fixed_prompt_scale": args.prompt_scale,
            "fixed_prompt_scale_uncertainty": args.prompt_scale_uncertainty,
        },
        "results": results,
        "caveats": [
            "Prompt DY is shape-normalized to the OS-isolated peak because the local subset has no luminosity normalization.",
            "The prompt normalization and closure factor are solved iteratively and are provisional.",
            "No top, diboson, W+jets, Z->tautau, or dedicated QCD simulation is included.",
            "The closure-factor uncertainty uses approximate independent-yield propagation and neglects category covariances.",
            "The local same-sign isolated sample is very sparse; the full Run2016H production is required for a decision.",
        ],
        "plots": sorted(path.name for path in output.glob("*.svg")),
    }
    (output / "closure_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {len(summary['plots'])} SVG plots and {output / 'closure_summary.json'}")
    for mode, result in results.items():
        closure = result["closure_factor"]
        print(
            f"{mode}: prompt scale={result['prompt_scale']:.6g}, "
            f"kappa={closure['value']}, uncertainty={closure['stat_uncertainty']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
