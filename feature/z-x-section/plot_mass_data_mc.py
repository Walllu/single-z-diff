#!/usr/bin/env python3
"""Plot the selected dimuon mass in data and absolutely normalized DY MC.

The distributed production already contains wide-window ABCD skims.  This
utility uses their OS, isolated events, applies the golden JSON to data, and
reconstructs the nominal muon-efficiency weight for MC.  The skim predates the
measurement-level momentum correction, so the resulting mass shape is a
validation diagnostic rather than an exact reproduction of the final yield.

A global background fraction fixes an integral but not a mass shape.  For the
requested background-added plots, the background is therefore distributed
proportionally to the DY template and is prominently labelled as a shape proxy.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot


HERE = Path(__file__).resolve().parent
DEFAULT_SF_DIR = HERE.parent / "muon-efficiency-sfs"
DEFAULT_GOLDEN_JSON = HERE.parent.parent / "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("selection_directory", type=Path,
                        help="Directory containing data_abcd_skim.root and mc_abcd_skim.root")
    result.add_argument("--cross-section-result", type=Path,
                        help="Final result JSON; supplies MC normalization and background fraction")
    result.add_argument("--mc-normalization-factor", type=float,
                        help="Override L*xsec/sum(genWeight) normalization")
    result.add_argument("--background-fraction", type=float,
                        help="Override the fraction of observed candidates assigned to background")
    result.add_argument("--golden-json", type=Path, default=DEFAULT_GOLDEN_JSON)
    result.add_argument("--sf-directory", type=Path, default=DEFAULT_SF_DIR)
    result.add_argument("--output-dir", type=Path, default=HERE / "plots")
    result.add_argument("--label", default="mass_data_mc")
    result.add_argument("--mass-range", type=float, nargs=2, default=(60.0, 120.0),
                        metavar=("LOW", "HIGH"))
    result.add_argument("--analysis-window", type=float, nargs=2,
                        default=(76.1876, 106.1876), metavar=("LOW", "HIGH"))
    result.add_argument("--bin-width", type=float, default=1.0)
    result.add_argument("--isolated-max", type=float, default=0.15)
    result.add_argument("--step-size", default="200 MB")
    return result


def read_inputs(args: argparse.Namespace) -> tuple[float, float, dict[str, Any] | None]:
    payload = None
    if args.cross_section_result:
        payload = json.loads(args.cross_section_result.read_text())
    factor = args.mc_normalization_factor
    if factor is None and payload is not None:
        diagnostic = payload.get("signal_mc_normalization_diagnostic") or {}
        factor = diagnostic.get("normalization_factor")
    if factor is None:
        raise ValueError("Supply --cross-section-result or --mc-normalization-factor")
    fraction = args.background_fraction
    if fraction is None and payload is not None:
        manual = (payload.get("prompt_background") or {}).get("manual_component") or {}
        fraction = manual.get("fraction_of_observed_candidates")
    if fraction is None:
        fraction = 0.0
    if factor <= 0.0 or not 0.0 <= fraction < 1.0:
        raise ValueError("MC normalization must be positive and background fraction in [0,1)")
    return float(factor), float(fraction), payload


def load_golden(filename: Path) -> dict[int, list[tuple[int, int]]]:
    raw = json.loads(filename.read_text())
    return {int(run): [(int(a), int(b)) for a, b in ranges]
            for run, ranges in raw.items()}


def certified_mask(runs: np.ndarray, lumis: np.ndarray,
                   golden: dict[int, list[tuple[int, int]]]) -> np.ndarray:
    result = np.zeros(len(runs), dtype=bool)
    for run in np.unique(runs):
        ranges = golden.get(int(run))
        if not ranges:
            continue
        same_run = runs == run
        selected = np.zeros(np.count_nonzero(same_run), dtype=bool)
        run_lumis = lumis[same_run]
        for low, high in ranges:
            selected |= (run_lumis >= low) & (run_lumis <= high)
        result[same_run] = selected
    return result


def load_sf_maps(directory: Path) -> dict[str, np.ndarray]:
    filenames = {
        "reco": "MuonReco_ScaleFactors_Nominal.txt",
        "iso": "MuonIso_ScaleFactors_Nominal.txt",
        "trigger": "MuonTrigger_ScaleFactors_Nominal.txt",
    }
    maps = {name: np.loadtxt(directory / filename, comments="#")
            for name, filename in filenames.items()}
    for name, values in maps.items():
        if values.shape[0] != 10:
            raise ValueError(f"Unexpected {name} SF shape {values.shape}")
        if values.shape[1] == 10:
            maps[name] = np.pad(values, ((0, 0), (0, 5)), mode="edge")
        elif values.shape[1] != 15:
            raise ValueError(f"Unexpected {name} SF shape {values.shape}")
    return maps


def lookup(values: np.ndarray, pt: np.ndarray, eta: np.ndarray) -> np.ndarray:
    eta_index = np.clip(np.searchsorted(np.linspace(-2.5, 2.5, 11), eta, side="right") - 1,
                        0, 9)
    pt_index = np.clip(np.searchsorted(np.linspace(25.0, 100.0, 16), pt, side="right") - 1,
                       0, 14)
    return values[eta_index, pt_index]


def selected(arrays: dict[str, np.ndarray], isolated_max: float) -> np.ndarray:
    return ((arrays["lead_charge"] * arrays["sublead_charge"] < 0)
            & (arrays["lead_rel_iso"] < isolated_max)
            & (arrays["sublead_rel_iso"] < isolated_max))


def accumulate_data(filename: Path, edges: np.ndarray, golden: dict[int, list[tuple[int, int]]],
                    isolated_max: float, analysis_window: tuple[float, float],
                    step_size: str) -> tuple[np.ndarray, np.ndarray, int, float]:
    content = np.zeros(len(edges) - 1)
    variance = np.zeros_like(content)
    accepted = 0
    window_yield = 0.0
    branches = ["run", "lumi", "lead_charge", "sublead_charge", "lead_rel_iso",
                "sublead_rel_iso", "dimuon_mass"]
    for chunk in uproot.iterate(f"{filename}:Events", branches, step_size=step_size, library="ak"):
        arrays = {name: ak.to_numpy(chunk[name]) for name in branches}
        mask = selected(arrays, isolated_max) & certified_mask(arrays["run"], arrays["lumi"], golden)
        masses = arrays["dimuon_mass"][mask]
        values, _ = np.histogram(masses, bins=edges)
        content += values
        variance += values
        accepted += int(len(masses))
        window_yield += float(np.count_nonzero(
            (masses > analysis_window[0]) & (masses < analysis_window[1])
        ))
    return content, variance, accepted, window_yield


def accumulate_mc(filename: Path, edges: np.ndarray, maps: dict[str, np.ndarray],
                  normalization: float, isolated_max: float,
                  analysis_window: tuple[float, float],
                  step_size: str) -> tuple[np.ndarray, np.ndarray, int, float]:
    content = np.zeros(len(edges) - 1)
    variance = np.zeros_like(content)
    accepted = 0
    window_yield = 0.0
    branches = ["event_weight", "lead_pt", "sublead_pt", "lead_eta", "sublead_eta",
                "lead_charge", "sublead_charge", "lead_trigger_matched",
                "sublead_trigger_matched", "lead_rel_iso", "sublead_rel_iso", "dimuon_mass"]
    for chunk in uproot.iterate(f"{filename}:Events", branches, step_size=step_size, library="ak"):
        arrays = {name: ak.to_numpy(chunk[name]) for name in branches}
        mask = selected(arrays, isolated_max)
        lead_reco = lookup(maps["reco"], arrays["lead_pt"], arrays["lead_eta"])
        sublead_reco = lookup(maps["reco"], arrays["sublead_pt"], arrays["sublead_eta"])
        lead_iso = lookup(maps["iso"], arrays["lead_pt"], arrays["lead_eta"])
        sublead_iso = lookup(maps["iso"], arrays["sublead_pt"], arrays["sublead_eta"])
        trigger_pt = np.where(arrays["lead_trigger_matched"], arrays["lead_pt"], arrays["sublead_pt"])
        trigger_eta = np.where(arrays["lead_trigger_matched"], arrays["lead_eta"], arrays["sublead_eta"])
        trigger = lookup(maps["trigger"], trigger_pt, trigger_eta)
        weights = (normalization * arrays["event_weight"] * lead_reco * sublead_reco
                   * lead_iso * sublead_iso * trigger)
        masses = arrays["dimuon_mass"][mask]
        weights = weights[mask]
        values, _ = np.histogram(masses, bins=edges, weights=weights)
        sumw2, _ = np.histogram(masses, bins=edges, weights=np.square(weights))
        content += values
        variance += sumw2
        accepted += int(len(masses))
        in_window = (masses > analysis_window[0]) & (masses < analysis_window[1])
        window_yield += float(weights[in_window].sum())
    return content, variance, accepted, window_yield


def ratio(data: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    return np.divide(data, prediction, out=np.full_like(data, np.nan), where=prediction > 0.0)


def draw(edges: np.ndarray, data: np.ndarray, data_variance: np.ndarray,
         dy: np.ndarray, dy_variance: np.ndarray, background: np.ndarray,
         analysis_window: tuple[float, float], logarithmic: bool,
         include_background: bool, output: Path) -> None:
    centers = (edges[:-1] + edges[1:]) / 2.0
    widths = np.diff(edges)
    prediction = dy + background if include_background else dy
    prediction_variance = dy_variance
    fig, (upper, lower) = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                                       gridspec_kw={"height_ratios": [3.2, 1], "hspace": 0.05})
    upper.errorbar(centers, data, xerr=widths / 2.0, yerr=np.sqrt(data_variance), fmt="o",
                   markersize=3.5, color="black", label="Data")
    hep.histplot(dy, edges, ax=upper, histtype="step", linewidth=1.8, color="#2563eb",
                 label="DY MC (absolute normalization)")
    if include_background:
        hep.histplot(background, edges, ax=upper, histtype="fill", alpha=0.28,
                     color="#f59e0b", label="Global-background shape proxy")
        hep.histplot(prediction, edges, ax=upper, histtype="step", linewidth=1.8,
                     color="#dc2626", label="DY + background proxy")
    for boundary in analysis_window:
        upper.axvline(boundary, color="0.45", linestyle="--", linewidth=1.0)
        lower.axvline(boundary, color="0.45", linestyle="--", linewidth=1.0)
    if logarithmic:
        upper.set_yscale("log")
        positive = np.concatenate((data[data > 0.0], prediction[prediction > 0.0]))
        if len(positive):
            upper.set_ylim(max(0.5, positive.min() * 0.4), positive.max() * 5.0)
    else:
        upper.set_ylim(0.0, max(data.max(initial=0.0), prediction.max(initial=0.0)) * 1.28)
    upper.set_ylabel(f"Events / {widths[0]:g} GeV")
    upper.legend(fontsize=11)
    hep.cms.label("Open Data", data=True, com=13, ax=upper)
    lower.errorbar(centers, ratio(data, prediction), xerr=widths / 2.0,
                   yerr=np.divide(np.sqrt(data_variance), prediction,
                                  out=np.zeros_like(data), where=prediction > 0.0),
                   fmt="o", markersize=3.2, color="black")
    relative_mc = np.divide(np.sqrt(prediction_variance), prediction,
                            out=np.zeros_like(prediction), where=prediction > 0.0)
    lower.fill_between(edges, np.r_[1.0-relative_mc, 1.0-relative_mc[-1]],
                       np.r_[1.0+relative_mc, 1.0+relative_mc[-1]], step="post",
                       color="#93c5fd", alpha=0.45, linewidth=0, label="MC stat.")
    lower.axhline(1.0, color="0.25", linewidth=1.0)
    finite = ratio(data, prediction)
    finite = finite[np.isfinite(finite)]
    lower_min = 0.5 if len(finite) == 0 else max(0.0, min(0.5, finite.min() * 0.88))
    lower_max = 1.5 if len(finite) == 0 else max(1.5, min(3.0, finite.max() * 1.12))
    lower.set_ylim(lower_min, lower_max)
    lower.set_ylabel("Data/pred.")
    lower.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    lower.grid(axis="y", alpha=0.2)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parser().parse_args()
    normalization, fraction, _source = read_inputs(args)
    if args.bin_width <= 0.0 or args.mass_range[1] <= args.mass_range[0]:
        raise ValueError("Invalid mass range or bin width")
    bins = round((args.mass_range[1] - args.mass_range[0]) / args.bin_width)
    edges = np.linspace(args.mass_range[0], args.mass_range[1], bins + 1)
    output = args.output_dir / args.label
    output.mkdir(parents=True, exist_ok=True)
    data_file = args.selection_directory / "data_abcd_skim.root"
    mc_file = args.selection_directory / "mc_abcd_skim.root"
    golden = load_golden(args.golden_json)
    maps = load_sf_maps(args.sf_directory)
    window = tuple(float(x) for x in args.analysis_window)
    data, data_variance, data_entries, data_window = accumulate_data(
        data_file, edges, golden, args.isolated_max, window, args.step_size
    )
    dy, dy_variance, mc_entries, dy_window = accumulate_mc(
        mc_file, edges, maps, normalization, args.isolated_max, window, args.step_size
    )
    background_target = fraction * data_window
    if dy_window <= 0.0:
        raise RuntimeError("DY integral in the analysis window is non-positive")
    background = dy * (background_target / dy_window)

    plots = []
    for include_background in (False, True):
        component = "with_global_background_proxy" if include_background else "dy_only"
        for logarithmic in (False, True):
            scale = "log" if logarithmic else "linear"
            filename = f"dimuon_mass_{component}_{scale}.svg"
            draw(edges, data, data_variance, dy, dy_variance, background, window,
                 logarithmic, include_background, output / filename)
            plots.append(filename)

    prediction_with = dy + background
    summary = {
        "schema_version": 1,
        "title": "Absolute DY-to-data dimuon-mass comparison",
        "inputs": {
            "selection_directory": str(args.selection_directory.resolve()),
            "cross_section_result": (None if args.cross_section_result is None
                                     else str(args.cross_section_result.resolve())),
            "golden_json": str(args.golden_json.resolve()),
            "muon_efficiency_sf_directory": str(args.sf_directory.resolve()),
        },
        "configuration": {
            "mass_range_gev": list(args.mass_range),
            "analysis_window_gev": list(window),
            "bin_width_gev": args.bin_width,
            "isolated_max": args.isolated_max,
            "mc_normalization_factor": normalization,
            "background_fraction_of_observed_candidates": fraction,
            "background_shape": "DY template scaled to the requested background integral",
        },
        "integrals_in_analysis_window": {
            "data": data_window,
            "dy_mc": dy_window,
            "global_background_proxy": background_target,
            "dy_plus_background_proxy": dy_window + background_target,
            "data_over_dy": (None if dy_window == 0.0 else data_window / dy_window),
            "data_over_dy_plus_background": (
                None if dy_window + background_target == 0.0
                else data_window / (dy_window + background_target)
            ),
        },
        "selected_skim_entries": {"data": data_entries, "mc": mc_entries},
        "plots": plots,
        "caveats": [
            "The fractional background fixes only an integral; its plotted DY-like shape is a proxy.",
            "The wide-window ABCD skim stores pre-momentum-correction masses and selections.",
            "Nominal reconstruction, isolation, and trigger SFs are reconstructed from skim kinematics.",
            "Pileup weighting is disabled in the nominal measurement configuration.",
            "Use a measurement-level histogram production for an exact post-correction comparison.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"output": str(output), **summary["integrals_in_analysis_window"],
                      "plots": len(plots)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
