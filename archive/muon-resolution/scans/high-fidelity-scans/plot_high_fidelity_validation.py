#!/usr/bin/env python3
"""Plot independent-split validation with seed envelopes, ratios, and pulls."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from array import array
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCANS = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(SCANS))

import ROOT  # noqa: E402
import high_fidelity_scan as hf  # noqa: E402
import scan_zmumu_scale_resolution as baseline  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

OBSERVABLES = ("mass", "lead_pt", "sublead_pt", "muon_eta", "z_pt", "nvertices")
SEED_COLORS = [ROOT.kAzure + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
               ROOT.kCyan + 2, ROOT.kViolet + 1, ROOT.kOrange + 1]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="High-fidelity summary.json")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--regions", nargs="+", default=["all"])
    parser.add_argument("--subset", choices=("evaluation", "training", "all"), default="evaluation")
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--mc-dir", type=Path, default=None)
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--render-only", action="store_true",
                        help="Redraw from an existing validation_histograms.root without reading events")
    parser.add_argument("--single-process-render", action="store_true", help=argparse.SUPPRESS)
    return parser.parse_args()


def observable_edges(region: str, name: str, mass_edges: list[float]) -> np.ndarray:
    if name == "mass":
        return np.asarray(mass_edges, dtype=float)
    bins, low, high, _ = baseline.histogram_config(region)[name]
    return np.linspace(low, high, bins + 1)


def observable_axis(region: str, name: str) -> str:
    if name == "mass":
        return "m_{#mu#mu} [GeV]"
    return baseline.histogram_config(region)[name][3]


def event_values(events: list[dict], scale: float = 0.0, resolution: float = 0.0,
                 seed: int | None = None) -> dict[str, np.ndarray]:
    result = {name: [] for name in OBSERVABLES}
    low, high = baseline.SELECTION["fit_mass_range_gev"]
    for event in events:
        if seed is None:
            g1 = g2 = 0.0
        else:
            g1, g2 = hf.keyed_gaussians(seed, event["event_key"])
        lead, sublead = event["lead"], event["sublead"]
        factor1 = 1.0 + scale + resolution * g1
        factor2 = 1.0 + scale + resolution * g2
        pt1, pt2 = lead[0] * factor1, sublead[0] * factor2
        mass = baseline.dimuon_mass(lead, sublead, factor1, factor2)
        if not low < mass < high:
            continue
        zpx = pt1 * math.cos(lead[2]) + pt2 * math.cos(sublead[2])
        zpy = pt1 * math.sin(lead[2]) + pt2 * math.sin(sublead[2])
        result["mass"].append(mass)
        result["lead_pt"].append(pt1); result["sublead_pt"].append(pt2)
        result["muon_eta"].extend((lead[1], sublead[1]))
        result["z_pt"].append(math.hypot(zpx, zpy))
        result["nvertices"].append(event["nvertices"])
    return {name: np.asarray(values, dtype=float) for name, values in result.items()}


def template(values: np.ndarray, edges: np.ndarray) -> hf.Template:
    counts, _ = np.histogram(values, bins=edges)
    content = counts.astype(float)
    return hf.Template(edges, content, content.copy(), float(content.sum()))


def normalize(data: hf.Template, mc: hf.Template) -> hf.Template:
    alpha = data.content.sum() / mc.content.sum() if mc.content.sum() > 0 else 0.0
    return hf.Template(mc.edges, alpha * mc.content, alpha * alpha * mc.variance,
                       alpha * mc.accepted)


def graph_band(edges: np.ndarray, central: np.ndarray, low: np.ndarray, high: np.ndarray,
               fill_color: int, fill_alpha: float = 0.28):
    centers = 0.5 * (edges[:-1] + edges[1:])
    halfwidths = 0.5 * np.diff(edges)
    graph = ROOT.TGraphAsymmErrors(len(centers), array("d", centers.tolist()),
                                   array("d", central.tolist()),
                                   array("d", halfwidths.tolist()), array("d", halfwidths.tolist()),
                                   array("d", np.maximum(0, central - low).tolist()),
                                   array("d", np.maximum(0, high - central).tolist()))
    # Solid pale fill avoids intermittent alpha-compositing corruption in the
    # macOS ROOT PNG backend during long batch rendering jobs.
    graph.SetFillColor(ROOT.kOrange - 9); graph.SetLineColor(fill_color)
    return graph


def display_hist(name: str, edges: np.ndarray, content: np.ndarray, errors: np.ndarray | None = None):
    hist = ROOT.TH1D(name, "", len(edges) - 1, array("d", edges.tolist()))
    hist.SetDirectory(0)
    for i, value in enumerate(content, 1):
        hist.SetBinContent(i, float(value))
        if errors is not None:
            hist.SetBinError(i, float(errors[i - 1]))
    return hist


def template_from_root(histogram) -> hf.Template:
    edges = np.asarray([histogram.GetBinLowEdge(i) for i in range(1, histogram.GetNbinsX() + 2)],
                       dtype=float)
    content = np.asarray([histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX() + 1)],
                         dtype=float)
    variance = np.asarray([histogram.GetBinError(i) ** 2 for i in range(1, histogram.GetNbinsX() + 1)],
                          dtype=float)
    return hf.Template(edges, content, variance, float(content.sum()))


def render_existing(root_path: Path, regions: list[str], seeds: list[int], output: Path,
                    subset: str) -> None:
    source = ROOT.TFile.Open(str(root_path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open render source {root_path}")
    try:
        for region in regions:
            region_output = output / baseline.abbreviated_region(region)
            region_output.mkdir(parents=True, exist_ok=True)
            for name in OBSERVABLES:
                data = template_from_root(source.Get(f"data_{region}_{name}"))
                uncorrected = template_from_root(source.Get(f"uncorrected_{region}_{name}"))
                corrected = {seed: template_from_root(source.Get(f"corrected_{seed}_{region}_{name}"))
                             for seed in seeds}
                draw(region, name, data, uncorrected, corrected, region_output, subset)
    finally:
        source.Close()


def render_isolated(summary: Path, regions: list[str], subset: str, output: Path) -> None:
    """Render each region in a fresh ROOT process to avoid macOS PNG corruption."""
    for region in regions:
        command = [sys.executable, str(Path(__file__).resolve()), str(summary.resolve()),
                   "--render-only", "--single-process-render", "--subset", subset,
                   "--regions", region, "--output-dir", str(output.resolve())]
        subprocess.run(command, check=True)


def safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan, dtype=float),
                     where=denominator > 0)


def finite_envelope(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return mean/min/max without warnings for bins undefined in every seed."""
    central = np.zeros(values.shape[1], dtype=float)
    low = np.zeros_like(central); high = np.zeros_like(central)
    for index in range(values.shape[1]):
        finite = values[:, index][np.isfinite(values[:, index])]
        if finite.size:
            central[index] = finite.mean(); low[index] = finite.min(); high[index] = finite.max()
    return central, low, high


def draw(region: str, name: str, data: hf.Template, uncorrected: hf.Template,
         corrected: dict[int, hf.Template], output: Path, subset: str) -> None:
    token = f"{region}_{name}"
    edges = data.edges; widths = np.diff(edges)
    unc = normalize(data, uncorrected)
    corr = {seed: normalize(data, value) for seed, value in corrected.items()}
    corrected_stack = np.vstack([value.content for value in corr.values()])
    central_counts = corrected_stack.mean(axis=0)
    lower_counts, upper_counts = corrected_stack.min(axis=0), corrected_stack.max(axis=0)
    has_seed_spread = bool(np.any(upper_counts - lower_counts > 1e-12))

    data_density = data.content / widths
    data_error_density = np.sqrt(data.variance) / widths
    unc_density = unc.content / widths
    corr_density = corrected_stack / widths
    central_density = central_counts / widths
    lower_density, upper_density = corr_density.min(axis=0), corr_density.max(axis=0)

    canvas = ROOT.TCanvas(f"hf_validation_{token}", "", 900, 920)
    canvas.SetFillColor(ROOT.kWhite)
    top = ROOT.TPad(f"top_{token}", "", 0, 0.48, 1, 1)
    ratio_pad = ROOT.TPad(f"ratio_{token}", "", 0, 0.24, 1, 0.48)
    pull_pad = ROOT.TPad(f"pull_{token}", "", 0, 0, 1, 0.24)
    top.SetBottomMargin(0.02); ratio_pad.SetTopMargin(0.03); ratio_pad.SetBottomMargin(0.04)
    pull_pad.SetTopMargin(0.03); pull_pad.SetBottomMargin(0.36)
    for pad in (top, ratio_pad, pull_pad):
        pad.SetLeftMargin(0.12); pad.SetRightMargin(0.04); pad.SetFillColor(ROOT.kWhite)
    top.Draw(); ratio_pad.Draw(); pull_pad.Draw(); top.cd()

    frame = display_hist(f"frame_{token}", edges, np.zeros_like(data.content))
    frame.SetTitle(""); frame.GetXaxis().SetLabelSize(0)
    frame.GetYaxis().SetTitle("Events / GeV" if name == "mass" else "Events / bin width")
    ymax = max(np.nanmax(data_density), np.nanmax(unc_density), np.nanmax(upper_density), 1.0)
    frame.SetMinimum(0); frame.SetMaximum(1.42 * ymax); frame.Draw("AXIS")
    corrected_band = None
    if has_seed_spread:
        corrected_band = graph_band(edges, central_density, lower_density, upper_density, ROOT.kOrange + 1)
        corrected_band.Draw("2 SAME")
    seed_hists = []
    for color, (seed, density) in zip(SEED_COLORS, zip(corr, corr_density)):
        hist = display_hist(f"seed_{seed}_{token}", edges, density)
        hist.SetLineColor(color); hist.SetLineWidth(1); hist.Draw("HIST SAME"); seed_hists.append(hist)
    central = display_hist(f"corrected_central_{token}", edges, central_density)
    central.SetLineColor(ROOT.kOrange + 7); central.SetLineWidth(3); central.Draw("HIST SAME")
    unc_hist = display_hist(f"draw_uncorrected_{token}", edges, unc_density)
    unc_hist.SetLineColor(ROOT.kGray + 2); unc_hist.SetLineStyle(2); unc_hist.SetLineWidth(3)
    unc_hist.Draw("HIST SAME")
    data_hist = display_hist(f"draw_data_{token}", edges, data_density, data_error_density)
    data_hist.SetMarkerStyle(20); data_hist.SetMarkerSize(0.65); data_hist.Draw("E1 SAME")
    legend = ROOT.TLegend(0.53, 0.62, 0.89, 0.89)
    legend.SetBorderSize(0); legend.SetFillStyle(1001); legend.SetFillColor(ROOT.kWhite)
    legend.SetTextColor(ROOT.kBlack); legend.SetTextSize(0.034)
    data_labels = {"evaluation": "evaluation data", "training": "training data",
                   "all": "full selected data"}
    legend.AddEntry(data_hist, data_labels[subset], "lep")
    legend.AddEntry(unc_hist, "uncorrected MC", "l")
    legend.AddEntry(central, "corrected seed mean" if has_seed_spread else "corrected MC (seed independent)", "l")
    if corrected_band:
        legend.AddEntry(corrected_band, "corrected seed envelope", "f")
    legend.Draw()
    label = ROOT.TLatex(); label.SetNDC(); label.SetTextSize(0.038)
    plot_labels = {"evaluation": "independent-split validation",
                   "training": "training-split diagnostic", "all": "full-sample diagnostic"}
    label.DrawLatex(0.13, 0.92, f"{baseline.abbreviated_region(region)} {plot_labels[subset]}")

    ratio_pad.cd()
    unc_ratio = safe_ratio(data.content, unc.content)
    corr_ratios = np.vstack([safe_ratio(data.content, value.content) for value in corr.values()])
    ratio_central, ratio_low, ratio_high = finite_envelope(corr_ratios)
    ratio_frame = display_hist(f"ratio_frame_{token}", edges, np.zeros_like(data.content))
    ratio_frame.SetMinimum(0.45); ratio_frame.SetMaximum(1.55); ratio_frame.GetXaxis().SetLabelSize(0)
    ratio_frame.GetYaxis().SetTitle("Data / MC"); ratio_frame.GetYaxis().SetNdivisions(505)
    ratio_frame.GetYaxis().SetTitleSize(0.10); ratio_frame.GetYaxis().SetTitleOffset(0.45)
    ratio_frame.GetYaxis().SetLabelSize(0.085); ratio_frame.Draw("AXIS")
    ratio_band = None
    if np.any(ratio_high - ratio_low > 1e-12):
        ratio_band = graph_band(edges, ratio_central, ratio_low, ratio_high, ROOT.kOrange + 1)
        ratio_band.Draw("2 SAME")
    ratio_valid = np.any(np.isfinite(corr_ratios), axis=0)
    ratio_centers = 0.5 * (edges[:-1] + edges[1:])
    ratio_corrected_graph = ROOT.TGraph(int(np.count_nonzero(ratio_valid)),
        array("d", ratio_centers[ratio_valid].tolist()),
        array("d", ratio_central[ratio_valid].tolist()))
    ratio_corrected_graph.SetLineColor(ROOT.kOrange + 7); ratio_corrected_graph.SetLineWidth(2)
    ratio_corrected_graph.SetMarkerColor(ROOT.kOrange + 7); ratio_corrected_graph.SetMarkerStyle(20)
    ratio_corrected_graph.SetMarkerSize(0.45); ratio_corrected_graph.Draw("P SAME")
    ratio_unc_hist = display_hist(f"ratio_unc_{token}", edges, np.nan_to_num(unc_ratio, nan=0.0))
    ratio_unc_hist.SetMarkerStyle(24); ratio_unc_hist.SetMarkerSize(0.5)
    ratio_unc_hist.SetLineColor(ROOT.kGray + 2); ratio_unc_hist.Draw("P SAME")
    ratio_line = ROOT.TLine(edges[0], 1, edges[-1], 1); ratio_line.SetLineStyle(2); ratio_line.Draw()

    pull_pad.cd()
    unc_pulls = np.asarray(hf.barlow_beeston(data, uncorrected)["signed_profile_pulls"])
    corr_pulls = np.vstack([hf.barlow_beeston(data, value)["signed_profile_pulls"]
                            for value in corrected.values()])
    pull_central = corr_pulls.mean(axis=0); pull_low = corr_pulls.min(axis=0); pull_high = corr_pulls.max(axis=0)
    pull_limit = max(3.5, float(np.nanmax(np.abs(np.concatenate((unc_pulls, corr_pulls.ravel()))))) * 1.15)
    pull_frame = display_hist(f"pull_frame_{token}", edges, np.zeros_like(data.content))
    pull_frame.SetMinimum(-pull_limit); pull_frame.SetMaximum(pull_limit)
    pull_frame.GetYaxis().SetTitle("BB pull"); pull_frame.GetYaxis().SetNdivisions(505)
    pull_frame.GetYaxis().SetTitleSize(0.10); pull_frame.GetYaxis().SetTitleOffset(0.45)
    pull_frame.GetYaxis().SetLabelSize(0.085)
    pull_frame.GetXaxis().SetTitle(observable_axis(region, name))
    pull_frame.GetXaxis().SetTitleSize(0.12); pull_frame.GetXaxis().SetLabelSize(0.10)
    pull_frame.Draw("AXIS")
    pull_band = None
    if np.any(pull_high - pull_low > 1e-12):
        pull_band = graph_band(edges, pull_central, pull_low, pull_high, ROOT.kOrange + 1)
        pull_band.Draw("2 SAME")
    pull_corrected = display_hist(f"pull_corrected_{token}", edges, pull_central)
    pull_corrected.SetLineColor(ROOT.kOrange + 7); pull_corrected.SetLineWidth(2)
    pull_corrected.Draw("HIST SAME")
    pull_unc = display_hist(f"pull_unc_{token}", edges, unc_pulls)
    pull_unc.SetMarkerStyle(24); pull_unc.SetMarkerSize(0.5); pull_unc.SetLineColor(ROOT.kGray + 2)
    pull_unc.Draw("P SAME")
    pull_zero = ROOT.TLine(edges[0], 0, edges[-1], 0); pull_zero.SetLineStyle(2); pull_zero.Draw()

    canvas.Modified()
    canvas.Update()
    raster = ROOT.TImage.Create()
    raster.FromPad(canvas)
    raster.WriteImage(str(output / f"{region}_{name}_validation.png"))
    canvas.SaveAs(str(output / f"{region}_{name}_validation.svg"))
    canvas.SaveAs(str(output / f"{region}_{name}_validation.pdf"))
    # ROOT keeps named canvases in a global list. Explicit closure prevents
    # later plots from inheriting/corrupting backing stores in long batch runs.
    canvas.Close()


def select_subset(train: dict, evaluation: dict, subset: str) -> dict:
    if subset == "training":
        return train
    if subset == "evaluation":
        return evaluation
    return {region: train[region] + evaluation[region] for region in train}


def main() -> int:
    start = time.perf_counter(); args = arguments()
    payload = json.loads(args.summary.read_text()); config = payload["configuration"]
    data_dir = args.data_dir or Path(config["data_directory"])
    mc_dir = args.mc_dir or Path(config["mc_directory"])
    max_files = config.get("max_files", -1) if args.max_files is None else args.max_files
    max_events = config.get("max_events", -1) if args.max_events is None else args.max_events
    available = list(payload["regions"])
    regions = available if args.regions == ["all"] else args.regions
    unknown = sorted(set(regions) - set(available))
    if unknown:
        raise SystemExit(f"Regions not present in fit summary: {unknown}")
    output = args.output_dir or args.summary.parent / f"validation-{args.subset}"
    output.mkdir(parents=True, exist_ok=True)

    seeds = [int(value) for value in config["seeds"]]
    if args.render_only:
        if len(regions) > 1 and not args.single_process_render:
            render_isolated(args.summary, regions, args.subset, output)
        else:
            render_existing(output / "validation_histograms.root", regions, seeds, output, args.subset)
        print(json.dumps({"output": str(output), "regions": len(regions), "seeds": seeds,
                          "render_only": True}, indent=2))
        return 0

    data_files = baseline.discover(data_dir, max_files); mc_files = baseline.discover(mc_dir, max_files)
    data, _, _, _ = baseline.read_events(data_files, max_events, False, seeds[0])
    data = {region: data[region] if region in regions else []
            for region in baseline.ETA_PAIR_REGIONS}
    mc, _, _, _ = baseline.read_events(mc_files, max_events, True, seeds[0])
    mc = {region: mc[region] if region in regions else []
          for region in baseline.ETA_PAIR_REGIONS}
    train_data, eval_data = hf.split_events(data, config["train_fraction"], config["split_seed"], "data")
    train_mc, eval_mc = hf.split_events(mc, config["train_fraction"], config["split_seed"], "mc")
    chosen_data = select_subset(train_data, eval_data, args.subset)
    chosen_mc = select_subset(train_mc, eval_mc, args.subset)

    numerical = {"source_summary": str(args.summary.resolve()), "subset": args.subset,
                 "seeds": seeds, "regions": {}}
    root_file = ROOT.TFile(str(output / "validation_histograms.root"), "RECREATE")
    for region in regions:
        region_output = output / baseline.abbreviated_region(region); region_output.mkdir(parents=True, exist_ok=True)
        best = payload["regions"][region]["best"]
        mass_edges = payload["regions"][region]["binning"]["edges_gev"]
        data_values = event_values(chosen_data[region])
        unc_values = event_values(chosen_mc[region])
        corrected_values = {seed: event_values(chosen_mc[region], best["scale"], best["resolution"], seed)
                            for seed in seeds}
        numerical["regions"][region] = {"observables": {}}
        for name in OBSERVABLES:
            edges = observable_edges(region, name, mass_edges)
            data_hist = template(data_values[name], edges)
            unc_hist = template(unc_values[name], edges)
            corr_hist = {seed: template(values[name], edges) for seed, values in corrected_values.items()}
            unc_metric = hf.barlow_beeston(data_hist, unc_hist)
            corr_metrics = {str(seed): hf.slim_metrics(hf.barlow_beeston(data_hist, hist))
                            for seed, hist in corr_hist.items()}
            numerical["regions"][region]["observables"][name] = {
                "data": data_hist.accepted, "uncorrected_mc": unc_hist.accepted,
                "uncorrected": hf.slim_metrics(unc_metric), "corrected_by_seed": corr_metrics,
            }
            root_file.cd()
            axis = observable_axis(region, name)
            hf.root_histogram(f"data_{region}_{name}", data_hist, axis).Write()
            hf.root_histogram(f"uncorrected_{region}_{name}", unc_hist, axis).Write()
            for seed, hist in corr_hist.items():
                hf.root_histogram(f"corrected_{seed}_{region}_{name}", hist, axis).Write()
    root_file.Close(); numerical["timing_seconds"] = time.perf_counter() - start
    (output / "validation_summary.json").write_text(json.dumps(numerical, indent=2) + "\n")
    render_isolated(args.summary, regions, args.subset, output)
    print(json.dumps({"output": str(output), "regions": len(regions), "seeds": seeds,
                      "timing_seconds": numerical["timing_seconds"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
