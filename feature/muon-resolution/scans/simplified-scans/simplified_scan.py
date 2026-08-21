#!/usr/bin/env python3
"""Simultaneous BB/BE/EE grid fit of per-muon barrel/endcap corrections.

The three unordered event categories share four physical parameters:
barrel/endcap fractional scale and added Gaussian resolution.  The objective is
the sum of shape-only, single-source Barlow--Beeston profile deviances.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import time
from array import array
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCANS = HERE.parent
HF = SCANS / "high-fidelity-scans"
sys.path[:0] = [str(SCANS), str(HF)]

import ROOT  # noqa: E402
import high_fidelity_scan as hf  # noqa: E402
import scan_zmumu_scale_resolution as baseline  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

REGIONS = ("BB", "BE", "EE")
ETA_CONFIGURATION = {"barrel_max_abs_eta": 1.4, "muon_max_abs_eta": 2.5}
RNG_SEEDS = [314159, 271828, 161803, 141421, 173205, 223607]

# A staged block-grid search avoids the 21^4 cost of a brute-force Cartesian
# scan while repeatedly profiling each barrel/endcap parameter pair against the
# other pair.  The final 3^4 local grid permits all four parameters to move.
GRID_PROFILES = {
    "production": {
        "coarse_scale": np.linspace(-0.005, 0.005, 11).round(7).tolist(),
        "coarse_resolution": np.linspace(0.0, 0.05, 11).round(7).tolist(),
        "fine_scale_step": 0.00025,
        "fine_scale_radius": 0.001,
        "fine_resolution_step": 0.00125,
        "fine_resolution_radius": 0.005,
        "block_iterations": 2,
    },
    "test": {
        "coarse_scale": [-0.004, -0.002, 0.0, 0.002, 0.004],
        "coarse_resolution": [0.0, 0.01, 0.02, 0.03, 0.04],
        "fine_scale_step": 0.001,
        "fine_scale_radius": 0.002,
        "fine_resolution_step": 0.005,
        "fine_resolution_radius": 0.01,
        "block_iterations": 1,
    },
    "smoke": {
        "coarse_scale": [-0.002, 0.0, 0.002],
        "coarse_resolution": [0.0, 0.01, 0.02],
        "fine_scale_step": 0.001,
        "fine_scale_radius": 0.001,
        "fine_resolution_step": 0.005,
        "fine_resolution_radius": 0.005,
        "block_iterations": 1,
    },
}

DYNAMIC_BINNING = {
    "range_gev": [70.0, 110.0],
    "segments": [[70.0, 82.0, 2.0], [82.0, 84.0, 1.0],
                 [84.0, 98.0, 0.5], [98.0, 100.0, 1.0],
                 [100.0, 110.0, 2.0]],
    "minimum_data_per_bin": 20,
    "minimum_mc_per_bin": 30,
}


def arguments() -> argparse.Namespace:
    raw = baseline.default_raw_data_directory()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=raw / "2016H")
    parser.add_argument("--mc-dir", type=Path,
                        default=raw / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8")
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--label", default=None)
    parser.add_argument("--grid-profile", choices=tuple(GRID_PROFILES), default="production")
    parser.add_argument("--replicas", type=int, default=3)
    parser.add_argument("--seeds", nargs="+", type=int, default=None)
    parser.add_argument("--max-files", type=int, default=-1)
    parser.add_argument("--max-events", type=int, default=-1,
                        help="Maximum entries per sample, not per file")
    parser.add_argument("--minimum-data-per-bin", type=int,
                        default=DYNAMIC_BINNING["minimum_data_per_bin"])
    parser.add_argument("--minimum-mc-per-bin", type=int,
                        default=DYNAMIC_BINNING["minimum_mc_per_bin"])
    return parser.parse_args()


def simple_region(event: dict) -> str:
    lead_barrel = abs(event["lead"][1]) < ETA_CONFIGURATION["barrel_max_abs_eta"]
    sub_barrel = abs(event["sublead"][1]) < ETA_CONFIGURATION["barrel_max_abs_eta"]
    if lead_barrel and sub_barrel:
        return "BB"
    if lead_barrel != sub_barrel:
        return "BE"
    return "EE"


def collapse(events: dict[str, list[dict]]) -> dict[str, list[dict]]:
    result = {region: [] for region in REGIONS}
    for values in events.values():
        for event in values:
            result[simple_region(event)].append(event)
    return result


def prepare(events: list[dict], seeds: list[int]) -> dict[str, np.ndarray]:
    values = hf.prepare_arrays(events, seeds)
    values["lead_barrel"] = np.fromiter(
        (abs(event["lead"][1]) < ETA_CONFIGURATION["barrel_max_abs_eta"] for event in events),
        dtype=bool, count=len(events))
    values["sublead_barrel"] = np.fromiter(
        (abs(event["sublead"][1]) < ETA_CONFIGURATION["barrel_max_abs_eta"] for event in events),
        dtype=bool, count=len(events))
    return values


def mass_matrix(values: dict[str, np.ndarray], parameters: dict[str, float],
                replica_indices: list[int] | None = None) -> np.ndarray:
    gaussian = values["gaussians"]
    indices = replica_indices if replica_indices is not None else list(range(gaussian.shape[1]))
    lead_scale = np.where(values["lead_barrel"], parameters["barrel_scale"], parameters["endcap_scale"])
    sub_scale = np.where(values["sublead_barrel"], parameters["barrel_scale"], parameters["endcap_scale"])
    lead_res = np.where(values["lead_barrel"], parameters["barrel_resolution"], parameters["endcap_resolution"])
    sub_res = np.where(values["sublead_barrel"], parameters["barrel_resolution"], parameters["endcap_resolution"])
    factor1 = 1.0 + lead_scale[:, None] + lead_res[:, None] * gaussian[:, indices, 0]
    factor2 = 1.0 + sub_scale[:, None] + sub_res[:, None] * gaussian[:, indices, 1]
    pt1 = values["lead_pt"][:, None] * factor1
    pt2 = values["sublead_pt"][:, None] * factor2
    eta1, eta2 = values["lead_eta"][:, None], values["sublead_eta"][:, None]
    phi1, phi2 = values["lead_phi"][:, None], values["sublead_phi"][:, None]
    m1, m2 = values["lead_mass"][:, None], values["sublead_mass"][:, None]
    px = pt1 * np.cos(phi1) + pt2 * np.cos(phi2)
    py = pt1 * np.sin(phi1) + pt2 * np.sin(phi2)
    pz = pt1 * np.sinh(eta1) + pt2 * np.sinh(eta2)
    energy = np.sqrt((pt1 * np.cosh(eta1)) ** 2 + m1**2)
    energy += np.sqrt((pt2 * np.cosh(eta2)) ** 2 + m2**2)
    return np.sqrt(np.maximum(0.0, energy**2 - px**2 - py**2 - pz**2))


def mc_template(values: dict[str, np.ndarray], edges: np.ndarray,
                parameters: dict[str, float], replica_indices: list[int] | None = None) -> hf.Template:
    masses = mass_matrix(values, parameters, replica_indices)
    replicas = masses.shape[1]
    bin_index = np.searchsorted(edges, masses, side="right") - 1
    valid = (masses > edges[0]) & (masses < edges[-1])
    bin_index[~valid] = -1
    content = np.zeros(len(edges) - 1)
    variance = np.zeros_like(content)
    for b in range(len(content)):
        fraction = np.count_nonzero(bin_index == b, axis=1) / replicas
        content[b] = fraction.sum()
        variance[b] = np.square(fraction).sum()
    return hf.Template(edges, content, variance, float(content.sum()))


def build_inputs(data: dict[str, list[dict]], mc: dict[str, list[dict]], seeds: list[int],
                 min_data: int, min_mc: int) -> tuple[dict, dict]:
    prepared, details = {}, {}
    for region in REGIONS:
        edges, binning = hf.adaptive_edges(data[region], mc[region], min_data, min_mc)
        prepared[region] = {
            "data": hf.data_template(data[region], edges),
            "mc": prepare(mc[region], seeds),
            "edges": edges,
        }
        details[region] = binning
    return prepared, details


def evaluate(inputs: dict, parameters: dict[str, float]) -> tuple[float, dict]:
    total, regional = 0.0, {}
    for region in REGIONS:
        template = mc_template(inputs[region]["mc"], inputs[region]["edges"], parameters)
        metrics = hf.barlow_beeston(inputs[region]["data"], template)
        regional[region] = {**hf.slim_metrics(metrics), "mc_accepted": template.accepted}
        total += metrics["objective"]
    return float(total), regional


def axis(center: float, radius: float, step: float, low: float, high: float) -> list[float]:
    start, stop = max(low, center - radius), min(high, center + radius)
    count = int(round((stop - start) / step))
    return [round(start + i * step, 8) for i in range(count + 1)]


def scan_plane(inputs: dict, current: dict[str, float], first: str, second: str,
               first_values: list[float], second_values: list[float], stage: str,
               cache: dict[tuple, tuple[float, dict]]) -> tuple[dict, list[dict]]:
    results, best = [], None
    for one in first_values:
        for two in second_values:
            parameters = dict(current); parameters[first] = one; parameters[second] = two
            key = tuple(parameters[name] for name in
                        ("barrel_scale", "barrel_resolution", "endcap_scale", "endcap_resolution"))
            if key not in cache:
                cache[key] = evaluate(inputs, parameters)
            objective, regional = cache[key]
            item = {**parameters, "objective": objective, "regional": regional, "stage": stage}
            results.append(item)
            if best is None or objective < best["objective"]:
                best = item
    assert best is not None
    return {key: best[key] for key in current}, results


def staged_scan(inputs: dict, profile: dict) -> tuple[dict, list[dict]]:
    current = {"barrel_scale": 0.0, "barrel_resolution": 0.0,
               "endcap_scale": 0.0, "endcap_resolution": 0.0}
    records: list[dict] = []
    cache: dict[tuple, tuple[float, dict]] = {}
    for iteration in range(profile["block_iterations"]):
        current, result = scan_plane(
            inputs, current, "barrel_scale", "endcap_scale",
            profile["coarse_scale"], profile["coarse_scale"], f"coarse_scale_{iteration + 1}", cache)
        records.extend(result)
        current, result = scan_plane(
            inputs, current, "barrel_resolution", "endcap_resolution",
            profile["coarse_resolution"], profile["coarse_resolution"],
            f"coarse_resolution_{iteration + 1}", cache)
        records.extend(result)

    fine_scales_b = axis(current["barrel_scale"], profile["fine_scale_radius"],
                         profile["fine_scale_step"], -0.005, 0.005)
    fine_scales_e = axis(current["endcap_scale"], profile["fine_scale_radius"],
                         profile["fine_scale_step"], -0.005, 0.005)
    current, result = scan_plane(inputs, current, "barrel_scale", "endcap_scale",
                                 fine_scales_b, fine_scales_e, "fine_scale", cache)
    records.extend(result)
    fine_res_b = axis(current["barrel_resolution"], profile["fine_resolution_radius"],
                      profile["fine_resolution_step"], 0.0, 0.05)
    fine_res_e = axis(current["endcap_resolution"], profile["fine_resolution_radius"],
                      profile["fine_resolution_step"], 0.0, 0.05)
    current, result = scan_plane(inputs, current, "barrel_resolution", "endcap_resolution",
                                 fine_res_b, fine_res_e, "fine_resolution", cache)
    records.extend(result)

    # Let all four coordinates move together over the adjacent fine-grid cells.
    neighbors = {
        "barrel_scale": axis(current["barrel_scale"], profile["fine_scale_step"],
                             profile["fine_scale_step"], -0.005, 0.005),
        "endcap_scale": axis(current["endcap_scale"], profile["fine_scale_step"],
                             profile["fine_scale_step"], -0.005, 0.005),
        "barrel_resolution": axis(current["barrel_resolution"], profile["fine_resolution_step"],
                                  profile["fine_resolution_step"], 0.0, 0.05),
        "endcap_resolution": axis(current["endcap_resolution"], profile["fine_resolution_step"],
                                  profile["fine_resolution_step"], 0.0, 0.05),
    }
    best_item = None
    for sb in neighbors["barrel_scale"]:
        for se in neighbors["endcap_scale"]:
            for rb in neighbors["barrel_resolution"]:
                for re in neighbors["endcap_resolution"]:
                    parameters = {"barrel_scale": sb, "barrel_resolution": rb,
                                  "endcap_scale": se, "endcap_resolution": re}
                    key = (sb, rb, se, re)
                    if key not in cache:
                        cache[key] = evaluate(inputs, parameters)
                    objective, regional = cache[key]
                    item = {**parameters, "objective": objective, "regional": regional,
                            "stage": "joint_local"}
                    records.append(item)
                    if best_item is None or objective < best_item["objective"]:
                        best_item = item
    assert best_item is not None
    return best_item, records


def normalized(template: hf.Template, factor: float) -> hf.Template:
    return hf.Template(template.edges, template.content * factor,
                       template.variance * factor**2, template.accepted)


def save_canvas(canvas, stem: Path, prefer_pad_raster: bool = False) -> None:
    """Export a stable raster plus vector copies and release ROOT backing stores."""
    canvas.Modified(); canvas.Update()
    png, svg, pdf = (stem.with_suffix(suffix) for suffix in (".png", ".svg", ".pdf"))
    # ROOT can append pages when reusing an existing PDF, so plot regeneration
    # explicitly replaces artifacts owned by this script.
    for filename in (png, svg, pdf):
        filename.unlink(missing_ok=True)
    if prefer_pad_raster:
        raster = ROOT.TImage.Create(); raster.FromPad(canvas); raster.WriteImage(str(png))
    canvas.SaveAs(str(svg)); canvas.SaveAs(str(pdf))
    # Cocoa ROOT occasionally corrupts the backing store of a direct PNG while
    # leaving vector output correct. On macOS, rasterize that vector output.
    sips = shutil.which("sips")
    if prefer_pad_raster:
        # The raster was deliberately captured before vector export above.
        # Do not overwrite it from a Cocoa pad whose backing store may have
        # changed during SVG/PDF printing.
        pass
    elif sips:
        subprocess.run([sips, "-s", "format", "png", str(pdf), "--out", str(png)],
                       check=True, stdout=subprocess.DEVNULL)
    else:
        raster = ROOT.TImage.Create(); raster.FromPad(canvas); raster.WriteImage(str(png))
    canvas.Close()


def draw_validation(region: str, data: hf.Template, nominal: hf.Template,
                    corrected: hf.Template, best: dict, output: Path,
                    observable: str = "mass", axis_title: str = "m_{#mu#mu} [GeV]") -> None:
    output.mkdir(parents=True, exist_ok=True)
    nominal_metrics = hf.barlow_beeston(data, nominal)
    corrected_metrics = hf.barlow_beeston(data, corrected)
    nominal_n = normalized(nominal, nominal_metrics["normalization"])
    corrected_n = normalized(corrected, corrected_metrics["normalization"])
    token = f"{region}_{observable}"
    d = hf.root_histogram(f"draw_data_{token}", data, axis_title)
    n = hf.root_histogram(f"draw_nominal_{token}", nominal_n, axis_title)
    c = hf.root_histogram(f"draw_corrected_{token}", corrected_n, axis_title)
    canvas = ROOT.TCanvas(f"validation_{token}", "", 850, 900)
    canvas.SetFillColor(ROOT.kWhite); canvas.SetFillStyle(1001)
    top = ROOT.TPad(f"top_{token}", "", 0, 0.43, 1, 1)
    ratio_pad = ROOT.TPad(f"ratio_{token}", "", 0, 0.215, 1, 0.43)
    pull_pad = ROOT.TPad(f"pull_{token}", "", 0, 0, 1, 0.215)
    for pad in (top, ratio_pad, pull_pad):
        pad.SetFillColor(ROOT.kWhite); pad.SetFillStyle(1001)
        pad.SetFrameFillColor(ROOT.kWhite); pad.SetFrameFillStyle(1001)
        pad.SetLeftMargin(0.10); pad.SetRightMargin(0.04)
    top.SetBottomMargin(0.02); ratio_pad.SetTopMargin(0.03); ratio_pad.SetBottomMargin(0.03)
    pull_pad.SetTopMargin(0.03); pull_pad.SetBottomMargin(0.37)
    top.Draw(); ratio_pad.Draw(); pull_pad.Draw(); top.cd()
    d.SetMarkerStyle(20); d.SetMarkerSize(0.8); d.SetLineColor(ROOT.kBlack)
    n.SetLineColor(ROOT.kGray + 2); n.SetLineStyle(2); n.SetLineWidth(2)
    c.SetLineColor(ROOT.kAzure + 2); c.SetLineWidth(3)
    maximum = 1.32 * max(d.GetMaximum(), n.GetMaximum(), c.GetMaximum(), 1)
    c.SetMaximum(maximum); c.SetMinimum(0); c.GetXaxis().SetLabelSize(0)
    c.GetYaxis().SetTitle("Events (shape normalized)")
    c.Draw("HIST"); n.Draw("HIST SAME"); d.Draw("E1 SAME")
    # Text-only key avoids intermittent Cocoa ROOT TLegend backing-store
    # corruption in long batch rasterization while retaining an explicit key.
    legend_texts = []
    for y, color, text_value in ((0.86, ROOT.kBlack, "Data (points)"),
                                 (0.81, ROOT.kGray + 2, "Uncorrected MC (dashed)"),
                                 (0.76, ROOT.kAzure + 2, "Best corrected MC (solid)")):
        entry = ROOT.TLatex(); entry.SetNDC(); entry.SetTextSize(0.030)
        entry.SetTextColor(color); entry.DrawLatex(0.60, y, text_value); legend_texts.append(entry)
    label = ROOT.TLatex(); label.SetNDC(); label.SetTextSize(0.033)
    label.DrawLatex(0.13, 0.92, f"{region}: shared barrel/endcap fit")
    label.DrawLatex(0.13, 0.87,
                    f"s_{{B}}={100*best['barrel_scale']:+.3f}%, r_{{B}}={100*best['barrel_resolution']:.3f}%")
    label.DrawLatex(0.13, 0.82,
                    f"s_{{E}}={100*best['endcap_scale']:+.3f}%, r_{{E}}={100*best['endcap_resolution']:.3f}%")

    ratio_pad.cd()
    ratio_c = d.Clone(f"ratio_c_{token}"); ratio_c.Divide(c)
    ratio_n = d.Clone(f"ratio_n_{token}"); ratio_n.Divide(n)
    ratio_c.SetTitle(""); ratio_c.GetYaxis().SetTitle("Data / MC")
    ratio_c.GetYaxis().SetRangeUser(0.55, 1.45); ratio_c.GetYaxis().SetNdivisions(505)
    ratio_c.GetYaxis().SetTitleSize(0.12); ratio_c.GetYaxis().SetTitleOffset(0.38)
    ratio_c.GetYaxis().SetLabelSize(0.095); ratio_c.GetXaxis().SetLabelSize(0)
    ratio_c.SetMarkerColor(ROOT.kAzure + 2); ratio_c.SetLineColor(ROOT.kAzure + 2)
    ratio_n.SetMarkerStyle(24); ratio_n.SetMarkerColor(ROOT.kGray + 2); ratio_n.SetLineColor(ROOT.kGray + 2)
    ratio_c.Draw("E1"); ratio_n.Draw("E1 SAME")
    line1 = ROOT.TLine(data.edges[0], 1, data.edges[-1], 1); line1.SetLineStyle(2); line1.Draw()

    pull_pad.cd()
    pulls = ROOT.TH1D(f"pulls_{token}", f";{axis_title};Profile pull",
                      len(data.edges) - 1, array("d", data.edges.tolist()))
    pulls.SetDirectory(0)
    for i, value in enumerate(corrected_metrics["signed_profile_pulls"], 1):
        pulls.SetBinContent(i, value)
    pulls.SetFillColor(ROOT.kAzure - 9); pulls.SetLineColor(ROOT.kAzure + 2)
    pulls.SetMinimum(-4); pulls.SetMaximum(4); pulls.GetYaxis().SetNdivisions(505)
    pulls.GetYaxis().SetTitleSize(0.12); pulls.GetYaxis().SetTitleOffset(0.38)
    pulls.GetYaxis().SetLabelSize(0.095); pulls.GetXaxis().SetTitleSize(0.14)
    pulls.GetXaxis().SetLabelSize(0.11); pulls.Draw("HIST")
    guide_lines, guide_labels = [], []
    for sigma, color, style in ((1, ROOT.kGreen + 2, 3),
                                (2, ROOT.kOrange + 7, 7),
                                (3, ROOT.kRed + 1, 9)):
        for sign in (-1, 1):
            line = ROOT.TLine(data.edges[0], sign * sigma, data.edges[-1], sign * sigma)
            line.SetLineColor(color); line.SetLineStyle(style); line.SetLineWidth(1); line.Draw()
            guide_lines.append(line)
        text = ROOT.TLatex(data.edges[-1] - 0.012 * (data.edges[-1] - data.edges[0]),
                           sigma + 0.08, f"#pm{sigma}#sigma")
        text.SetTextAlign(31); text.SetTextSize(0.075); text.SetTextColor(color); text.Draw()
        guide_labels.append(text)
    zero = ROOT.TLine(data.edges[0], 0, data.edges[-1], 0)
    zero.SetLineColor(ROOT.kBlack); zero.SetLineStyle(2); zero.Draw()
    save_canvas(canvas, output / f"{region}_{observable}_ratio_pull", prefer_pad_raster=True)


def draw_profile_planes(records: list[dict], best: dict, output: Path,
                        stages: tuple[str, ...] = ("fine_scale", "fine_resolution")) -> None:
    for stage, xkey, ykey, xlabel, ylabel in (
        ("fine_scale", "barrel_scale", "endcap_scale", "barrel scale", "endcap scale"),
        ("fine_resolution", "barrel_resolution", "endcap_resolution",
         "barrel added resolution", "endcap added resolution"),
    ):
        if stage not in stages:
            continue
        selected = [item for item in records if item["stage"] == stage]
        if not selected:
            continue
        xs, ys = sorted({i[xkey] for i in selected}), sorted({i[ykey] for i in selected})
        if len(xs) < 2 or len(ys) < 2:
            continue
        hist = ROOT.TH2D(f"surface_{stage}", f";{xlabel};{ylabel};#Delta deviance",
                         len(xs), xs[0] - (xs[1]-xs[0])/2, xs[-1] + (xs[1]-xs[0])/2,
                         len(ys), ys[0] - (ys[1]-ys[0])/2, ys[-1] + (ys[1]-ys[0])/2)
        hist.SetDirectory(0)
        minimum = min(item["objective"] for item in selected)
        for item in selected:
            hist.Fill(item[xkey], item[ykey], item["objective"] - minimum)
        hist.SetMinimum(0); hist.SetMaximum(30)
        canvas = ROOT.TCanvas(f"canvas_{stage}", "", 850, 700)
        canvas.SetFillColor(ROOT.kWhite); canvas.SetRightMargin(0.17)
        ROOT.gPad.SetFillColor(ROOT.kWhite)
        hist.Draw("COLZ")
        marker = ROOT.TMarker(best[xkey], best[ykey], 20); marker.SetMarkerColor(ROOT.kRed + 1)
        marker.SetMarkerSize(1.5); marker.Draw("SAME")
        save_canvas(canvas, output / f"{stage}_landscape", prefer_pad_raster=True)


def main() -> int:
    start = time.perf_counter(); args = arguments()
    seeds = args.seeds if args.seeds else RNG_SEEDS[:args.replicas]
    if not seeds:
        raise SystemExit("At least one seed/replica is required")
    label = args.label or f"full_{args.grid_profile}_replicas{len(seeds)}"
    output = args.output_dir / label; output.mkdir(parents=True, exist_ok=True)

    data_files = baseline.discover(args.data_dir, args.max_files)
    mc_files = baseline.discover(args.mc_dir, args.max_files)
    extraction_start = time.perf_counter()
    raw_data, data_flow, data_entries, data_processed = baseline.read_events(
        data_files, args.max_events, False, seeds[0])
    raw_mc, mc_flow, mc_entries, mc_processed = baseline.read_events(
        mc_files, args.max_events, True, seeds[0])
    data, mc = collapse(raw_data), collapse(raw_mc)
    del raw_data, raw_mc
    inputs, binning = build_inputs(data, mc, seeds, args.minimum_data_per_bin,
                                   args.minimum_mc_per_bin)
    extraction_seconds = time.perf_counter() - extraction_start

    zero = {"barrel_scale": 0.0, "barrel_resolution": 0.0,
            "endcap_scale": 0.0, "endcap_resolution": 0.0}
    zero_objective, zero_regional = evaluate(inputs, zero)
    scan_start = time.perf_counter()
    best, records = staged_scan(inputs, GRID_PROFILES[args.grid_profile])
    scan_seconds = time.perf_counter() - scan_start

    validation = {}
    root_file = ROOT.TFile(str(output / "simplified_histograms.root"), "RECREATE")
    for region in REGIONS:
        nominal = mc_template(inputs[region]["mc"], inputs[region]["edges"], zero)
        corrected = mc_template(inputs[region]["mc"], inputs[region]["edges"], best)
        data_hist = inputs[region]["data"]
        nominal_metrics = hf.barlow_beeston(data_hist, nominal)
        corrected_metrics = hf.barlow_beeston(data_hist, corrected)
        validation[region] = {
            "counts": {"data": len(data[region]), "mc": len(mc[region])},
            "nominal": hf.slim_metrics(nominal_metrics),
            "corrected": hf.slim_metrics(corrected_metrics),
            "delta_vs_nominal": nominal_metrics["objective"] - corrected_metrics["objective"],
        }
        draw_validation(region, data_hist, nominal, corrected, best, output / region)
        for name, template in (("data", data_hist), ("nominal_mc", nominal),
                               ("corrected_mc", corrected)):
            hf.root_histogram(f"{name}_{region}_mass", template).Write()
    root_file.Close()
    draw_profile_planes(records, best, output)

    summary = {
        "schema_version": 1,
        "configuration": {
            "fit": "simultaneous unordered BB/BE/EE shape fit",
            "parameters": ["barrel_scale", "barrel_resolution", "endcap_scale", "endcap_resolution"],
            "metric": "single-source analytic Barlow-Beeston profile deviance",
            "search": "iterated two-dimensional block grids plus local four-dimensional grid",
            "grid_profile": args.grid_profile,
            "grid": GRID_PROFILES[args.grid_profile],
            "eta": ETA_CONFIGURATION,
            "selection": baseline.SELECTION,
            "dynamic_binning": DYNAMIC_BINNING,
            "seeds": seeds,
            "replicas_averaged": len(seeds),
            "data_directory": str(args.data_dir.resolve()),
            "mc_directory": str(args.mc_dir.resolve()),
            "max_files": args.max_files, "max_events": args.max_events,
        },
        "inputs": {
            "data_files": len(data_files), "mc_files": len(mc_files),
            "data_entries": data_entries, "mc_entries": mc_entries,
            "data_processed": data_processed, "mc_processed": mc_processed,
        },
        "cutflows": {"data": data_flow, "mc": mc_flow},
        "binning": binning,
        "zero_fit": {"objective": zero_objective, "regional": zero_regional},
        "best_fit": {key: best[key] for key in
                     ("barrel_scale", "barrel_resolution", "endcap_scale", "endcap_resolution", "objective")},
        "delta_vs_zero": zero_objective - best["objective"],
        "best_regional": best["regional"],
        "validation": validation,
        "scan_points_evaluated": len(records),
        "scan_records": records,
        "timing_seconds": {"event_extraction": extraction_seconds, "scan": scan_seconds,
                           "total": time.perf_counter() - start},
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "calibration.json").write_text(json.dumps({
        "barrel_scale": best["barrel_scale"],
        "barrel_resolution": best["barrel_resolution"],
        "endcap_scale": best["endcap_scale"],
        "endcap_resolution": best["endcap_resolution"],
    }, indent=2) + "\n")
    print(json.dumps({"output": str(output), "inputs": summary["inputs"],
                      "best_fit": summary["best_fit"], "delta_vs_zero": summary["delta_vs_zero"],
                      "validation": validation, "timing_seconds": summary["timing_seconds"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
