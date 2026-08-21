#!/usr/bin/env python3
"""Derive regional nVertices weights and validate their kinematic effects."""

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
SIMPLIFIED = HERE.parent
sys.path.insert(0, str(SIMPLIFIED))

import simplified_scan as scan  # noqa: E402
import plot_simplified_observables as observables  # noqa: E402


WEIGHT_CONFIGURATION = {
    "overflow_vertex": 50,
    "minimum_data_per_group": 100,
    "minimum_mc_per_group": 300,
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="simplified scan summary.json")
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--label", default="full_production_r3")
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--mc-dir", type=Path, default=None)
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--minimum-data-per-group", type=int,
                        default=WEIGHT_CONFIGURATION["minimum_data_per_group"])
    parser.add_argument("--minimum-mc-per-group", type=int,
                        default=WEIGHT_CONFIGURATION["minimum_mc_per_group"])
    parser.add_argument("--regions", nargs="+", choices=scan.REGIONS, default=list(scan.REGIONS))
    parser.add_argument("--observables", nargs="+", choices=tuple(observables.OBSERVABLES),
                        default=list(observables.OBSERVABLES))
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--render-weights", action="store_true")
    return parser.parse_args()


def vertex_counts(events: list[dict], overflow: int) -> np.ndarray:
    result = np.zeros(overflow + 1)
    for event in events:
        result[min(max(0, int(event["nvertices"])), overflow)] += 1
    return result


def derive_groups(data_counts: np.ndarray, mc_counts: np.ndarray,
                  minimum_data: int, minimum_mc: int) -> list[dict]:
    groups, start = [], 0
    accumulated_data = accumulated_mc = 0.0
    for value in range(len(data_counts)):
        accumulated_data += data_counts[value]
        accumulated_mc += mc_counts[value]
        if accumulated_data >= minimum_data and accumulated_mc >= minimum_mc:
            groups.append({"minimum": start, "maximum": value,
                           "data_count": accumulated_data, "mc_count": accumulated_mc})
            start = value + 1; accumulated_data = accumulated_mc = 0.0
    if start < len(data_counts):
        if groups:
            groups[-1]["maximum"] = len(data_counts) - 1
            groups[-1]["data_count"] += accumulated_data
            groups[-1]["mc_count"] += accumulated_mc
        else:
            groups.append({"minimum": 0, "maximum": len(data_counts) - 1,
                           "data_count": accumulated_data, "mc_count": accumulated_mc})
    total_data, total_mc = float(data_counts.sum()), float(mc_counts.sum())
    for group in groups:
        data_probability = group["data_count"] / total_data if total_data else 0.0
        mc_probability = group["mc_count"] / total_mc if total_mc else 0.0
        group["weight"] = data_probability / mc_probability if mc_probability else 0.0
    return groups


def lookup(groups: list[dict], overflow: int) -> np.ndarray:
    result = np.ones(overflow + 1)
    for group in groups:
        result[group["minimum"]:group["maximum"] + 1] = group["weight"]
    return result


def event_weights(events: list[dict], values: np.ndarray, overflow: int) -> np.ndarray:
    return np.fromiter((values[min(max(0, int(event["nvertices"])), overflow)]
                        for event in events), dtype=float, count=len(events))


def weighted_template(value_matrix: np.ndarray, mass_matrix: np.ndarray,
                      edges: np.ndarray, weights: np.ndarray) -> scan.hf.Template:
    replicas = value_matrix.shape[1]
    accepted_mass = (mass_matrix > 70.0) & (mass_matrix < 110.0)
    indices = np.searchsorted(edges, value_matrix, side="right") - 1
    valid = accepted_mass & (value_matrix >= edges[0]) & (value_matrix < edges[-1])
    indices[~valid] = -1
    content = np.zeros(len(edges) - 1); variance = np.zeros_like(content)
    for b in range(len(content)):
        fraction = np.count_nonzero(indices == b, axis=1) / replicas
        contribution = weights * fraction
        content[b] = contribution.sum(); variance[b] = np.square(contribution).sum()
    return scan.hf.Template(edges, content, variance, float(content.sum()))


def normalized(template: scan.hf.Template, factor: float) -> scan.hf.Template:
    return scan.hf.Template(template.edges, template.content * factor,
                            template.variance * factor**2, template.accepted)


def draw_comparison(region: str, observable: str, data: scan.hf.Template,
                    before: scan.hf.Template, after: scan.hf.Template,
                    output: Path, *, title_text: str | None = None,
                    before_label: str = "Muon-corrected MC (dashed)",
                    after_label: str = "+ vertex reweighting (solid)",
                    filename_suffix: str = "pileup",
                    axis_title: str | None = None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    before_metric = scan.hf.barlow_beeston(data, before)
    after_metric = scan.hf.barlow_beeston(data, after)
    before_n = normalized(before, before_metric["normalization"])
    after_n = normalized(after, after_metric["normalization"])
    axis_title = axis_title or observables.OBSERVABLES[observable][1]
    token = f"{region}_{observable}_pileup"
    d = scan.hf.root_histogram(f"draw_data_{token}", data, axis_title)
    b = scan.hf.root_histogram(f"draw_before_{token}", before_n, axis_title)
    a = scan.hf.root_histogram(f"draw_after_{token}", after_n, axis_title)
    canvas = scan.ROOT.TCanvas(f"canvas_{token}", "", 850, 900)
    canvas.SetFillColor(scan.ROOT.kWhite); canvas.SetFillStyle(1001)
    top = scan.ROOT.TPad(f"top_{token}", "", 0, .43, 1, 1)
    ratio = scan.ROOT.TPad(f"ratio_{token}", "", 0, .215, 1, .43)
    pull = scan.ROOT.TPad(f"pull_{token}", "", 0, 0, 1, .215)
    for pad in (top, ratio, pull):
        pad.SetFillColor(scan.ROOT.kWhite); pad.SetFillStyle(1001)
        pad.SetFrameFillColor(scan.ROOT.kWhite); pad.SetLeftMargin(.10); pad.SetRightMargin(.04)
    top.SetBottomMargin(.02); ratio.SetTopMargin(.03); ratio.SetBottomMargin(.03)
    pull.SetTopMargin(.03); pull.SetBottomMargin(.37)
    top.Draw(); ratio.Draw(); pull.Draw(); top.cd()
    d.SetMarkerStyle(20); d.SetMarkerSize(.75); d.SetLineColor(scan.ROOT.kBlack)
    b.SetLineColor(scan.ROOT.kGray + 2); b.SetLineStyle(2); b.SetLineWidth(2)
    a.SetLineColor(scan.ROOT.kMagenta + 2); a.SetLineWidth(3)
    a.SetMinimum(0); a.SetMaximum(1.35 * max(d.GetMaximum(), b.GetMaximum(), a.GetMaximum(), 1))
    a.GetXaxis().SetLabelSize(0); a.GetYaxis().SetTitle("Events (shape normalized)")
    a.Draw("HIST"); b.Draw("HIST SAME"); d.Draw("E1 SAME")
    key = []
    for y, color, text in ((.86, scan.ROOT.kBlack, "Data (points)"),
                           (.81, scan.ROOT.kGray + 2, before_label),
                           (.76, scan.ROOT.kMagenta + 2, after_label)):
        item = scan.ROOT.TLatex(); item.SetNDC(); item.SetTextColor(color); item.SetTextSize(.029)
        item.DrawLatex(.57, y, text); key.append(item)
    title = scan.ROOT.TLatex(); title.SetNDC(); title.SetTextSize(.034)
    title.DrawLatex(.13, .92, title_text or f"{region}: reconstructed-vertex pileup reweighting")

    ratio.cd(); ratio_after = d.Clone(f"ratio_after_{token}"); ratio_after.Divide(a)
    ratio_before = d.Clone(f"ratio_before_{token}"); ratio_before.Divide(b)
    ratio_after.SetTitle(""); ratio_after.GetYaxis().SetTitle("Data / MC")
    ratio_after.GetYaxis().SetRangeUser(.55, 1.45); ratio_after.GetYaxis().SetNdivisions(505)
    ratio_after.GetYaxis().SetTitleSize(.12); ratio_after.GetYaxis().SetTitleOffset(.38)
    ratio_after.GetYaxis().SetLabelSize(.095); ratio_after.GetXaxis().SetLabelSize(0)
    ratio_after.SetMarkerColor(scan.ROOT.kMagenta + 2); ratio_after.SetLineColor(scan.ROOT.kMagenta + 2)
    ratio_before.SetMarkerStyle(24); ratio_before.SetMarkerColor(scan.ROOT.kGray + 2)
    ratio_before.SetLineColor(scan.ROOT.kGray + 2)
    ratio_after.Draw("E1"); ratio_before.Draw("E1 SAME")
    unity = scan.ROOT.TLine(data.edges[0], 1, data.edges[-1], 1); unity.SetLineStyle(2); unity.Draw()

    pull.cd(); pulls = scan.ROOT.TH1D(f"pull_{token}", f";{axis_title};Profile pull",
                                      len(data.edges) - 1, array("d", data.edges.tolist()))
    pulls.SetDirectory(0)
    for i, value in enumerate(after_metric["signed_profile_pulls"], 1): pulls.SetBinContent(i, value)
    pulls.SetFillColor(scan.ROOT.kMagenta - 9); pulls.SetLineColor(scan.ROOT.kMagenta + 2)
    pulls.SetMinimum(-4); pulls.SetMaximum(4); pulls.GetYaxis().SetNdivisions(505)
    pulls.GetYaxis().SetTitleSize(.12); pulls.GetYaxis().SetTitleOffset(.38)
    pulls.GetYaxis().SetLabelSize(.095); pulls.GetXaxis().SetTitleSize(.14)
    pulls.GetXaxis().SetLabelSize(.11); pulls.Draw("HIST")
    guides = []
    for sigma, color, style in ((1, scan.ROOT.kGreen + 2, 3),
                                (2, scan.ROOT.kOrange + 7, 7),
                                (3, scan.ROOT.kRed + 1, 9)):
        for sign in (-1, 1):
            line = scan.ROOT.TLine(data.edges[0], sign*sigma, data.edges[-1], sign*sigma)
            line.SetLineColor(color); line.SetLineStyle(style); line.Draw(); guides.append(line)
    zero = scan.ROOT.TLine(data.edges[0], 0, data.edges[-1], 0); zero.SetLineStyle(2); zero.Draw()
    scan.save_canvas(canvas, output / f"{region}_{observable}_{filename_suffix}_ratio_pull",
                     prefer_pad_raster=True)


def from_root(histogram) -> scan.hf.Template:
    bins = histogram.GetNbinsX()
    edges = np.asarray([histogram.GetXaxis().GetBinLowEdge(i) for i in range(1, bins + 2)])
    content = np.asarray([histogram.GetBinContent(i) for i in range(1, bins + 1)])
    variance = np.asarray([histogram.GetBinError(i)**2 for i in range(1, bins + 1)])
    return scan.hf.Template(edges, content, variance, float(content.sum()))


def draw_weight_curve(region: str, calibration: dict, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    groups = calibration["regions"][region]; overflow = calibration["configuration"]["overflow_vertex"]
    hist = scan.ROOT.TH1D(f"weights_{region}", f";number of reconstructed vertices;event weight",
                          overflow + 1, -.5, overflow + .5); hist.SetDirectory(0)
    for group in groups:
        for value in range(group["minimum"], group["maximum"] + 1):
            hist.SetBinContent(value + 1, group["weight"])
    hist.SetLineColor(scan.ROOT.kMagenta + 2); hist.SetLineWidth(3); hist.SetMinimum(0)
    hist.SetMaximum(1.20 * max(hist.GetMaximum(), 1));
    canvas = scan.ROOT.TCanvas(f"weight_canvas_{region}", "", 850, 650)
    canvas.SetFillColor(scan.ROOT.kWhite); hist.Draw("HIST")
    label = scan.ROOT.TLatex(); label.SetNDC(); label.SetTextSize(.037)
    label.DrawLatex(.13, .92, f"{region}: regional reconstructed-vertex weights")
    unity = scan.ROOT.TLine(-.5, 1, overflow + .5, 1); unity.SetLineStyle(2); unity.Draw()
    scan.save_canvas(canvas, output / f"{region}_pileup_weights", prefer_pad_raster=True)


def render(args: argparse.Namespace, output: Path) -> int:
    calibration = json.loads((output / "pileup_weights.json").read_text())
    if args.render_weights:
        for region in args.regions: draw_weight_curve(region, calibration, output / region)
        return 0
    root_file = scan.ROOT.TFile.Open(str(output / "pileup_histograms.root"), "READ")
    for region in args.regions:
        for observable in args.observables:
            objects = [root_file.Get(f"{kind}_{region}_{observable}")
                       for kind in ("data", "before", "after")]
            if any(not item for item in objects):
                raise SystemExit(f"Missing {region}/{observable} pileup histograms")
            draw_comparison(region, observable, *(from_root(item) for item in objects), output / region)
    root_file.Close(); return 0


def main() -> int:
    start = time.perf_counter(); args = arguments(); output = args.output_dir / args.label
    if args.render_only:
        return render(args, output)
    payload = json.loads(args.summary.read_text()); config = payload["configuration"]
    data_dir = args.data_dir or Path(config["data_directory"])
    mc_dir = args.mc_dir or Path(config["mc_directory"])
    max_files = config.get("max_files", -1) if args.max_files is None else args.max_files
    max_events = config.get("max_events", -1) if args.max_events is None else args.max_events
    seeds = [int(seed) for seed in config["seeds"]]; best = payload["best_fit"]
    data_files = scan.baseline.discover(data_dir, max_files); mc_files = scan.baseline.discover(mc_dir, max_files)
    raw_data, data_flow, data_entries, data_processed = scan.baseline.read_events(data_files, max_events, False, seeds[0])
    raw_mc, mc_flow, mc_entries, mc_processed = scan.baseline.read_events(mc_files, max_events, True, seeds[0])
    data, mc = scan.collapse(raw_data), scan.collapse(raw_mc); del raw_data, raw_mc
    overflow = WEIGHT_CONFIGURATION["overflow_vertex"]
    calibration = {"schema_version": 1, "method": "regional reconstructed-nVertices probability ratio",
                   "configuration": {"overflow_vertex": overflow,
                       "minimum_data_per_group": args.minimum_data_per_group,
                       "minimum_mc_per_group": args.minimum_mc_per_group}, "regions": {}}
    for region in args.regions:
        calibration["regions"][region] = derive_groups(
            vertex_counts(data[region], overflow), vertex_counts(mc[region], overflow),
            args.minimum_data_per_group, args.minimum_mc_per_group)
    output.mkdir(parents=True, exist_ok=True)
    (output / "pileup_weights.json").write_text(json.dumps(calibration, indent=2) + "\n")

    numerical = {"source_summary": str(args.summary.resolve()), "configuration": calibration["configuration"],
                 "inputs": {"data_files": len(data_files), "mc_files": len(mc_files),
                            "data_entries": data_entries, "mc_entries": mc_entries,
                            "data_processed": data_processed, "mc_processed": mc_processed},
                 "cutflows": {"data": data_flow, "mc": mc_flow}, "regions": {}}
    root_file = scan.ROOT.TFile(str(output / "pileup_histograms.root"), "RECREATE")
    for region in args.regions:
        data_values = observables.values(data[region], [seeds[0]], None)
        corrected_values = observables.values(mc[region], seeds, best)
        table = lookup(calibration["regions"][region], overflow)
        pu_weights = event_weights(mc[region], table, overflow)
        unit_weights = np.ones(len(data[region])); unit_mc = np.ones(len(mc[region]))
        numerical["regions"][region] = {"weight_groups": calibration["regions"][region],
            "mc_weight_sum": float(pu_weights.sum()), "mc_weight_sum2": float(np.square(pu_weights).sum()),
            "effective_mc_events": float(pu_weights.sum()**2 / np.square(pu_weights).sum()),
            "observables": {}}
        for observable in args.observables:
            edges = observables.edges(payload, region, observable)
            data_hist = weighted_template(data_values[observable], data_values["mass"], edges, unit_weights)
            before_hist = weighted_template(corrected_values[observable], corrected_values["mass"], edges, unit_mc)
            after_hist = weighted_template(corrected_values[observable], corrected_values["mass"], edges, pu_weights)
            before_metric = scan.hf.barlow_beeston(data_hist, before_hist)
            after_metric = scan.hf.barlow_beeston(data_hist, after_hist)
            numerical["regions"][region]["observables"][observable] = {
                "before": scan.hf.slim_metrics(before_metric), "after": scan.hf.slim_metrics(after_metric),
                "delta_vs_before": before_metric["objective"] - after_metric["objective"]}
            root_file.cd()
            for kind, hist in (("data", data_hist), ("before", before_hist), ("after", after_hist)):
                scan.hf.root_histogram(f"{kind}_{region}_{observable}", hist,
                                       observables.OBSERVABLES[observable][1]).Write()
        del data_values, corrected_values
    root_file.Close(); numerical["timing_seconds"] = {"derivation_and_histograms": time.perf_counter() - start}
    (output / "pileup_validation_summary.json").write_text(json.dumps(numerical, indent=2) + "\n")

    for region in args.regions:
        subprocess.run([sys.executable, str(Path(__file__).resolve()), str(args.summary.resolve()),
                        "--output-dir", str(args.output_dir), "--label", args.label,
                        "--render-only", "--render-weights", "--regions", region], check=True)
        for observable in args.observables:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), str(args.summary.resolve()),
                            "--output-dir", str(args.output_dir), "--label", args.label,
                            "--render-only", "--regions", region, "--observables", observable], check=True)
    print(json.dumps({"output": str(output), "regions": args.regions,
                      "observables": args.observables, "timing_seconds": time.perf_counter()-start}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
