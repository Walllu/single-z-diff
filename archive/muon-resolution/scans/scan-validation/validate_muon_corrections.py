#!/usr/bin/env python3
"""Validate fitted muon corrections against data across eta-pair regions."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCANS = HERE.parent
sys.path.insert(0, str(SCANS))

import ROOT  # noqa: E402
import scan_zmumu_scale_resolution as scan  # noqa: E402
from muon_pair_correction import DEFAULT_PARAMETERS_FILE, Muon, correct_muon_pair  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

DEFAULT_SEEDS = [314159, 271828, 161803, 141421]
COLORS = [ROOT.kAzure + 2, ROOT.kOrange + 7, ROOT.kGreen + 2, ROOT.kMagenta + 1,
          ROOT.kCyan + 2, ROOT.kRed + 1]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    raw = SCANS.parents[2] / "HackathonDataRaw"
    parser.add_argument("--data-dir", type=Path, default=raw / "2016H")
    parser.add_argument("--mc-dir", type=Path,
                        default=raw / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8")
    parser.add_argument("--parameters", type=Path, default=DEFAULT_PARAMETERS_FILE)
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--regions", nargs="+", default=["all"])
    parser.add_argument("--seeds", nargs="+", type=int, default=DEFAULT_SEEDS)
    parser.add_argument("--max-files", type=int, default=-1)
    parser.add_argument("--max-events", type=int, default=-1)
    return parser.parse_args()


def corrected_histograms(events: list[dict], region: str, seed: int, parameters: Path):
    hist = scan.make_histograms(f"corrected_seed{seed}", region)
    accepted = 0
    for event in events:
        lead = Muon(*event["lead"]); sublead = Muon(*event["sublead"])
        result = correct_muon_pair(lead, sublead, parameters, apply_correction=True,
                                   rng_seed=seed, event_key=event["event_key"])
        first, second = result["muons"]
        lead4 = scan.p4(first["pt"], first["eta"], first["phi"], first["mass"])
        sublead4 = scan.p4(second["pt"], second["eta"], second["phi"], second["mass"])
        z = lead4 + sublead4
        low, high = scan.SELECTION["fit_mass_range_gev"]
        if not low < z.M() < high:
            continue
        accepted += 1
        hist["mass"].Fill(z.M()); hist["lead_pt"].Fill(lead4.Pt())
        hist["sublead_pt"].Fill(sublead4.Pt())
        hist["muon_eta"].Fill(first["eta"]); hist["muon_eta"].Fill(second["eta"])
        hist["z_pt"].Fill(z.Pt()); hist["nvertices"].Fill(event["nvertices"])
    return hist, accepted


def normalized(data, source, name: str):
    result = source.Clone(name); result.SetDirectory(0)
    if result.Integral() > 0:
        result.Scale(data.Integral() / result.Integral())
    return result


def draw_validation(data, uncorrected, corrected: dict[int, object], name: str,
                    region: str, output: Path) -> None:
    canvas = ROOT.TCanvas("validation_" + region + name, "", 850, 820)
    canvas.SetFillColor(ROOT.kWhite)
    top = ROOT.TPad("top_" + region + name, "", 0, 0.31, 1, 1)
    bottom = ROOT.TPad("bottom_" + region + name, "", 0, 0, 1, 0.31)
    for pad in (top, bottom):
        pad.SetFillColor(ROOT.kWhite)
    top.SetBottomMargin(0.02); bottom.SetTopMargin(0.03); bottom.SetBottomMargin(0.32)
    top.Draw(); bottom.Draw(); top.cd()

    d = data.Clone("draw_data_" + region + name); d.SetDirectory(0)
    u = normalized(d, uncorrected, "draw_uncorrected_" + region + name)
    cs = {seed: normalized(d, value, f"draw_corrected_{seed}_{region}{name}")
          for seed, value in corrected.items()}
    d.SetMarkerStyle(20); d.SetMarkerSize(0.7); d.SetLineColor(ROOT.kBlack)
    u.SetLineColor(ROOT.kGray + 2); u.SetLineStyle(2); u.SetLineWidth(3)
    maximum = max([d.GetMaximum(), u.GetMaximum()] + [h.GetMaximum() for h in cs.values()] + [1.0]) * 1.38
    u.SetMaximum(maximum); u.SetMinimum(0); u.GetXaxis().SetLabelSize(0)
    u.GetYaxis().SetTitle("Events (each MC shape normalized)"); u.Draw("HIST")
    for color, (seed, hist) in zip(COLORS, cs.items()):
        hist.SetLineColor(color); hist.SetLineWidth(2); hist.Draw("HIST SAME")
    d.Draw("E1 SAME")
    legend = ROOT.TLegend(0.56, 0.57, 0.89, 0.88); legend.SetBorderSize(0); legend.SetFillStyle(0)
    legend.AddEntry(d, "Data", "lep"); legend.AddEntry(u, "Uncorrected MC", "l")
    for seed, hist in cs.items():
        legend.AddEntry(hist, f"Corrected MC, seed {seed}", "l")
    legend.Draw()
    label = ROOT.TLatex(); label.SetNDC(); label.SetTextSize(0.038)
    label.DrawLatex(0.13, 0.92, f"{scan.abbreviated_region(region)} correction validation")

    bottom.cd()
    ratios = []
    ratio_u = d.Clone("ratio_uncorrected_" + region + name); ratio_u.SetDirectory(0); ratio_u.Divide(u)
    ratio_u.SetTitle(""); ratio_u.GetYaxis().SetTitle("Data / MC"); ratio_u.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio_u.GetYaxis().SetNdivisions(505); ratio_u.GetYaxis().SetTitleSize(0.09)
    ratio_u.GetYaxis().SetTitleOffset(0.50); ratio_u.GetYaxis().SetLabelSize(0.08)
    ratio_u.GetXaxis().SetTitleSize(0.105); ratio_u.GetXaxis().SetLabelSize(0.09)
    ratio_u.SetLineColor(ROOT.kGray + 2); ratio_u.SetMarkerColor(ROOT.kGray + 2)
    ratio_u.SetMarkerStyle(24); ratio_u.SetMarkerSize(0.55); ratio_u.Draw("E1")
    ratios.append(ratio_u)
    for color, (seed, hist) in zip(COLORS, cs.items()):
        ratio = d.Clone(f"ratio_corrected_{seed}_{region}{name}"); ratio.SetDirectory(0); ratio.Divide(hist)
        ratio.SetLineColor(color); ratio.SetMarkerColor(color); ratio.SetMarkerStyle(20); ratio.SetMarkerSize(0.45)
        ratio.Draw("E1 SAME"); ratios.append(ratio)
    line = ROOT.TLine(ratio_u.GetXaxis().GetXmin(), 1.0, ratio_u.GetXaxis().GetXmax(), 1.0)
    line.SetLineStyle(2); line.Draw()
    canvas.SaveAs(str(output / f"{region}_{name}_validation.png"))
    canvas.SaveAs(str(output / f"{region}_{name}_validation.pdf"))


def main() -> int:
    start = time.perf_counter(); args = arguments()
    payload = json.loads(args.parameters.read_text())
    requested = list(payload["regions"]) if args.regions == ["all"] else args.regions
    unknown = sorted(set(requested) - set(scan.ETA_PAIR_REGIONS))
    if unknown:
        raise SystemExit(f"Unknown regions: {unknown}")
    label = f"{payload['configuration'].get('histogram_configuration_label', 'scan')}_seeds_{len(args.seeds)}"
    output = args.output_dir / label; output.mkdir(parents=True, exist_ok=True)
    data_files = scan.discover(args.data_dir, args.max_files)
    mc_files = scan.discover(args.mc_dir, args.max_files)
    data_events, data_flow, data_total, data_read = scan.read_events(
        data_files, args.max_events, False, scan.RNG_SEEDS["nominal"])
    mc_events, mc_flow, mc_total, mc_read = scan.read_events(
        mc_files, args.max_events, True, scan.RNG_SEEDS["nominal"])
    summary = {"configuration": {"parameters": str(args.parameters.resolve()), "seeds": args.seeds,
                                  "regions": requested, "normalization": "each MC histogram to data integral"},
               "inputs": {"data_files": len(data_files), "mc_files": len(mc_files),
                          "data_entries": data_total, "mc_entries": mc_total,
                          "data_processed": data_read, "mc_processed": mc_read},
               "cutflows": {"data": data_flow, "mc": mc_flow}, "regions": {}}
    root_output = ROOT.TFile(str(output / "validation_histograms.root"), "RECREATE")
    for region in requested:
        region_output = output / scan.abbreviated_region(region); region_output.mkdir(parents=True, exist_ok=True)
        data_hist, data_count = scan.fill(data_events[region], "validation_data", region, 0, 0, False)
        uncorr_hist, uncorr_count = scan.fill(mc_events[region], "validation_uncorrected", region, 0, 0, True)
        corrected, accepted = {}, {}
        for seed in args.seeds:
            corrected[seed], accepted[seed] = corrected_histograms(mc_events[region], region, seed, args.parameters)
        for name in data_hist:
            draw_validation(data_hist[name], uncorr_hist[name],
                            {seed: hist[name] for seed, hist in corrected.items()},
                            name, region, region_output)
        metrics = {"uncorrected": dict(zip(("deviance", "ndof", "normalization"),
                   scan.objective(data_hist["mass"], uncorr_hist["mass"], "poisson")))}
        for seed, hist in corrected.items():
            metrics[str(seed)] = dict(zip(("deviance", "ndof", "normalization"),
                                      scan.objective(data_hist["mass"], hist["mass"], "poisson")))
        summary["regions"][region] = {"data": data_count, "uncorrected_mc": uncorr_count,
                                      "corrected_mc": accepted, "mass_metrics": metrics}
        root_output.cd()
        for collection in [data_hist, uncorr_hist, *corrected.values()]:
            for hist in collection.values(): hist.Write()
    root_output.Close(); summary["timing_seconds"] = time.perf_counter() - start
    (output / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"output": str(output), "timing_seconds": summary["timing_seconds"],
                      "regions": len(summary["regions"]), "seeds": args.seeds}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
