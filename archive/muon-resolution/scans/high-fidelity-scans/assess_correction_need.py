#!/usr/bin/env python3
"""Assess whether each fitted correction is supported beyond the zero point."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from array import array
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCANS = HERE.parent
sys.path.insert(0, str(SCANS))

import ROOT  # noqa: E402
import scan_zmumu_scale_resolution as baseline  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="High-fidelity summary.json")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--minimum-training-delta", type=float, default=4.0,
                        help="Default AIC-like cost for adding scale and resolution")
    parser.add_argument("--minimum-evaluation-delta", type=float, default=0.0)
    parser.add_argument("--minimum-seed-fraction", type=float, default=0.60)
    return parser.parse_args()


def classification(train_delta: float, eval_delta: float, fraction: float,
                   args: argparse.Namespace) -> tuple[str, list[str]]:
    tests = {
        "training_improvement": train_delta > args.minimum_training_delta,
        "independent_evaluation_improvement": eval_delta > args.minimum_evaluation_delta,
        "seed_robustness": fraction >= args.minimum_seed_fraction,
    }
    failed = [name for name, passed in tests.items() if not passed]
    if not failed:
        return "supported", failed
    if not tests["independent_evaluation_improvement"] and not tests["seed_robustness"]:
        return "prefer_identity", failed
    return "inconclusive", failed


def draw_delta_grid(rows: dict[str, dict], output: Path) -> None:
    names = list(baseline.ETA_REGIONS)
    labels = [baseline.ETA_ABBREVIATIONS[name] for name in names]
    canvas = ROOT.TCanvas("hf_need_delta_grid", "", 1300, 620)
    canvas.Divide(2, 1)
    keepalive = []
    for pad_index, key, title in (
        (1, "training_delta", "Training #DeltaD = D(0)-D(best)"),
        (2, "evaluation_delta", "Evaluation #DeltaD = D(0)-D(best)"),
    ):
        canvas.cd(pad_index); ROOT.gPad.SetRightMargin(0.16); ROOT.gPad.SetLeftMargin(0.13)
        hist = ROOT.TH2D("need_" + key, f"{title};subleading-muon region;leading-muon region",
                         5, 0, 5, 5, 0, 5)
        hist.SetDirectory(0)
        for i, label in enumerate(labels, 1):
            hist.GetXaxis().SetBinLabel(i, label); hist.GetYaxis().SetBinLabel(i, label)
        values = [row[key] for row in rows.values()]
        extent = max([abs(value) for value in values] + [1.0])
        hist.SetMinimum(-extent); hist.SetMaximum(extent)
        for iy, lead in enumerate(names, 1):
            for ix, sublead in enumerate(names, 1):
                region = f"{lead}__{sublead}"
                if region in rows:
                    value = rows[region][key]
                    hist.SetBinContent(ix, iy, value if value != 0 else 1e-12)
        hist.Draw("COLZ")
        boxes, texts = [], []
        for iy, lead in enumerate(names, 1):
            for ix, sublead in enumerate(names, 1):
                region = f"{lead}__{sublead}"
                if region not in rows:
                    box = ROOT.TBox(ix - 1, iy - 1, ix, iy)
                    box.SetFillColor(ROOT.kGray + 1); box.Draw("SAME"); boxes.append(box)
                    label = "X"
                else:
                    label = f"{rows[region][key]:.1f}"
                text = ROOT.TLatex(ix - 0.5, iy - 0.5, label)
                text.SetTextAlign(22); text.SetTextSize(0.035); text.SetTextFont(62); text.Draw("SAME")
                texts.append(text)
        keepalive.extend([hist, boxes, texts])
    canvas.SaveAs(str(output / "correction_delta_grid.png"))
    canvas.SaveAs(str(output / "correction_delta_grid.pdf"))


def draw_region_summary(rows: dict[str, dict], output: Path) -> None:
    ordered = sorted(rows, key=baseline.abbreviated_region)
    x = array("d", [i + 0.5 for i in range(len(ordered))])
    train = array("d", [rows[r]["training_delta"] for r in ordered])
    evaluation = array("d", [rows[r]["evaluation_delta"] for r in ordered])
    canvas = ROOT.TCanvas("hf_need_summary", "", 1100, 650)
    frame = ROOT.TH1D("need_frame", ";eta-pair region;#DeltaD (positive favors correction)",
                      len(ordered), 0, len(ordered))
    frame.SetDirectory(0)
    for i, region in enumerate(ordered, 1):
        frame.GetXaxis().SetBinLabel(i, baseline.abbreviated_region(region))
    values = list(train) + list(evaluation)
    low, high = min(values + [0]), max(values + [0])
    margin = 0.15 * max(high - low, 1.0)
    frame.SetMinimum(low - margin); frame.SetMaximum(high + margin); frame.Draw()
    train_graph = ROOT.TGraph(len(ordered), x, train)
    train_graph.SetMarkerStyle(20); train_graph.SetMarkerColor(ROOT.kAzure + 2)
    evaluation_graph = ROOT.TGraph(len(ordered), x, evaluation)
    evaluation_graph.SetMarkerStyle(24); evaluation_graph.SetMarkerColor(ROOT.kOrange + 7)
    train_graph.Draw("P SAME"); evaluation_graph.Draw("P SAME")
    zero = ROOT.TLine(0, 0, len(ordered), 0); zero.SetLineStyle(2); zero.Draw()
    legend = ROOT.TLegend(0.68, 0.75, 0.89, 0.88); legend.SetBorderSize(0)
    legend.AddEntry(train_graph, "training", "p"); legend.AddEntry(evaluation_graph, "evaluation", "p")
    legend.Draw()
    canvas.SaveAs(str(output / "correction_delta_by_region.png"))
    canvas.SaveAs(str(output / "correction_delta_by_region.pdf"))


def main() -> int:
    args = arguments()
    payload = json.loads(args.summary.read_text())
    output = args.output_dir or args.summary.parent / "correction-need"
    output.mkdir(parents=True, exist_ok=True)
    rows = {}
    for region, result in payload["regions"].items():
        train_delta = float(result["training_delta_vs_zero"])
        eval_delta = float(result["evaluation_at_best"]["delta_vs_zero"])
        seed_deltas = {seed: float(metrics["delta_vs_zero"])
                       for seed, metrics in result["evaluation_by_seed"].items()}
        fraction = sum(value > 0 for value in seed_deltas.values()) / max(1, len(seed_deltas))
        status, failed = classification(train_delta, eval_delta, fraction, args)
        best = result["best"]
        rows[region] = {
            "region_abbreviation": baseline.abbreviated_region(region),
            "scale": float(best["scale"]), "resolution": float(best["resolution"]),
            "training_delta": train_delta, "evaluation_delta": eval_delta,
            "seed_deltas": seed_deltas, "seed_fraction_improved": fraction,
            "classification": status, "failed_criteria": failed,
            "aic_like_delta": train_delta - 4.0,
            "note": ("The AIC-like diagnostic subtracts 2k for k=2 fitted parameters; "
                     "grid boundaries and the non-negative resolution make it descriptive, not a p-value."),
        }
    assessment = {
        "source_summary": str(args.summary.resolve()),
        "criteria": {
            "minimum_training_delta": args.minimum_training_delta,
            "minimum_evaluation_delta": args.minimum_evaluation_delta,
            "minimum_seed_fraction": args.minimum_seed_fraction,
            "interpretation": "supported requires all three; prefer_identity fails evaluation and seed tests",
        },
        "regions": rows,
    }
    (output / "correction_need.json").write_text(json.dumps(assessment, indent=2) + "\n")
    with (output / "correction_need.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["region", "scale", "resolution", "training_delta", "evaluation_delta",
                         "seed_fraction_improved", "classification"])
        for region, row in sorted(rows.items(), key=lambda item: item[1]["region_abbreviation"]):
            writer.writerow([row["region_abbreviation"], row["scale"], row["resolution"],
                             row["training_delta"], row["evaluation_delta"],
                             row["seed_fraction_improved"], row["classification"]])
    draw_delta_grid(rows, output)
    draw_region_summary(rows, output)
    print(json.dumps({row["region_abbreviation"]: {
        "classification": row["classification"], "training_delta": row["training_delta"],
        "evaluation_delta": row["evaluation_delta"],
        "seed_fraction": row["seed_fraction_improved"],
    } for row in rows.values()}, indent=2))
    print(f"Wrote correction-need assessment to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
