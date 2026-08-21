#!/usr/bin/env python3
"""Recommend per-region dimuon-mass binning from a selected data subset."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import ROOT

import scan_zmumu_scale_resolution as scan

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=scan.HERE.parents[2] / "HackathonDataRaw/2016H")
    parser.add_argument("--output-dir", type=Path, default=scan.HERE / "binning-tests")
    parser.add_argument("--regions", nargs="+", default=["all"], help="Region names or 'all'")
    parser.add_argument("--max-files", type=int, default=3)
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("--minimum-events-per-bin", type=int, default=20)
    parser.add_argument("--minimum-bins", type=int, default=5)
    parser.add_argument("--maximum-bins", type=int, default=100)
    return parser.parse_args()


def rule_bins(values: np.ndarray, low: float, high: float) -> dict[str, int]:
    n = len(values)
    if n < 2:
        return {"freedman_diaconis": 1, "scott": 1, "sturges": 1}
    iqr = float(np.percentile(values, 75) - np.percentile(values, 25))
    sigma = float(np.std(values, ddof=1))
    fd_width = 2.0 * iqr * n ** (-1.0 / 3.0)
    scott_width = 3.5 * sigma * n ** (-1.0 / 3.0)
    return {
        "freedman_diaconis": max(1, math.ceil((high - low) / fd_width)) if fd_width > 0 else 1,
        "scott": max(1, math.ceil((high - low) / scott_width)) if scott_width > 0 else 1,
        "sturges": max(1, math.ceil(math.log2(n) + 1.0)),
    }


def recommendation(values: np.ndarray, low: float, high: float, args) -> dict:
    rules = rule_bins(values, low, high)
    statistics_choice = rules["freedman_diaconis"]
    population_cap = max(1, len(values) // args.minimum_events_per_bin)
    recommended = min(statistics_choice, population_cap, args.maximum_bins)
    recommended = max(args.minimum_bins, recommended) if len(values) >= args.minimum_bins else max(1, len(values))
    return {**rules, "population_cap": population_cap, "recommended": recommended,
            "recommended_width_gev": (high - low) / recommended if recommended else None}


def draw(values: np.ndarray, region: str, current: int, recommended: int,
         low: float, high: float, output: Path) -> None:
    canvas = ROOT.TCanvas("binning_" + region, "", 1200, 500); canvas.Divide(2, 1)
    histograms = []
    for pad, bins, title in ((1, current, f"Current: {current} bins"),
                             (2, recommended, f"Recommended: {recommended} bins")):
        canvas.cd(pad)
        hist = ROOT.TH1D(f"{region}_{pad}", f"{title};m_{{#mu#mu}} [GeV];Data events", bins, low, high)
        hist.SetDirectory(0); hist.SetLineColor(ROOT.kBlack); hist.SetMarkerStyle(20); hist.SetMarkerSize(0.6)
        for value in values:
            hist.Fill(float(value))
        hist.Draw("E1"); histograms.append(hist)
    canvas.SaveAs(str(output / f"{scan.abbreviated_region(region)}_mass_binning.png"))
    canvas.SaveAs(str(output / f"{scan.abbreviated_region(region)}_mass_binning.pdf"))


def main() -> int:
    args = arguments(); args.output_dir.mkdir(parents=True, exist_ok=True)
    requested = scan.ETA_PAIR_REGIONS if args.regions == ["all"] else args.regions
    unknown = sorted(set(requested) - set(scan.ETA_PAIR_REGIONS))
    if unknown:
        raise SystemExit(f"Unknown regions: {unknown}")
    regions = [r for r in requested if scan.abbreviated_region(r) not in scan.EXCLUDED_ETA_PAIR_ABBREVIATIONS]
    files = scan.discover(args.data_dir, args.max_files)
    events, cutflow, total, processed = scan.read_events(files, args.max_events, False, scan.RNG_SEEDS["nominal"])
    results = {}
    for region in regions:
        values = np.asarray([scan.dimuon_mass(e["lead"], e["sublead"], 1.0, 1.0)
                             for e in events[region]], dtype=float)
        current, low, high, _ = scan.histogram_config(region)["mass"]
        result = recommendation(values, low, high, args)
        result.update({"events": len(values), "current_bins": current,
                       "current_width_gev": (high - low) / current, "range_gev": [low, high]})
        results[region] = result
        if len(values):
            draw(values, region, current, result["recommended"], low, high, args.output_dir)
    summary = {"method": {"primary_rule": "Freedman-Diaconis",
                           "population_constraint": (f"soft cap targeting >={args.minimum_events_per_bin} "
                                                     "average events/bin, subject to minimum-bins floor"),
                           "note": "Recommendation uses selected data only; validate against fit stability."},
               "input": {"directory": str(args.data_dir), "files": len(files),
                         "available_entries": total, "processed_entries": processed, "cutflow": cutflow},
               "configuration": {"minimum_events_per_bin": args.minimum_events_per_bin,
                                 "minimum_bins": args.minimum_bins, "maximum_bins": args.maximum_bins},
               "regions": results}
    (args.output_dir / "binning_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("REGION EVENTS CURRENT FD SCOTT STURGES CAP RECOMMENDED")
    for region, value in results.items():
        print(f"{scan.abbreviated_region(region):5s} {value['events']:6d} {value['current_bins']:7d} "
              f"{value['freedman_diaconis']:2d} {value['scott']:5d} {value['sturges']:7d} "
              f"{value['population_cap']:3d} {value['recommended']:11d}")
    print(f"Wrote binning study to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
