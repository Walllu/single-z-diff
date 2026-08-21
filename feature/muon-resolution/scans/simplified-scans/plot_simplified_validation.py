#!/usr/bin/env python3
"""Regenerate BB/BE/EE best-fit mass ratio and profile-pull plots."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import simplified_scan as scan


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="simplified scan summary.json")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--regions", nargs="+", choices=scan.REGIONS, default=list(scan.REGIONS))
    parser.add_argument("--landscape-only", choices=("fine_scale", "fine_resolution"), default=None)
    return parser.parse_args()


def template(histogram) -> scan.hf.Template:
    bins = histogram.GetNbinsX()
    edges = np.asarray([histogram.GetXaxis().GetBinLowEdge(i) for i in range(1, bins + 2)])
    content = np.asarray([histogram.GetBinContent(i) for i in range(1, bins + 1)])
    variance = np.asarray([histogram.GetBinError(i) ** 2 for i in range(1, bins + 1)])
    return scan.hf.Template(edges, content, variance, float(content.sum()))


def main() -> int:
    args = arguments()
    summary = json.loads(args.summary.read_text())
    source = args.summary.parent / "simplified_histograms.root"
    output = args.output_dir or args.summary.parent
    root_file = scan.ROOT.TFile.Open(str(source), "READ")
    if not root_file or root_file.IsZombie():
        raise SystemExit(f"Could not open {source}")
    best = summary["best_fit"]
    if args.landscape_only:
        root_file.Close()
        scan.draw_profile_planes(summary["scan_records"], best, output, (args.landscape_only,))
        print(f"Wrote {args.landscape_only} landscape under {output}")
        return 0
    for region in args.regions:
        objects = [root_file.Get(f"{name}_{region}_mass") for name in
                   ("data", "nominal_mc", "corrected_mc")]
        if any(not obj for obj in objects):
            raise SystemExit(f"Missing {region} histograms in {source}")
        scan.draw_validation(region, *(template(obj) for obj in objects), best, output / region)
    root_file.Close()
    if "scan_records" in summary and set(args.regions) == set(scan.REGIONS):
        scan.draw_profile_planes(summary["scan_records"], best, output)
    print(f"Wrote regenerated validation plots under {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
