#!/usr/bin/env python3
"""First-pass Z->mumu data/MC comparison using the pfExtractor ROOT trees."""

from __future__ import annotations

import argparse
import json
import math
import sys
from array import array
from pathlib import Path

try:
    import ROOT
except ImportError as exc:
    raise SystemExit(
        "PyROOT is required. On the workshop machine run with:\n"
        "  conda run -n hep-play python analyze_zmumu.py [options]"
    ) from exc

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0)

DEFAULT_DATA = Path("../../HackathonDataRaw/2016H")
DEFAULT_MC = Path("../../HackathonDataRaw/ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8")
Z_MASS = 91.1876


def arguments() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=(here / DEFAULT_DATA).resolve())
    parser.add_argument("--mc-dir", type=Path, default=(here / DEFAULT_MC).resolve())
    parser.add_argument("--output-dir", type=Path, default=here / "plots")
    parser.add_argument("--max-events", type=int, default=-1,
                        help="Maximum entries per sample; -1 processes all entries")
    parser.add_argument("--max-files", type=int, default=-1,
                        help="Maximum files per sample; -1 uses all files")
    parser.add_argument("--label", default="CMS Open Data 2016H (subset)")
    return parser.parse_args()


def root_files(directory: Path, maximum: int) -> list[Path]:
    files = sorted(directory.glob("*.root"))
    if maximum >= 0:
        files = files[:maximum]
    if not files:
        raise FileNotFoundError(f"No ROOT files found in {directory}")
    usable = []
    for filename in files:
        source = ROOT.TFile.Open(str(filename), "READ")
        tree = source.Get("pfExtractor/pfTree") if source and not source.IsZombie() else None
        if tree and tree.InheritsFrom("TTree"):
            usable.append(filename)
        else:
            print(f"WARNING: skipping file without pfExtractor/pfTree: {filename}", file=sys.stderr)
        if source:
            source.Close()
    if not usable:
        raise RuntimeError(f"None of the ROOT files in {directory} contains pfExtractor/pfTree")
    return usable


def chain(files: list[Path]) -> ROOT.TChain:
    result = ROOT.TChain("pfExtractor/pfTree")
    for filename in files:
        if result.Add(str(filename)) == 0:
            raise OSError(f"Could not add {filename}")
    # These ntuples contain large PF-candidate and jet collections. Reading only
    # the fields used below makes full-sample running substantially faster and
    # avoids retaining irrelevant collection buffers.
    result.SetBranchStatus("*", 0)
    for branch in (
        "nMuons", "nVertices", "IsoMu24", "muonPt", "muonEta", "muonPhi",
        "muonMass", "muonCharge", "muonDxy", "muonDz", "muonIsoCharged",
        "muonIsoNeutral", "muonIsoPhoton", "muonIsoPU", "muonIsTight", "muonIsPF"
    ):
        result.SetBranchStatus(branch, 1)
    return result


def histograms(prefix: str) -> dict[str, ROOT.TH1D]:
    specs = {
        "mass": (60, 60.0, 120.0, "m_{#mu#mu} [GeV]"),
        "lead_pt": (50, 20.0, 120.0, "leading muon p_{T} [GeV]"),
        "sublead_pt": (45, 10.0, 100.0, "subleading muon p_{T} [GeV]"),
        "muon_eta": (48, -2.4, 2.4, "muon #eta"),
        "z_pt": (50, 0.0, 100.0, "Z candidate p_{T} [GeV]"),
        "nvertices": (50, 0.0, 50.0, "number of vertices"),
    }
    return {
        name: ROOT.TH1D(f"{prefix}_{name}", f";{axis};Events", bins, low, high)
        for name, (bins, low, high, axis) in specs.items()
    }


def relative_isolation(tree, index: int) -> float:
    neutral = tree.muonIsoNeutral[index] + tree.muonIsoPhoton[index]
    corrected = tree.muonIsoCharged[index] + max(0.0, neutral - 0.5 * tree.muonIsoPU[index])
    return corrected / tree.muonPt[index] if tree.muonPt[index] > 0 else math.inf


def selected_muons(tree) -> list[int]:
    selected = []
    for i in range(tree.nMuons):
        if (
            tree.muonPt[i] > 15.0
            and abs(tree.muonEta[i]) < 2.4
            and bool(tree.muonIsTight[i])
            and bool(tree.muonIsPF[i])
            and abs(tree.muonDxy[i]) < 0.05
            and abs(tree.muonDz[i]) < 0.10
            and relative_isolation(tree, i) < 0.15
        ):
            selected.append(i)
    return selected


def four_vector(tree, i: int) -> ROOT.Math.PtEtaPhiMVector:
    return ROOT.Math.PtEtaPhiMVector(
        tree.muonPt[i], tree.muonEta[i], tree.muonPhi[i], tree.muonMass[i]
    )


def best_pair(tree, indices: list[int]):
    candidates = []
    for position, first in enumerate(indices):
        for second in indices[position + 1:]:
            if tree.muonCharge[first] * tree.muonCharge[second] >= 0:
                continue
            p4 = four_vector(tree, first) + four_vector(tree, second)
            candidates.append((abs(p4.M() - Z_MASS), first, second, p4))
    return min(candidates, default=None, key=lambda item: item[0])


def process(files: list[Path], sample: str, max_events: int):
    tree = chain(files)
    hist = histograms(sample)
    cutflow = {name: 0 for name in (
        "all", "trigger", ">=2 selected muons", "opposite-sign pair",
        "leading pT > 26 GeV", "60 < mass < 120 GeV"
    )}
    entries = tree.GetEntries()
    limit = entries if max_events < 0 else min(entries, max_events)

    for entry in range(limit):
        tree.GetEntry(entry)
        cutflow["all"] += 1
        if not bool(tree.IsoMu24):
            continue
        cutflow["trigger"] += 1
        muons = selected_muons(tree)
        if len(muons) < 2:
            continue
        cutflow[">=2 selected muons"] += 1
        pair = best_pair(tree, muons)
        if pair is None:
            continue
        cutflow["opposite-sign pair"] += 1
        _, first, second, z = pair
        if tree.muonPt[first] < tree.muonPt[second]:
            first, second = second, first
        if tree.muonPt[first] <= 26.0:
            continue
        cutflow["leading pT > 26 GeV"] += 1
        if not 60.0 < z.M() < 120.0:
            continue
        cutflow["60 < mass < 120 GeV"] += 1
        hist["mass"].Fill(z.M())
        hist["lead_pt"].Fill(tree.muonPt[first])
        hist["sublead_pt"].Fill(tree.muonPt[second])
        hist["muon_eta"].Fill(tree.muonEta[first])
        hist["muon_eta"].Fill(tree.muonEta[second])
        hist["z_pt"].Fill(z.Pt())
        hist["nvertices"].Fill(tree.nVertices)

    return hist, cutflow, entries, limit


def draw_comparison(data, mc, name: str, output: Path, label: str) -> None:
    canvas = ROOT.TCanvas(f"c_{name}", "", 800, 800)
    top = ROOT.TPad("top", "", 0, 0.29, 1, 1)
    bottom = ROOT.TPad("bottom", "", 0, 0, 1, 0.29)
    top.SetBottomMargin(0.02); bottom.SetTopMargin(0.03); bottom.SetBottomMargin(0.34)
    top.Draw(); bottom.Draw(); top.cd()

    d = data.Clone(f"draw_data_{name}")
    m = mc.Clone(f"draw_mc_{name}")
    if m.Integral() > 0:
        m.Scale(d.Integral() / m.Integral())
    d.SetMarkerStyle(20); d.SetMarkerSize(0.8); d.SetLineColor(ROOT.kBlack)
    m.SetLineColor(ROOT.kAzure + 2); m.SetLineWidth(2); m.SetFillColorAlpha(ROOT.kAzure - 9, 0.45)
    maximum = max(d.GetMaximum(), m.GetMaximum()) * 1.35
    m.SetMaximum(maximum); m.SetMinimum(0); m.GetXaxis().SetLabelSize(0)
    m.GetYaxis().SetTitle("Events (MC normalized to data)")
    m.Draw("HIST"); d.Draw("E1 SAME")
    legend = ROOT.TLegend(0.62, 0.70, 0.88, 0.87)
    legend.SetBorderSize(0); legend.AddEntry(d, "Data", "lep"); legend.AddEntry(m, "Z#rightarrow#mu#mu simulation", "lf")
    legend.Draw()
    text = ROOT.TLatex(); text.SetNDC(); text.SetTextSize(0.038); text.DrawLatex(0.14, 0.92, label)

    bottom.cd()
    ratio = d.Clone(f"ratio_{name}"); ratio.Divide(m)
    ratio.SetTitle(""); ratio.GetYaxis().SetTitle("Data / MC"); ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetRangeUser(0.5, 1.5); ratio.GetYaxis().SetTitleSize(0.105)
    ratio.GetYaxis().SetTitleOffset(0.45); ratio.GetYaxis().SetLabelSize(0.09)
    ratio.GetXaxis().SetTitleSize(0.12); ratio.GetXaxis().SetLabelSize(0.10)
    ratio.Draw("E1")
    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2); line.Draw()
    canvas.SaveAs(str(output / f"{name}.png"))
    canvas.SaveAs(str(output / f"{name}.pdf"))


def main() -> int:
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data_files = root_files(args.data_dir, args.max_files)
    mc_files = root_files(args.mc_dir, args.max_files)
    data_hist, data_flow, data_entries, data_read = process(data_files, "data", args.max_events)
    mc_hist, mc_flow, mc_entries, mc_read = process(mc_files, "mc", args.max_events)

    if data_flow["60 < mass < 120 GeV"] == 0 or mc_flow["60 < mass < 120 GeV"] == 0:
        raise RuntimeError("Selection produced an empty final sample")
    for name in data_hist:
        draw_comparison(data_hist[name], mc_hist[name], name, args.output_dir, args.label)

    root_output = ROOT.TFile(str(args.output_dir / "histograms.root"), "RECREATE")
    for collection in (data_hist, mc_hist):
        for histogram in collection.values():
            histogram.Write()
    root_output.Close()
    summary = {
        "configuration": {
            "data_directory": str(args.data_dir), "mc_directory": str(args.mc_dir),
            "max_events_per_sample": args.max_events, "max_files_per_sample": args.max_files,
            "tree": "pfExtractor/pfTree", "mc_normalization": "unit/shape normalized to selected data yield"
        },
        "data": {"files": len(data_files), "available_entries": data_entries,
                 "processed_entries": data_read, "cutflow": data_flow},
        "mc": {"files": len(mc_files), "available_entries": mc_entries,
               "processed_entries": mc_read, "cutflow": mc_flow},
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Wrote plots and results to {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
