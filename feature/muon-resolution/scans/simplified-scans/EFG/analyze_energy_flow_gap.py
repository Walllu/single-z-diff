#!/usr/bin/env python3
"""Compute and plot Energy Gap Flow observables in selected BB/BE/EE events."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SIMPLIFIED = HERE.parent
PILEUP = SIMPLIFIED / "pileup-reweighting"
sys.path[:0] = [str(SIMPLIFIED), str(PILEUP)]

import simplified_scan as scan  # noqa: E402
import derive_and_validate_pileup as pu  # noqa: E402
from pileup_weights import PileupCalibration  # noqa: E402
from energy_flow_gap import event_observables  # noqa: E402


EFG_CONFIGURATION = {
    # These candidate-level choices follow Analysis.ipynb in this directory.
    "eta_min": -2.5, "eta_max": 2.5, "pt_min_gev": 0.5,
    "beta": 0.1, "q0_gev": 2.0,
    "candidate_from_pv": 3,
    "muon_veto_delta_r": 0.01,
}

OBSERVABLES = {
    "efg_score": ([50, 0.0, 1.0], "max(G_{-}, G_{+})"),
    "efg_score_inclusive": ([50, 0.0, 1.0], "inclusive max(G_{-}, G_{+})"),
    "efg_effective_gap": ([50, 0.0, 5.0], "#Delta#eta^{EFG}"),
    "forward_gap": ([50, 0.0, 5.0], "traditional #Delta#eta^{F}"),
    "efg_forward": ([50, 0.0, 1.0], "G_{+}"),
    "efg_backward": ([50, 0.0, 1.0], "G_{-}"),
    "charged_multiplicity": ([60, 0.0, 300.0], "fiducial charged multiplicity"),
    "charged_sum_pt": ([60, 0.0, 300.0], "fiducial charged #Sigma p_{T} [GeV]"),
}

ROC_DIRECTION = {
    "efg_score": "high", "efg_score_inclusive": "high",
    "efg_effective_gap": "high", "forward_gap": "high",
    "efg_forward": "high", "efg_backward": "high",
    "charged_multiplicity": "low", "charged_sum_pt": "low",
}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    raw = scan.baseline.default_raw_data_directory()
    parser.add_argument("--data-dir", type=Path, default=raw / "2016H")
    parser.add_argument("--mc-dir", type=Path,
                        default=raw / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8")
    parser.add_argument("--signal-dir", type=Path,
                        default=raw / "DYToMuMu_pomflux_Pt-30_TuneCP5_13TeV-pythia8")
    parser.add_argument("--skip-signal", action="store_true",
                        help="Do not process the optional diffractive signal sample")
    parser.add_argument("--pileup-weights", type=Path,
                        default=SIMPLIFIED.parents[1] / "calibrations/experimental/pileup_nvertices_2016H_bb_be_ee.json")
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--label", default="full_production")
    parser.add_argument("--max-files", type=int, default=-1)
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("--regions", nargs="+", choices=scan.REGIONS, default=list(scan.REGIONS))
    parser.add_argument("--observables", nargs="+", choices=tuple(OBSERVABLES), default=list(OBSERVABLES))
    parser.add_argument("--render-only", action="store_true")
    return parser.parse_args()


def remove_z_muon_overlaps(eta, phi, muons):
    """Notebook-compatible veto of every candidate within dR of either Z muon."""
    keep = np.ones(len(eta), dtype=bool)
    for muon in muons:
        deta = eta - muon[1]
        dphi = np.arctan2(np.sin(phi - muon[2]), np.cos(phi - muon[2]))
        keep &= np.hypot(deta, dphi) >= EFG_CONFIGURATION["muon_veto_delta_r"]
    return keep


def read_events(files, maximum_events, regions):
    tree = scan.baseline.make_chain(files)
    for branch in ("candPt", "candEta", "candPhi", "candFromPV"):
        tree.SetBranchStatus(branch, 1)
    total = tree.GetEntries(); limit = total if maximum_events < 0 else min(total, maximum_events)
    result = {region: [] for region in scan.REGIONS}
    flow = {key: 0 for key in ("all", "trigger", "two_muons", "opposite_sign",
                                "leading_pt", "fit_mass", "eta_region", "selected")}
    removed_counts = {"zero": 0, "one": 0, "two_or_more": 0}
    empty_candidate_events = 0
    for entry in range(limit):
        if entry and entry % 500_000 == 0:
            print(f"  processed {entry:,}/{limit:,} entries; selected {flow['selected']:,}",
                  flush=True)
        tree.GetEntry(entry); flow["all"] += 1
        if scan.baseline.SELECTION["require_trigger"] and not bool(getattr(tree, scan.baseline.SELECTION["trigger_branch"])):
            continue
        flow["trigger"] += 1
        selected = scan.baseline.selected_muons(tree)
        if len(selected) < 2: continue
        flow["two_muons"] += 1
        pair = scan.baseline.best_pair(tree, selected)
        if pair is None: continue
        flow["opposite_sign"] += 1
        _, lead, sublead, mass = pair
        if tree.muonPt[lead] < tree.muonPt[sublead]: lead, sublead = sublead, lead
        if tree.muonPt[lead] <= scan.baseline.SELECTION["leading_min_pt_gev"]: continue
        flow["leading_pt"] += 1
        if not 70.0 < mass < 110.0: continue
        flow["fit_mass"] += 1
        lead_muon = (float(tree.muonPt[lead]), float(tree.muonEta[lead]), float(tree.muonPhi[lead]), float(tree.muonMass[lead]))
        sub_muon = (float(tree.muonPt[sublead]), float(tree.muonEta[sublead]), float(tree.muonPhi[sublead]), float(tree.muonMass[sublead]))
        event = {"lead": lead_muon, "sublead": sub_muon}
        region = scan.simple_region(event)
        flow["eta_region"] += 1
        if region not in regions: continue
        eta = np.asarray(tree.candEta, dtype=np.float64); phi = np.asarray(tree.candPhi, dtype=np.float64)
        pt = np.asarray(tree.candPt, dtype=np.float64)
        from_pv = np.asarray(tree.candFromPV, dtype=np.int32)
        base = ((eta > EFG_CONFIGURATION["eta_min"])
                & (eta < EFG_CONFIGURATION["eta_max"])
                & (pt > EFG_CONFIGURATION["pt_min_gev"])
                & (from_pv == EFG_CONFIGURATION["candidate_from_pv"]))
        inclusive = event_observables(eta[base], pt[base], eta_min=EFG_CONFIGURATION["eta_min"],
            eta_max=EFG_CONFIGURATION["eta_max"], pt_min=EFG_CONFIGURATION["pt_min_gev"],
            beta=EFG_CONFIGURATION["beta"], q0=EFG_CONFIGURATION["q0_gev"])
        nonoverlap = remove_z_muon_overlaps(eta, phi, (lead_muon, sub_muon))
        clean = base & nonoverlap
        removed = int(np.count_nonzero(base & ~nonoverlap))
        exclusive = event_observables(eta[clean], pt[clean], eta_min=EFG_CONFIGURATION["eta_min"],
            eta_max=EFG_CONFIGURATION["eta_max"], pt_min=EFG_CONFIGURATION["pt_min_gev"],
            beta=EFG_CONFIGURATION["beta"], q0=EFG_CONFIGURATION["q0_gev"])
        if not np.any(clean): empty_candidate_events += 1
        result[region].append({
            "efg_score": exclusive["efg_score"], "efg_score_inclusive": inclusive["efg_score"],
            "efg_effective_gap": exclusive["efg_effective_gap"], "forward_gap": exclusive["forward_gap"],
            "efg_forward": exclusive["efg_forward"], "efg_backward": exclusive["efg_backward"],
            "charged_multiplicity": float(np.count_nonzero(clean)),
            "charged_sum_pt": float(pt[clean].sum()), "nvertices": int(tree.nVertices),
        })
        removed_counts[("zero", "one", "two_or_more")[min(removed, 2)]] += 1; flow["selected"] += 1
    return result, flow, removed_counts, empty_candidate_events, total, limit


def template(values, edges, weights):
    content, _ = np.histogram(values, bins=edges, weights=weights)
    variance, _ = np.histogram(values, bins=edges, weights=np.square(weights))
    return scan.hf.Template(edges, content.astype(float), variance.astype(float), float(content.sum()))


def from_root(histogram):
    bins = histogram.GetNbinsX(); edges = np.asarray([histogram.GetXaxis().GetBinLowEdge(i) for i in range(1, bins + 2)])
    content = np.asarray([histogram.GetBinContent(i) for i in range(1, bins + 1)])
    variance = np.asarray([histogram.GetBinError(i)**2 for i in range(1, bins + 1)])
    return scan.hf.Template(edges, content, variance, float(content.sum()))


def draw_signal_diagnostics(region, observable, background, signal, output):
    """Draw the notebook-inspired diffractive/inclusive shape and ROC checks."""
    output.mkdir(parents=True, exist_ok=True)
    bsum, ssum = background.content.sum(), signal.content.sum()
    if bsum <= 0 or ssum <= 0:
        return
    background = pu.normalized(background, 1.0 / bsum)
    signal = pu.normalized(signal, 1.0 / ssum)
    axis_title = OBSERVABLES[observable][1]
    token = f"signal_{region}_{observable}"
    b = scan.hf.root_histogram(f"background_{token}", background, axis_title)
    s = scan.hf.root_histogram(f"signal_{token}", signal, axis_title)
    canvas = scan.ROOT.TCanvas(f"canvas_{token}", "", 850, 700)
    b.SetLineColor(scan.ROOT.kBlue + 1); b.SetLineWidth(3)
    s.SetLineColor(scan.ROOT.kRed + 1); s.SetLineWidth(3)
    b.SetMaximum(1.25 * max(b.GetMaximum(), s.GetMaximum(), 1e-9))
    b.GetYaxis().SetTitle("Normalized events")
    b.Draw("HIST"); s.Draw("HIST SAME")
    legend = scan.ROOT.TLegend(.56, .75, .91, .89)
    legend.SetBorderSize(0); legend.SetFillStyle(0)
    legend.AddEntry(b, "inclusive DY simulation", "l")
    legend.AddEntry(s, "diffractive pomflux simulation", "l"); legend.Draw()
    title = scan.ROOT.TLatex(); title.SetNDC(); title.SetTextSize(.037)
    title.DrawLatex(.13, .92, f"{region}: diffractive EFG diagnostic")
    scan.save_canvas(canvas, output / f"{region}_{observable}_signal_shape",
                     prefer_pad_raster=True)

    direction = ROC_DIRECTION[observable]
    if direction == "high":
        signal_efficiency = np.r_[np.cumsum(signal.content[::-1])[::-1], 0.0]
        background_efficiency = np.r_[np.cumsum(background.content[::-1])[::-1], 0.0]
    else:
        signal_efficiency = np.r_[0.0, np.cumsum(signal.content)]
        background_efficiency = np.r_[0.0, np.cumsum(background.content)]
    roc = scan.ROOT.TGraph(len(signal_efficiency), background_efficiency, signal_efficiency)
    roc.SetLineColor(scan.ROOT.kRed + 1); roc.SetLineWidth(3)
    roc.SetTitle(";inclusive DY efficiency;diffractive signal efficiency")
    roc_canvas = scan.ROOT.TCanvas(f"roc_canvas_{token}", "", 750, 700)
    roc.Draw("AL"); roc.GetXaxis().SetLimits(0, 1); roc.SetMinimum(0); roc.SetMaximum(1)
    diagonal = scan.ROOT.TLine(0, 0, 1, 1); diagonal.SetLineStyle(2); diagonal.Draw()
    title.DrawLatex(.13, .92, f"{region}: {direction}-{observable} ROC")
    scan.save_canvas(roc_canvas, output / f"{region}_{observable}_signal_roc",
                     prefer_pad_raster=True)


def render(args, output):
    source = scan.ROOT.TFile.Open(str(output / "efg_histograms.root"), "READ")
    for region in args.regions:
        for observable in args.observables:
            objects = [source.Get(f"{kind}_{region}_{observable}") for kind in ("data", "before", "after")]
            if any(not item for item in objects): raise SystemExit(f"Missing {region}/{observable} EFG histograms")
            pu.draw_comparison(region, observable, *(from_root(item) for item in objects), output / region,
                title_text=f"{region}: Energy Flow Gap validation",
                before_label="MC without vertex weight (dashed)",
                after_label="MC + vertex weight (solid)", filename_suffix="efg",
                axis_title=OBSERVABLES[observable][1])
            signal = source.Get(f"signal_{region}_{observable}")
            if signal:
                draw_signal_diagnostics(region, observable, from_root(objects[1]),
                                        from_root(signal), output / region)
    source.Close(); return 0


def main():
    start = time.perf_counter(); args = arguments(); output = args.output_dir / args.label
    if args.render_only: return render(args, output)
    data_files = scan.baseline.discover(args.data_dir, args.max_files)
    mc_files = scan.baseline.discover(args.mc_dir, args.max_files)
    signal_files = [] if args.skip_signal else scan.baseline.discover(args.signal_dir, args.max_files)
    data, data_flow, data_removed, data_empty, data_total, data_processed = read_events(data_files, args.max_events, set(args.regions))
    mc, mc_flow, mc_removed, mc_empty, mc_total, mc_processed = read_events(mc_files, args.max_events, set(args.regions))
    if signal_files:
        signal, signal_flow, signal_removed, signal_empty, signal_total, signal_processed = read_events(
            signal_files, args.max_events, set(args.regions))
    else:
        signal = {region: [] for region in scan.REGIONS}
        signal_flow = signal_removed = {}; signal_empty = signal_total = signal_processed = 0
    pileup = PileupCalibration.from_json(args.pileup_weights)
    output.mkdir(parents=True, exist_ok=True)
    summary = {"configuration": {"energy_flow_gap": EFG_CONFIGURATION,
                "z_muon_treatment": "geometrically matched candidates removed for primary observables",
                "reference_gist": str((HERE / "gistfile1.txt").resolve()),
                "reference_notebook": str((HERE / "Analysis.ipynb").resolve()),
                "z_selection_note": "Established calibration selection retained; notebook candidate cleaning adopted",
                "pileup_weights": str(args.pileup_weights.resolve())},
               "inputs": {"data_files": len(data_files), "mc_files": len(mc_files), "signal_files": len(signal_files),
                 "data_entries": data_total, "mc_entries": mc_total, "signal_entries": signal_total,
                 "data_processed": data_processed, "mc_processed": mc_processed,
                 "signal_processed": signal_processed},
               "cutflows": {"data": data_flow, "mc": mc_flow, "signal": signal_flow},
               "z_muon_candidate_overlaps_removed": {"data": data_removed, "mc": mc_removed, "signal": signal_removed},
               "selected_events_without_clean_candidates": {"data": data_empty, "mc": mc_empty, "signal": signal_empty},
               "regions": {}}
    root_file = scan.ROOT.TFile(str(output / "efg_histograms.root"), "RECREATE")
    for region in args.regions:
        data_weights = np.ones(len(data[region])); before_weights = np.ones(len(mc[region]))
        after_weights = np.asarray([pileup.weight(region, event["nvertices"]) for event in mc[region]])
        summary["regions"][region] = {"counts": {"data": len(data[region]), "mc": len(mc[region]),
                                                   "signal": len(signal[region])}, "observables": {}}
        for observable in args.observables:
            bins, low, high = OBSERVABLES[observable][0]; edges = np.linspace(low, high, bins + 1)
            data_hist = template([x[observable] for x in data[region]], edges, data_weights)
            before_hist = template([x[observable] for x in mc[region]], edges, before_weights)
            after_hist = template([x[observable] for x in mc[region]], edges, after_weights)
            signal_hist = template([x[observable] for x in signal[region]], edges,
                                   np.ones(len(signal[region])))
            before_metric = scan.hf.barlow_beeston(data_hist, before_hist); after_metric = scan.hf.barlow_beeston(data_hist, after_hist)
            summary["regions"][region]["observables"][observable] = {
                "before": scan.hf.slim_metrics(before_metric), "after": scan.hf.slim_metrics(after_metric),
                "delta_vs_before": before_metric["objective"] - after_metric["objective"]}
            root_file.cd()
            for kind, hist in (("data", data_hist), ("before", before_hist), ("after", after_hist)):
                scan.hf.root_histogram(f"{kind}_{region}_{observable}", hist, OBSERVABLES[observable][1]).Write()
            if signal_files:
                scan.hf.root_histogram(f"signal_{region}_{observable}", signal_hist,
                                       OBSERVABLES[observable][1]).Write()
    root_file.Close(); summary["timing_seconds"] = {"extraction_histogramming": time.perf_counter() - start}
    (output / "efg_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for region in args.regions:
        for observable in args.observables:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), "--output-dir", str(args.output_dir),
                "--label", args.label, "--render-only", "--regions", region, "--observables", observable], check=True)
    print(json.dumps({"output": str(output), "regions": args.regions, "observables": args.observables,
                      "timing_seconds": time.perf_counter() - start}, indent=2)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
