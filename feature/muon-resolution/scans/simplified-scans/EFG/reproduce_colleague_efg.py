#!/usr/bin/env python3
"""Reproduce the colleague EFG notebook selection and plots with robust PyROOT code.

Muon scale/resolution corrections are optional.  The notebook's selection and
candidate definition are kept separate from the calibration selection used by
analyze_energy_flow_gap.py.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from array import array
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SIMPLIFIED = HERE.parent
sys.path.insert(0, str(SIMPLIFIED))

import simplified_scan as scan  # noqa: E402
from energy_flow_gap import event_observables  # noqa: E402
from muon_corrections import Calibration, correct_muon  # noqa: E402


SELECTION = {
    "trigger": "IsoMu24",
    "muon_min_pt_gev": 22.0,
    "require_muon_is_loose": True,
    "number_of_muons": 2,
    "dimuon_mass_gev": [80.0, 100.0],
    "dimuon_max_pt_gev": 15.0,
    "require_opposite_sign": False,
    "use_massless_muons": True,
    "apply_muon_scale_or_smearing": False,
}

CANDIDATES = {
    "min_pt_gev": 0.5,
    "from_pv": 3,
    "muon_overlap_delta_r": 0.01,
    "efg_eta_range": [-2.5, 2.5],
    "efg_beta": 0.1,
    "efg_q0_gev": 2.0,
}

SPECS = {
    # The edges reproduce the notebook's np.linspace calls exactly.
    "track_pt": (np.linspace(0.0, 2.0, 30), "cleaned charged-candidate p_{T} [GeV]", False, "Candidates"),
    "track_eta": (np.linspace(-4.0, 4.0, 30), "cleaned charged-candidate #eta", False, "Candidates"),
    "rapidity_gap": (np.linspace(0.0, 5.0, 20), "maximum forward rapidity gap #Delta#eta", True, "Events"),
    "efg_score": (np.linspace(0.0, 1.0, 20), "Energy Gap Flow max(G_{-}, G_{+})", True, "Events"),
}


def arguments():
    raw = scan.baseline.default_raw_data_directory()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=raw / "2016H")
    parser.add_argument("--background-dir", type=Path,
                        default=raw / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8")
    parser.add_argument("--signal-dir", type=Path,
                        default=raw / "DYToMuMu_pomflux_Pt-30_TuneCP5_13TeV-pythia8")
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots/colleague-reproduction")
    parser.add_argument("--label", default="full_samples_no_muon_corrections")
    parser.add_argument("--max-events", type=int, default=-1,
                        help="Maximum entries per sample; -1 means all")
    parser.add_argument("--notebook-file-subset", action="store_true",
                        help="Use the notebook's 2 data, 3 DY, and 3 signal file counts")
    parser.add_argument("--apply-muon-corrections", action="store_true",
                        help="Apply fitted barrel/endcap scale and resolution to both MC samples")
    parser.add_argument("--calibration", type=Path,
                        default=SIMPLIFIED.parents[1] / "calibrations/muon_momentum_2016H_bb_be_ee.json")
    parser.add_argument("--correction-seeds", nargs="+", type=int,
                        default=[314159, 271828, 161803])
    return parser.parse_args()


def histogram(values, edges):
    content, _ = np.histogram(values, bins=edges)
    return content.astype(float)


def clean_candidates(eta, phi, pt, from_pv, muons):
    keep = (pt > CANDIDATES["min_pt_gev"]) & (from_pv == CANDIDATES["from_pv"])
    for mu_eta, mu_phi in muons:
        dphi = np.arctan2(np.sin(phi - mu_phi), np.cos(phi - mu_phi))
        keep &= np.hypot(eta - mu_eta, dphi) >= CANDIDATES["muon_overlap_delta_r"]
    return keep


def conventional_gap(eta, pt):
    accepted = ((np.abs(eta) < CANDIDATES["efg_eta_range"][1])
                & (pt > CANDIDATES["min_pt_gev"]))
    if not np.any(accepted):
        # This is the notebook convention, intentionally distinct from the gist.
        return 0.0
    selected_eta = eta[accepted]
    return float(max(CANDIDATES["efg_eta_range"][1] - selected_eta.max(),
                     CANDIDATES["efg_eta_range"][1] + selected_eta.min()))


def process(files, maximum_events, sample, calibration, correction_seeds):
    tree = scan.baseline.make_chain(files)
    for branch in ("run", "lumi", "event", "muonIsLoose", "candPt", "candEta",
                   "candPhi", "candFromPV"):
        tree.SetBranchStatus(branch, 1)
    entries = int(tree.GetEntries())
    limit = entries if maximum_events < 0 else min(entries, maximum_events)
    counts = {name: np.zeros(len(spec[0]) - 1) for name, spec in SPECS.items()}
    event_scores = {"rapidity_gap": [], "efg_score": []}
    flow = {key: 0.0 for key in ("all", "trigger", "two_loose_muons", "mass_window",
                                  "z_pt", "selected")}
    ordering_violations = 0.0
    cleaned_candidates = 0.0
    out_of_calibration = 0
    apply_correction = calibration is not None and sample != "data"
    seeds = correction_seeds if apply_correction else [None]
    replica_weight = 1.0 / len(seeds)
    for entry in range(limit):
        if entry and entry % 500_000 == 0:
            print(f"{sample}: {entry:,}/{limit:,}, selected {flow['selected']:,}", flush=True)
        tree.GetEntry(entry); flow["all"] += 1
        if not bool(getattr(tree, SELECTION["trigger"])):
            continue
        flow["trigger"] += 1
        event_key = f"{int(tree.run)}:{int(tree.lumi)}:{int(tree.event)}"
        passing = []
        for seed in seeds:
            muons = []
            for i in range(int(tree.nMuons)):
                if not bool(tree.muonIsLoose[i]):
                    continue
                raw = (float(tree.muonPt[i]), float(tree.muonEta[i]),
                       float(tree.muonPhi[i]), float(tree.muonMass[i]))
                if apply_correction and abs(raw[1]) < calibration.muon_max_abs_eta:
                    corrected = correct_muon(raw, calibration, seed=seed,
                                             event_key=event_key, muon_key=i)
                    values = (corrected.pt, corrected.eta, corrected.phi)
                else:
                    values = raw[:3]
                    if apply_correction and abs(raw[1]) >= calibration.muon_max_abs_eta:
                        out_of_calibration += 1
                if values[0] > SELECTION["muon_min_pt_gev"]:
                    muons.append((i, *values))
            if len(muons) < 2:
                continue
            flow["two_loose_muons"] += replica_weight
            first, second = muons[:2]
            _, pt1, eta1, phi1 = first; _, pt2, eta2, phi2 = second
            if pt1 < pt2:
                ordering_violations += replica_weight
            dphi = math.atan2(math.sin(phi1 - phi2), math.cos(phi1 - phi2))
            mass2 = 2.0 * pt1 * pt2 * (math.cosh(eta1 - eta2) - math.cos(dphi))
            mass = math.sqrt(max(mass2, 0.0))
            if not SELECTION["dimuon_mass_gev"][0] < mass < SELECTION["dimuon_mass_gev"][1]:
                continue
            flow["mass_window"] += replica_weight
            z_px = pt1 * math.cos(phi1) + pt2 * math.cos(phi2)
            z_py = pt1 * math.sin(phi1) + pt2 * math.sin(phi2)
            if math.hypot(z_px, z_py) >= SELECTION["dimuon_max_pt_gev"]:
                continue
            flow["z_pt"] += replica_weight
            passing.append((first[0], second[0], eta1, phi1, eta2, phi2))
        if not passing:
            continue
        eta = np.asarray(tree.candEta, dtype=np.float64)
        phi = np.asarray(tree.candPhi, dtype=np.float64)
        pt = np.asarray(tree.candPt, dtype=np.float64)
        from_pv = np.asarray(tree.candFromPV, dtype=np.int32)
        groups = {}
        for first, second, eta1, phi1, eta2, phi2 in passing:
            groups.setdefault((first, second, eta1, phi1, eta2, phi2), 0)
            groups[(first, second, eta1, phi1, eta2, phi2)] += 1
        for (_, _, eta1, phi1, eta2, phi2), replicas in groups.items():
            weight = replicas * replica_weight
            keep = clean_candidates(eta, phi, pt, from_pv, ((eta1, phi1), (eta2, phi2)))
            clean_eta, clean_pt = eta[keep], pt[keep]
            cleaned_candidates += weight * len(clean_pt)
            counts["track_pt"] += weight * histogram(clean_pt, SPECS["track_pt"][0])
            counts["track_eta"] += weight * histogram(clean_eta, SPECS["track_eta"][0])
            gap = conventional_gap(clean_eta, clean_pt)
            efg = event_observables(clean_eta, clean_pt, eta_min=-2.5, eta_max=2.5, pt_min=0.5,
                                    beta=0.1, q0=2.0)["efg_score"]
            counts["rapidity_gap"] += weight * histogram([gap], SPECS["rapidity_gap"][0])
            counts["efg_score"] += weight * histogram([efg], SPECS["efg_score"][0])
            event_scores["rapidity_gap"].extend([gap] * replicas)
            event_scores["efg_score"].extend([efg] * replicas)
            flow["selected"] += weight
    return {"histograms": counts, "scores": event_scores, "cutflow": flow,
            "entries": entries, "processed": limit,
            "cleaned_candidates": cleaned_candidates,
            "first_two_not_pt_ordered": ordering_violations,
            "muon_replicas_outside_calibration_eta": out_of_calibration}


def scaled(content, target):
    integral = content.sum()
    factor = target / integral if integral > 0 else 0.0
    return content * factor, factor


def root_hist(name, edges, content, title, variance=None):
    hist = scan.ROOT.TH1D(name, f";{title};", len(edges) - 1, array("d", edges.tolist()))
    hist.SetDirectory(0)
    variance = np.asarray(content if variance is None else variance, dtype=float)
    for index, (value, error2) in enumerate(zip(content, variance), 1):
        hist.SetBinContent(index, value); hist.SetBinError(index, math.sqrt(max(error2, 0.0)))
    return hist


def cms_labels(extra="Open Data"):
    cms = scan.ROOT.TLatex(); cms.SetNDC(); cms.SetTextFont(62); cms.SetTextSize(.052)
    cms.DrawLatex(.14, .955, "CMS")
    qualifier = scan.ROOT.TLatex(); qualifier.SetNDC(); qualifier.SetTextFont(52); qualifier.SetTextSize(.039)
    qualifier.DrawLatex(.255, .955, extra)
    energy = scan.ROOT.TLatex(); energy.SetNDC(); energy.SetTextAlign(31); energy.SetTextSize(.034)
    energy.DrawLatex(.95, .955, "2016, #sqrt{s} = 13 TeV")
    return cms, qualifier, energy


def save_svg(canvas, stem):
    """Write the requested vector artifact and remove obsolete raster/PDF copies."""
    canvas.Modified(); canvas.Update()
    svg = stem.with_suffix(".svg")
    for obsolete in (stem.with_suffix(".png"), stem.with_suffix(".pdf")):
        obsolete.unlink(missing_ok=True)
    svg.unlink(missing_ok=True)
    canvas.SaveAs(str(svg))
    canvas.Close()


def draw_distribution(name, samples, output, corrected):
    edges, axis_title, logarithmic, unit = SPECS[name]
    data_raw = samples["data"]["histograms"][name]
    dy_raw = samples["inclusive_dy"]["histograms"][name]
    sd_raw = samples["diffractive"]["histograms"][name]
    target = data_raw.sum()
    dy_content, dy_scale = scaled(dy_raw, target)
    sd_content, sd_scale = scaled(sd_raw, target)
    data = root_hist(f"data_{name}", edges, data_raw, axis_title)
    dy = root_hist(f"dy_{name}", edges, dy_content, axis_title, dy_raw * dy_scale**2)
    sd = root_hist(f"sd_{name}", edges, sd_content, axis_title, sd_raw * sd_scale**2)
    canvas = scan.ROOT.TCanvas(f"canvas_{name}", "", 850, 850)
    canvas.SetFillColor(scan.ROOT.kWhite)
    top = scan.ROOT.TPad(f"top_{name}", "", 0, .30, 1, 1)
    ratio = scan.ROOT.TPad(f"ratio_{name}", "", 0, 0, 1, .30)
    for pad in (top, ratio):
        pad.SetFillColor(scan.ROOT.kWhite); pad.SetLeftMargin(.12); pad.SetRightMargin(.04)
    top.SetBottomMargin(.02); ratio.SetTopMargin(.025); ratio.SetBottomMargin(.34)
    if logarithmic: top.SetLogy()
    top.Draw(); ratio.Draw(); top.cd()
    data.SetMarkerStyle(20); data.SetMarkerSize(.8); data.SetLineColor(scan.ROOT.kBlack)
    dy.SetLineColor(scan.ROOT.kAzure + 2); dy.SetLineWidth(3)
    sd.SetLineColor(scan.ROOT.kRed + 1); sd.SetLineStyle(7); sd.SetLineWidth(3)
    dy.GetXaxis().SetLabelSize(0); dy.GetYaxis().SetTitle(f"{unit} (MC normalized to data)")
    dy.GetYaxis().SetTitleOffset(1.35)
    maximum = max(data.GetMaximum(), dy.GetMaximum(), sd.GetMaximum(), 1.0)
    if logarithmic:
        positive = np.r_[data_raw[data_raw > 0], dy_content[dy_content > 0], sd_content[sd_content > 0]]
        dy.SetMinimum(max(0.5, positive.min() * .4 if len(positive) else .5)); dy.SetMaximum(maximum * 20)
    else:
        dy.SetMinimum(0); dy.SetMaximum(maximum * 1.35)
    dy.Draw("HIST"); sd.Draw("HIST SAME"); data.Draw("E1 SAME")
    legend = scan.ROOT.TLegend(.54, .69, .94, .88); legend.SetBorderSize(0); legend.SetFillStyle(0)
    legend.AddEntry(data, "Data", "ep")
    suffix = " + #mu correction" if corrected else ""
    legend.AddEntry(dy, f"Inclusive DY MC{suffix}", "l")
    legend.AddEntry(sd, f"Single-diffractive MC{suffix} (shape)", "l"); legend.Draw()
    labels = cms_labels()

    ratio.cd()
    ratio_hist = data.Clone(f"ratio_{name}"); ratio_hist.Divide(dy)
    ratio_hist.SetTitle(""); ratio_hist.GetYaxis().SetTitle("Data / DY")
    ratio_hist.GetYaxis().SetRangeUser(0.0, 2.0); ratio_hist.GetYaxis().SetNdivisions(505)
    ratio_hist.GetYaxis().SetTitleSize(.12); ratio_hist.GetYaxis().SetTitleOffset(.45)
    ratio_hist.GetYaxis().SetLabelSize(.095); ratio_hist.GetXaxis().SetTitleSize(.13)
    ratio_hist.GetXaxis().SetLabelSize(.105); ratio_hist.SetMarkerStyle(20); ratio_hist.SetMarkerSize(.7)
    ratio_hist.Draw("E1")
    unity = scan.ROOT.TLine(edges[0], 1, edges[-1], 1); unity.SetLineStyle(2); unity.Draw()
    save_svg(canvas, output / f"{name}_cms_data_mc")
    return {"data_in_range": target, "dy_in_range": float(dy_raw.sum()),
            "signal_in_range": float(sd_raw.sum()), "dy_scale_to_data": dy_scale,
            "signal_scale_to_data": sd_scale}


def roc_curve(signal_scores, background_scores):
    signal_scores = np.asarray(signal_scores, dtype=float)
    background_scores = np.asarray(background_scores, dtype=float)
    scores = np.r_[signal_scores, background_scores]
    labels = np.r_[np.ones(len(signal_scores)), np.zeros(len(background_scores))]
    order = np.argsort(scores, kind="mergesort")[::-1]
    scores, labels = scores[order], labels[order]
    distinct = np.where(np.diff(scores))[0]
    thresholds = np.r_[distinct, len(scores) - 1]
    true_positive = np.cumsum(labels)[thresholds]
    false_positive = 1 + thresholds - true_positive
    tpr = np.r_[0.0, true_positive / max(len(signal_scores), 1)]
    fpr = np.r_[0.0, false_positive / max(len(background_scores), 1)]
    auc = float(np.trapz(tpr, fpr))
    return fpr, tpr, auc


def draw_roc(samples, output, corrected):
    canvas = scan.ROOT.TCanvas("roc_canvas", "", 800, 720)
    canvas.SetFillColor(scan.ROOT.kWhite); canvas.SetLeftMargin(.13); canvas.SetBottomMargin(.12)
    canvas.SetLogy()
    frame = scan.ROOT.TH1D("roc_frame", ";Single-diffractive efficiency;Inclusive DY efficiency", 1, 0, 1)
    frame.SetDirectory(0); frame.SetMinimum(.001); frame.SetMaximum(1.2); frame.Draw("AXIS")
    colors = {"rapidity_gap": scan.ROOT.kGreen + 2, "efg_score": scan.ROOT.kMagenta + 2}
    graphs, metrics = [], {}
    legend = scan.ROOT.TLegend(.43, .68, .93, .88); legend.SetBorderSize(0); legend.SetFillStyle(0)
    for observable in ("rapidity_gap", "efg_score"):
        fpr, tpr, auc = roc_curve(samples["diffractive"]["scores"][observable],
                                  samples["inclusive_dy"]["scores"][observable])
        graph = scan.ROOT.TGraph(len(tpr), tpr, fpr); graph.SetLineColor(colors[observable]); graph.SetLineWidth(3)
        graph.Draw("L SAME"); graphs.append(graph)
        label = "Rapidity gap" if observable == "rapidity_gap" else "Energy Gap Flow"
        legend.AddEntry(graph, f"{label}, AUC = {auc:.3f}", "l")
        metrics[observable] = {"auc": auc,
            "dy_efficiency_at_50pct_signal": float(np.interp(.5, tpr, fpr)),
            "dy_efficiency_at_80pct_signal": float(np.interp(.8, tpr, fpr)),
            "dy_efficiency_at_90pct_signal": float(np.interp(.9, tpr, fpr))}
    diagonal = scan.ROOT.TLine(.001, .001, 1, 1); diagonal.SetLineColor(scan.ROOT.kGray + 1)
    diagonal.SetLineStyle(2); diagonal.Draw(); legend.AddEntry(diagonal, "Random ordering", "l")
    legend.Draw(); labels = cms_labels("Simulation")
    note = scan.ROOT.TLatex(); note.SetNDC(); note.SetTextSize(.026)
    correction_note = "; corrected muons" if corrected else ""
    note.DrawLatex(.15, .17, f"ROC uses MC truth labels{correction_note}; data cannot enter this curve")
    save_svg(canvas, output / "roc_cms_simulation")
    return metrics


def draw_efficiencies(samples, output, corrected):
    """Data-inclusive threshold efficiencies corresponding to both ROC scores."""
    for observable, color in (("rapidity_gap", scan.ROOT.kGreen + 2),
                              ("efg_score", scan.ROOT.kMagenta + 2)):
        canvas = scan.ROOT.TCanvas(f"eff_{observable}", "", 800, 700)
        canvas.SetFillColor(scan.ROOT.kWhite); canvas.SetLeftMargin(.12); canvas.SetBottomMargin(.16)
        graphs = []
        legend = scan.ROOT.TLegend(.53, .69, .93, .88); legend.SetBorderSize(0); legend.SetFillStyle(0)
        xmax = 5.0 if observable == "rapidity_gap" else 1.0
        thresholds = np.linspace(0, xmax, 201)
        suffix = " + #mu correction" if corrected else ""
        for key, label, style, line_color in (
                ("data", "Data", 1, scan.ROOT.kBlack),
                ("inclusive_dy", f"Inclusive DY MC{suffix}", 1, scan.ROOT.kAzure + 2),
                ("diffractive", f"Single-diffractive MC{suffix}", 7, scan.ROOT.kRed + 1)):
            values = np.asarray(samples[key]["scores"][observable])
            efficiency = np.asarray([np.count_nonzero(values >= cut) / max(len(values), 1)
                                     for cut in thresholds])
            graph = scan.ROOT.TGraph(len(thresholds), thresholds, efficiency)
            graph.SetLineColor(line_color); graph.SetLineStyle(style); graph.SetLineWidth(3)
            option = "AL" if not graphs else "L SAME"; graph.Draw(option)
            if not graphs:
                graph.SetTitle(f";minimum {SPECS[observable][1]};fraction above threshold")
                graph.SetMinimum(0); graph.SetMaximum(1.08); graph.GetXaxis().SetLimits(0, xmax)
                graph.GetXaxis().SetTitleSize(.045); graph.GetXaxis().SetTitleOffset(1.25)
            graphs.append(graph); legend.AddEntry(graph, label, "l")
        legend.Draw(); labels = cms_labels()
        save_svg(canvas, output / f"{observable}_threshold_efficiency_data_mc")


def write_root(samples, output):
    target = scan.ROOT.TFile(str(output / "colleague_reproduction_histograms.root"), "RECREATE")
    for sample, result in samples.items():
        for observable, content in result["histograms"].items():
            hist = root_hist(f"{sample}_{observable}", SPECS[observable][0], content,
                             SPECS[observable][1]); hist.Write()
    target.Close()


def main():
    start = time.perf_counter(); args = arguments()
    file_limits = {"data": 2, "inclusive_dy": 3, "diffractive": 3} if args.notebook_file_subset else {
        "data": -1, "inclusive_dy": -1, "diffractive": -1}
    directories = {"data": args.data_dir, "inclusive_dy": args.background_dir,
                   "diffractive": args.signal_dir}
    files = {key: scan.baseline.discover(directory, file_limits[key])
             for key, directory in directories.items()}
    calibration = Calibration.from_summary(args.calibration) if args.apply_muon_corrections else None
    samples = {key: process(files[key], args.max_events, key, calibration,
                            args.correction_seeds) for key in directories}
    output = args.output_dir / args.label; output.mkdir(parents=True, exist_ok=True)
    plot_integrals = {name: draw_distribution(name, samples, output,
                                              args.apply_muon_corrections) for name in SPECS}
    roc = draw_roc(samples, output, args.apply_muon_corrections)
    draw_efficiencies(samples, output, args.apply_muon_corrections); write_root(samples, output)
    summary = {
        "purpose": ("Independent reproduction of Analysis.ipynb with fitted MC muon corrections"
                    if args.apply_muon_corrections else
                    "Independent reproduction of Analysis.ipynb without muon corrections"),
        "configuration": {"selection": {**SELECTION,
                                          "apply_muon_scale_or_smearing": args.apply_muon_corrections},
                          "candidates": CANDIDATES,
                          "notebook_file_subset": args.notebook_file_subset,
                          "muon_correction": {"applied_to": ["inclusive_dy", "diffractive"]
                                              if args.apply_muon_corrections else [],
                                              "calibration_file": str(args.calibration.resolve()),
                                              "parameters": calibration.to_dict() if calibration else None,
                                              "seeds": args.correction_seeds if calibration else []},
                          "normalization": "Each MC shape independently normalized to the data in-range integral"},
        "samples": {key: {"role": {"data": "recorded 2016H data",
                                     "inclusive_dy": "non-diffractive Z/DY background simulation",
                                     "diffractive": "single-diffractive pomflux signal simulation"}[key],
                                  "files": len(files[key]), "entries": value["entries"],
                                  "processed": value["processed"], "cutflow": value["cutflow"],
                                  "cleaned_candidates": value["cleaned_candidates"],
                                  "first_two_not_pt_ordered": value["first_two_not_pt_ordered"],
                                  "muon_replicas_outside_calibration_eta":
                                      value["muon_replicas_outside_calibration_eta"]}
                    for key, value in samples.items()},
        "plot_integrals": plot_integrals, "roc": roc,
        "timing_seconds": time.perf_counter() - start,
    }
    (output / "colleague_reproduction_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"output": str(output), "selected": {k: v["cutflow"]["selected"] for k, v in samples.items()},
                      "roc": roc, "timing_seconds": summary["timing_seconds"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
