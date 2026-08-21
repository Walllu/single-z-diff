#!/usr/bin/env python3
"""Scan muon scale/resolution hypotheses in signed leading-subleading eta regions.

    python scan_zmumu_scale_resolution.py --mode point --regions all --scale 0 --resolution 0 --region-subdirectories --output-dir plots

"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import sys
import time
from array import array
from pathlib import Path

try:
    import ROOT
except ImportError as exc:
    raise SystemExit("PyROOT is required. Run with: conda run -n mg5 python scan_zmumu_scale_resolution.py") from exc

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)

HERE = Path(__file__).resolve().parent
Z_MASS = 91.1876
HISTOGRAM_CONFIGURATION_LABEL = "tiered_v1"


def default_raw_data_directory() -> Path:
    """Locate the workspace-level HackathonDataRaw directory."""
    for parent in (HERE, *HERE.parents):
        candidate = parent / "HackathonDataRaw"
        if candidate.is_dir():
            return candidate
    return HERE.parents[3] / "HackathonDataRaw"

# -----------------------------------------------------------------------------
# User configuration: selection, eta regions, scan grid, histograms, RNG seeds
# -----------------------------------------------------------------------------
SELECTION = {
    "require_trigger": True,
    "trigger_branch": "IsoMu24",
    "muon_min_pt_gev": 15.0,
    "leading_min_pt_gev": 26.0,
    "muon_max_abs_eta": 2.5,
    "require_tight_id": True,
    "require_pf_muon": True,
    "max_abs_dxy_cm": 0.05,
    "max_abs_dz_cm": 0.10,
    "max_delta_beta_relative_isolation": 0.15,
    "fit_mass_range_gev": [70.0, 110.0],
}

# Five signed regions produce 25 ordered (nominal leading, nominal subleading)
# eta-eta categories. Intervals are [low, high), except the final upper edge.
ETA_REGIONS = {
    "neg_endcap": [-2.5, -2.1],
    "neg_transition": [-2.1, -1.4],
    "barrel": [-1.4, 1.4],
    "pos_transition": [1.4, 2.1],
    "pos_endcap": [2.1, 2.5],
}
ETA_PAIR_REGIONS = [f"{lead}__{sublead}" for lead, sublead in itertools.product(ETA_REGIONS, repeat=2)]
ETA_ABBREVIATIONS = {
    "neg_endcap": "NE", "neg_transition": "NT", "barrel": "B",
    "pos_transition": "PT", "pos_endcap": "PE",
}

# Eta-pair categories excluded from normal point/scan running. These are the
# eight ordered positive-side + negative-side combinations without a barrel
# muon. Use --include-excluded-regions to override this list for diagnostics.
EXCLUDED_ETA_PAIR_ABBREVIATIONS = {
    "NE_PT", "NE_PE", "NT_PT", "NT_PE",
    "PT_NE", "PT_NT", "PE_NE", "PE_NT",
}

# Fine common scale shift (-0.5% to +0.5% in 0.05% steps) and wider added
# fractional resolution (0% to 5% in 0.25% steps). Each simulated muon receives
# an independent N(0,1) draw. There are 21 x 21 = 441 hypotheses per region.
SCAN_GRID = {
    "scale_shifts": [-0.0050, -0.0045, -0.0040, -0.0035, -0.0030,
                     -0.0025, -0.0020, -0.0015, -0.0010, -0.0005, 0.0000,
                      0.0005,  0.0010,  0.0015,  0.0020,  0.0025,
                      0.0030,  0.0035,  0.0040,  0.0045,  0.0050],
    "resolution_smearings": [0.0000, 0.0025, 0.0050, 0.0075, 0.0100,
                              0.0125, 0.0150, 0.0175, 0.0200, 0.0225,
                              0.0250, 0.0275, 0.0300, 0.0325, 0.0350,
                              0.0375, 0.0400, 0.0425, 0.0450, 0.0475, 0.0500],
}

DEFAULT_HISTOGRAMS = {
    "mass": [80, 70.0, 110.0, "m_{#mu#mu} [GeV]"],
    "lead_pt": [50, 15.0, 115.0, "leading muon p_{T} [GeV]"],
    "sublead_pt": [40, 15.0, 95.0, "subleading muon p_{T} [GeV]"],
    "muon_eta": [50, -2.5, 2.5, "muon #eta"],
    "z_pt": [50, 0.0, 100.0, "Z candidate p_{T} [GeV]"],
    "nvertices": [50, 0.0, 50.0, "number of vertices"],
}
# Add entries such as "neg_endcap__barrel": {"mass": [80,70,110,"..."]}.
REGION_HISTOGRAM_OVERRIDES = {
    #"barrel__barrel": {"mass": [40, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},
    #"barrel__pos_transition": {"mass": [20, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},
    
    "neg_transition__neg_transition": {"mass": [40, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},  # NT_NT
    "pos_transition__pos_transition": {"mass": [40, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},  # PT_PT
    "barrel__neg_endcap": {"mass": [40, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},              # B_NE
    "neg_endcap__barrel": {"mass": [40, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},              # NE_B
    "barrel__pos_endcap": {"mass": [40, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},              # B_PE
    "pos_endcap__barrel": {"mass": [40, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},              # PE_B

    "pos_endcap__pos_endcap": {"mass": [20, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},      # PE_PE
    "pos_endcap__pos_transition": {"mass": [20, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},  # PE_PT
    "pos_transition__pos_endcap": {"mass": [20, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},  # PT_PE
    "neg_endcap__neg_endcap": {"mass": [20, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},      # NE_NE
    "neg_endcap__neg_transition": {"mass": [20, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},  # NE_NT
    "neg_transition__neg_endcap": {"mass": [20, 70.0, 110.0, "m_{#mu#mu} [GeV]"]},  # NT_NE
}

# One named seed is selected with --seed-name. The table makes alternatives
# explicit and permits repeatable robustness tests without editing algorithms.
RNG_SEEDS = {
    "nominal": 314159,
    "validation_a": 271828,
    "validation_b": 161803,
    "validation_c": 141421,
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    raw = default_raw_data_directory()
    parser.add_argument("--data-dir", type=Path, default=raw / "2016H")
    parser.add_argument("--mc-dir", type=Path,
                        default=raw / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8")
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--mode", choices=("point", "scan"), default="point")
    parser.add_argument("--regions", nargs="+", default=["barrel__barrel"],
                        help="One or more ordered eta-pair names, or 'all'")
    parser.add_argument("--metric", choices=("chi2", "poisson"), default="chi2")
    parser.add_argument("--scale", type=float, default=0.0,
                        help="Fractional common scale shift in point mode")
    parser.add_argument("--resolution", type=float, default=0.0,
                        help="Additional fractional Gaussian width in point mode")
    parser.add_argument("--seed-name", choices=tuple(RNG_SEEDS), default="nominal")
    parser.add_argument("--max-events", type=int, default=-1, help="Per sample; -1 means all")
    parser.add_argument("--max-files", type=int, default=-1, help="Per sample; -1 means all")
    parser.add_argument("--region-subdirectories", action="store_true",
                        help="Nest outputs in a configuration directory and NE_NT-style region directories")
    parser.add_argument("--include-excluded-regions", action="store_true",
                        help="Process requested regions even if listed in EXCLUDED_ETA_PAIR_ABBREVIATIONS")
    return parser.parse_args()


def filesystem_number(value: float) -> str:
    """Format a signed decimal as a compact, shell-friendly directory token."""
    sign = "p" if value >= 0 else "m"
    return sign + f"{abs(value):.3f}".replace(".", "p")


def configuration_directory(args: argparse.Namespace) -> str:
    if args.mode == "point":
        return f"scale_{filesystem_number(args.scale)}_resolution_{filesystem_number(args.resolution)}"
    scales, resolutions = SCAN_GRID["scale_shifts"], SCAN_GRID["resolution_smearings"]
    return (f"scan_{args.metric}_seed_{args.seed_name}_"
            f"s{len(scales)}_{filesystem_number(min(scales))}_{filesystem_number(max(scales))}_"
            f"r{len(resolutions)}_{filesystem_number(min(resolutions))}_{filesystem_number(max(resolutions))}_"
            f"bins_{HISTOGRAM_CONFIGURATION_LABEL}")


def abbreviated_region(region: str) -> str:
    lead, sublead = region.split("__")
    return f"{ETA_ABBREVIATIONS[lead]}_{ETA_ABBREVIATIONS[sublead]}"


def discover(directory: Path, maximum: int) -> list[Path]:
    files = sorted(directory.glob("*.root"))
    if maximum >= 0:
        files = files[:maximum]
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
        raise RuntimeError(f"No usable pfExtractor/pfTree inputs in {directory}")
    return usable


def make_chain(files: list[Path]):
    tree = ROOT.TChain("pfExtractor/pfTree")
    for filename in files:
        tree.Add(str(filename))
    tree.SetBranchStatus("*", 0)
    for branch in (
        "run", "lumi", "event", "nMuons", "nVertices", "IsoMu24", "muonPt",
        "muonEta", "muonPhi", "muonMass", "muonCharge", "muonDxy", "muonDz",
        "muonIsoCharged", "muonIsoNeutral", "muonIsoPhoton", "muonIsoPU",
        "muonIsTight", "muonIsPF"
    ):
        tree.SetBranchStatus(branch, 1)
    return tree


def relative_isolation(tree, i: int) -> float:
    neutral = tree.muonIsoNeutral[i] + tree.muonIsoPhoton[i]
    numerator = tree.muonIsoCharged[i] + max(0.0, neutral - 0.5 * tree.muonIsoPU[i])
    return numerator / tree.muonPt[i] if tree.muonPt[i] > 0 else math.inf


def selected_muons(tree) -> list[int]:
    result = []
    for i in range(tree.nMuons):
        if tree.muonPt[i] <= SELECTION["muon_min_pt_gev"]:
            continue
        if abs(tree.muonEta[i]) >= SELECTION["muon_max_abs_eta"]:
            continue
        if SELECTION["require_tight_id"] and not bool(tree.muonIsTight[i]):
            continue
        if SELECTION["require_pf_muon"] and not bool(tree.muonIsPF[i]):
            continue
        if abs(tree.muonDxy[i]) >= SELECTION["max_abs_dxy_cm"]:
            continue
        if abs(tree.muonDz[i]) >= SELECTION["max_abs_dz_cm"]:
            continue
        if relative_isolation(tree, i) >= SELECTION["max_delta_beta_relative_isolation"]:
            continue
        result.append(i)
    return result


def p4(pt: float, eta: float, phi: float, mass: float):
    return ROOT.Math.PtEtaPhiMVector(pt, eta, phi, mass)


def dimuon_mass(lead: tuple, sublead: tuple, lead_factor: float, sublead_factor: float) -> float:
    """Exact scalar equivalent of adding two PtEtaPhiM four-vectors."""
    lead_pt, lead_eta, lead_phi, lead_mass = lead
    sub_pt, sub_eta, sub_phi, sub_mass = sublead
    lead_pt *= lead_factor; sub_pt *= sublead_factor
    lead_px, lead_py = lead_pt * math.cos(lead_phi), lead_pt * math.sin(lead_phi)
    sub_px, sub_py = sub_pt * math.cos(sub_phi), sub_pt * math.sin(sub_phi)
    lead_pz, sub_pz = lead_pt * math.sinh(lead_eta), sub_pt * math.sinh(sub_eta)
    lead_energy = math.sqrt((lead_pt * math.cosh(lead_eta)) ** 2 + lead_mass ** 2)
    sub_energy = math.sqrt((sub_pt * math.cosh(sub_eta)) ** 2 + sub_mass ** 2)
    mass2 = ((lead_energy + sub_energy) ** 2 - (lead_px + sub_px) ** 2
             - (lead_py + sub_py) ** 2 - (lead_pz + sub_pz) ** 2)
    return math.sqrt(max(0.0, mass2))


def best_pair(tree, muons: list[int]):
    pairs = []
    for pos, first in enumerate(muons):
        for second in muons[pos + 1:]:
            if tree.muonCharge[first] * tree.muonCharge[second] >= 0:
                continue
            z = p4(tree.muonPt[first], tree.muonEta[first], tree.muonPhi[first], tree.muonMass[first])
            z += p4(tree.muonPt[second], tree.muonEta[second], tree.muonPhi[second], tree.muonMass[second])
            pairs.append((abs(z.M() - Z_MASS), first, second, z.M()))
    return min(pairs, key=lambda item: item[0]) if pairs else None


def eta_region(eta: float) -> str | None:
    names = list(ETA_REGIONS)
    for position, name in enumerate(names):
        low, high = ETA_REGIONS[name]
        if low <= eta < high or (position == len(names) - 1 and eta == high):
            return name
    return None


def read_events(files: list[Path], max_events: int, is_mc: bool, seed: int):
    tree = make_chain(files)
    total = tree.GetEntries()
    limit = total if max_events < 0 else min(total, max_events)
    rng = random.Random(seed)
    events = {region: [] for region in ETA_PAIR_REGIONS}
    flow = {key: 0 for key in ("all", "trigger", "two_muons", "opposite_sign",
                                "leading_pt", "fit_mass", "eta_region")}
    for entry in range(limit):
        tree.GetEntry(entry); flow["all"] += 1
        if SELECTION["require_trigger"] and not bool(getattr(tree, SELECTION["trigger_branch"])):
            continue
        flow["trigger"] += 1
        muons = selected_muons(tree)
        if len(muons) < 2:
            continue
        flow["two_muons"] += 1
        pair = best_pair(tree, muons)
        if pair is None:
            continue
        flow["opposite_sign"] += 1
        _, lead, sublead, nominal_mass = pair
        if tree.muonPt[lead] < tree.muonPt[sublead]:
            lead, sublead = sublead, lead
        if tree.muonPt[lead] <= SELECTION["leading_min_pt_gev"]:
            continue
        flow["leading_pt"] += 1
        low_mass, high_mass = SELECTION["fit_mass_range_gev"]
        if not low_mass < nominal_mass < high_mass:
            continue
        flow["fit_mass"] += 1
        region = f"{eta_region(tree.muonEta[lead])}__{eta_region(tree.muonEta[sublead])}"
        if "None" in region:
            continue
        flow["eta_region"] += 1
        values = {
            "lead": (float(tree.muonPt[lead]), float(tree.muonEta[lead]),
                     float(tree.muonPhi[lead]), float(tree.muonMass[lead])),
            "sublead": (float(tree.muonPt[sublead]), float(tree.muonEta[sublead]),
                        float(tree.muonPhi[sublead]), float(tree.muonMass[sublead])),
            "nvertices": int(tree.nVertices),
            "event_key": f"{int(tree.run)}:{int(tree.lumi)}:{int(tree.event)}",
            "gaussians": (rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)) if is_mc else (0.0, 0.0),
        }
        events[region].append(values)
    return events, flow, total, limit


def histogram_config(region: str) -> dict:
    config = {key: list(value) for key, value in DEFAULT_HISTOGRAMS.items()}
    config.update({key: list(value) for key, value in REGION_HISTOGRAM_OVERRIDES.get(region, {}).items()})
    return config


def make_histograms(prefix: str, region: str):
    result = {name: ROOT.TH1D(f"{prefix}_{region}_{name}", f";{axis};Events", bins, low, high)
              for name, (bins, low, high, axis) in histogram_config(region).items()}
    for histogram in result.values():
        histogram.SetDirectory(0)
    return result


def fill(events: list[dict], prefix: str, region: str, scale: float, resolution: float, is_mc: bool):
    hist = make_histograms(prefix, region)
    accepted = 0
    for event in events:
        lead, sublead = event["lead"], event["sublead"]
        g_lead, g_sublead = event["gaussians"]
        lead_factor = 1.0 + scale + resolution * g_lead if is_mc else 1.0
        sublead_factor = 1.0 + scale + resolution * g_sublead if is_mc else 1.0
        lead4 = p4(lead[0] * lead_factor, lead[1], lead[2], lead[3])
        sublead4 = p4(sublead[0] * sublead_factor, sublead[1], sublead[2], sublead[3])
        z = lead4 + sublead4
        low_mass, high_mass = SELECTION["fit_mass_range_gev"]
        if not low_mass < z.M() < high_mass:
            continue
        accepted += 1
        hist["mass"].Fill(z.M())
        hist["lead_pt"].Fill(lead4.Pt())
        hist["sublead_pt"].Fill(sublead4.Pt())
        hist["muon_eta"].Fill(lead[1]); hist["muon_eta"].Fill(sublead[1])
        hist["z_pt"].Fill(z.Pt()); hist["nvertices"].Fill(event["nvertices"])
    return hist, accepted


def fill_mass(events: list[dict], prefix: str, region: str, scale: float, resolution: float, is_mc: bool):
    """Fill only the fit observable, avoiding five unnecessary fills per grid point."""
    bins, low, high, axis = histogram_config(region)["mass"]
    hist = ROOT.TH1D(f"{prefix}_{region}_mass", f";{axis};Events", bins, low, high)
    hist.SetDirectory(0)
    accepted = 0
    for event in events:
        lead, sublead = event["lead"], event["sublead"]
        g_lead, g_sublead = event["gaussians"]
        lead_factor = 1.0 + scale + resolution * g_lead if is_mc else 1.0
        sublead_factor = 1.0 + scale + resolution * g_sublead if is_mc else 1.0
        mass = dimuon_mass(lead, sublead, lead_factor, sublead_factor)
        low_mass, high_mass = SELECTION["fit_mass_range_gev"]
        if low_mass < mass < high_mass:
            hist.Fill(mass); accepted += 1
    return hist, accepted


def normalized_mc(data, mc):
    result = mc.Clone(mc.GetName() + "_normalized")
    alpha = data.Integral() / mc.Integral() if mc.Integral() > 0 else 0.0
    result.Scale(alpha)
    return result, alpha


def objective(data, mc, metric: str):
    normalized, alpha = normalized_mc(data, mc)
    value, bins_used = 0.0, 0
    for b in range(1, data.GetNbinsX() + 1):
        observed, expected = data.GetBinContent(b), normalized.GetBinContent(b)
        if metric == "chi2":
            variance = data.GetBinError(b) ** 2 + (alpha * mc.GetBinError(b)) ** 2
            if variance <= 0:
                continue
            value += (observed - expected) ** 2 / variance
        else:
            expected = max(expected, 1e-12)
            value += 2.0 * (expected - observed + (observed * math.log(observed / expected) if observed > 0 else 0.0))
        bins_used += 1
    # One degree of freedom is consumed by shape normalization.
    return value, max(0, bins_used - 1), alpha


def draw_comparison(data, mc, name: str, output: Path, label: str):
    canvas = ROOT.TCanvas("canvas_" + name, "", 800, 800)
    top = ROOT.TPad("top_" + name, "", 0, 0.29, 1, 1)
    bottom = ROOT.TPad("bottom_" + name, "", 0, 0, 1, 0.29)
    top.SetBottomMargin(0.02); bottom.SetTopMargin(0.03); bottom.SetBottomMargin(0.34)
    top.Draw(); bottom.Draw(); top.cd()
    d = data.Clone("draw_data_" + name); m, _ = normalized_mc(d, mc)
    d.SetMarkerStyle(20); d.SetMarkerSize(0.75); m.SetLineColor(ROOT.kAzure + 2)
    m.SetLineWidth(2); m.SetFillColorAlpha(ROOT.kAzure - 9, 0.45)
    m.SetMaximum(1.35 * max(d.GetMaximum(), m.GetMaximum(), 1.0)); m.SetMinimum(0)
    m.GetYaxis().SetTitle("Events (MC normalized to data)"); m.GetXaxis().SetLabelSize(0)
    m.Draw("HIST"); d.Draw("E1 SAME")
    legend = ROOT.TLegend(0.61, 0.72, 0.88, 0.87); legend.SetBorderSize(0)
    legend.AddEntry(d, "Data", "lep"); legend.AddEntry(m, "smeared simulation", "lf"); legend.Draw()
    text = ROOT.TLatex(); text.SetNDC(); text.SetTextSize(0.033); text.DrawLatex(0.13, 0.92, label)
    bottom.cd(); ratio = d.Clone("ratio_" + name); ratio.Divide(m)
    ratio.SetTitle(""); ratio.GetYaxis().SetTitle("Data / MC"); ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio.GetYaxis().SetNdivisions(505); ratio.GetYaxis().SetTitleSize(0.105)
    ratio.GetYaxis().SetTitleOffset(0.45); ratio.GetYaxis().SetLabelSize(0.09)
    ratio.GetXaxis().SetTitleSize(0.12); ratio.GetXaxis().SetLabelSize(0.10); ratio.Draw("E1")
    line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2); line.Draw()
    canvas.SaveAs(str(output / f"{name}.png")); canvas.SaveAs(str(output / f"{name}.pdf"))


def draw_scan_surface(region: str, results: list[dict], best: dict, output: Path) -> None:
    """Draw a readable delta-objective landscape from stored grid results."""
    scales, resolutions = SCAN_GRID["scale_shifts"], SCAN_GRID["resolution_smearings"]
    scale_half_step = 0.5 * (scales[1] - scales[0])
    resolution_half_step = 0.5 * (resolutions[1] - resolutions[0])
    surface = ROOT.TH2D("surface_" + region, ";scale shift;added resolution", len(scales),
                        min(scales) - scale_half_step, max(scales) + scale_half_step,
                        len(resolutions), min(resolutions) - resolution_half_step,
                        max(resolutions) + resolution_half_step)
    for item in results:
        surface.Fill(item["scale"], item["resolution"], item["objective"] - best["objective"])
    surface.SetTitle(f"{abbreviated_region(region)};scale shift;added resolution;#Delta objective")
    surface.SetMinimum(0.0); surface.SetMaximum(30.0)
    canvas = ROOT.TCanvas("surface_canvas_" + region, "", 850, 700); canvas.SetFillColor(ROOT.kWhite)
    ROOT.gPad.SetFillColor(ROOT.kWhite)
    canvas.SetRightMargin(0.18); surface.GetZaxis().SetTitleOffset(1.25); surface.Draw("COLZ")
    contours = surface.Clone("contours_" + region)
    contours.SetContour(3, array("d", [2.30, 6.18, 11.83]))
    contours.SetLineColor(ROOT.kBlack); contours.SetLineWidth(2); contours.Draw("CONT3 SAME")
    marker = ROOT.TGraph(1, array("d", [best["scale"]]), array("d", [best["resolution"]]))
    marker.SetMarkerStyle(20); marker.SetMarkerSize(1.5); marker.SetMarkerColor(ROOT.kRed + 1); marker.Draw("P SAME")
    legend = ROOT.TLegend(0.13, 0.77, 0.42, 0.89); legend.SetBorderSize(0); legend.SetFillStyle(0)
    legend.AddEntry(marker, "best grid point", "p"); legend.AddEntry(contours, "#Delta = 2.30, 6.18, 11.83", "l"); legend.Draw()
    canvas.SaveAs(str(output / f"{region}_scan.png")); canvas.SaveAs(str(output / f"{region}_scan.pdf"))


def scan_region(region: str, data_events, mc_events, metric: str, output: Path):
    data_hist, _ = fill(data_events, "data", region, 0.0, 0.0, False)
    results, best = [], None
    for scale in SCAN_GRID["scale_shifts"]:
        for resolution in SCAN_GRID["resolution_smearings"]:
            mc_mass, accepted = fill_mass(mc_events, "mc", region, scale, resolution, True)
            value, ndof, alpha = objective(data_hist["mass"], mc_mass, metric)
            item = {"scale": scale, "resolution": resolution, "objective": value,
                    "ndof": ndof, "normalization": alpha, "mc_accepted": accepted}
            results.append(item)
            if best is None or value < best["objective"]:
                best = item
    best_hist, _ = fill(mc_events, "best_mc", region, best["scale"], best["resolution"], True)
    for name in data_hist:
        draw_comparison(data_hist[name], best_hist[name], f"{region}_{name}", output,
                        f"{region}: scale={best['scale']:+.3f}, resolution={best['resolution']:.3f}")
    draw_scan_surface(region, results, best, output)
    return results, best, data_hist, best_hist


def draw_parameter_summary(summary: dict, output: Path) -> None:
    """Draw 5x5 ordered eta-pair maps, masking excluded or unavailable cells."""
    names = list(ETA_REGIONS); labels = [ETA_ABBREVIATIONS[name] for name in names]
    canvas = ROOT.TCanvas("parameter_summary", "", 1300, 620); canvas.SetFillColor(ROOT.kWhite); canvas.Divide(2, 1)
    keepalive = []
    for pad, key, title in ((1, "scale", "Best scale shift [%]"),
                            (2, "resolution", "Best added resolution [%]")):
        canvas.cd(pad); ROOT.gPad.SetFillColor(ROOT.kWhite)
        ROOT.gPad.SetRightMargin(0.16); ROOT.gPad.SetLeftMargin(0.13)
        hist = ROOT.TH2D(f"summary_{key}", f"{title};subleading-muon region;leading-muon region",
                         5, 0, 5, 5, 0, 5); hist.SetDirectory(0); hist.SetMarkerSize(1.35)
        for i, label in enumerate(labels, 1):
            hist.GetXaxis().SetBinLabel(i, label); hist.GetYaxis().SetBinLabel(i, label)
        values = []
        for iy, lead in enumerate(names, 1):
            for ix, sublead in enumerate(names, 1):
                region = f"{lead}__{sublead}"
                if region in summary["regions"]:
                    value = 100.0 * summary["regions"][region]["best"][key]
                    # ROOT leaves exactly-zero TH2 bins unpainted. Store a
                    # negligible epsilon for display and annotate the true value.
                    hist.SetBinContent(ix, iy, value if value != 0 else 1e-12); values.append(value)
        if key == "scale" and values:
            extent = max(abs(min(values)), abs(max(values)), 0.01); hist.SetMinimum(-extent); hist.SetMaximum(extent)
        elif values:
            hist.SetMinimum(0.0); hist.SetMaximum(max(values) * 1.15)
        hist.Draw("COLZ")
        boxes, texts = [], []
        for iy, lead in enumerate(names, 1):
            for ix, sublead in enumerate(names, 1):
                region = f"{lead}__{sublead}"
                if region in summary["regions"]:
                    value = 100.0 * summary["regions"][region]["best"][key]
                    text = ROOT.TLatex(ix - 0.5, iy - 0.5, f"{value:g}")
                    text.SetTextAlign(22); text.SetTextSize(0.032); text.SetTextFont(62)
                    text.SetTextColor(ROOT.kBlack); text.Draw("SAME")
                    texts.append(text)
                    continue
                box = ROOT.TBox(ix - 1, iy - 1, ix, iy); box.SetFillColor(ROOT.kGray + 1); box.Draw("SAME")
                text = ROOT.TLatex(ix - 0.5, iy - 0.5, "X"); text.SetTextAlign(22); text.SetTextSize(0.06)
                text.SetTextColor(ROOT.kBlack); text.Draw("SAME")
                boxes.append(box); texts.append(text)
        keepalive.extend([hist, boxes, texts])
    canvas.Update(); canvas.SaveAs(str(output / "best_fit_parameter_grid.png"))
    canvas.SaveAs(str(output / "best_fit_parameter_grid.pdf"))


def main() -> int:
    start_time = time.perf_counter()
    args = arguments()
    requested_regions = ETA_PAIR_REGIONS if args.regions == ["all"] else args.regions
    unknown = sorted(set(requested_regions) - set(ETA_PAIR_REGIONS))
    if unknown:
        raise SystemExit(f"Unknown regions: {unknown}. Valid choices: {ETA_PAIR_REGIONS}")
    excluded_regions = [region for region in requested_regions
                        if abbreviated_region(region) in EXCLUDED_ETA_PAIR_ABBREVIATIONS]
    regions = requested_regions if args.include_excluded_regions else [
        region for region in requested_regions if region not in excluded_regions]
    if args.resolution < 0:
        raise SystemExit("--resolution must be non-negative")
    output_root = args.output_dir / configuration_directory(args) if args.region_subdirectories else args.output_dir
    output_root.mkdir(parents=True, exist_ok=True)
    data_files, mc_files = discover(args.data_dir, args.max_files), discover(args.mc_dir, args.max_files)
    seed = RNG_SEEDS[args.seed_name]
    extraction_start = time.perf_counter()
    data_events, data_flow, data_total, data_read = read_events(data_files, args.max_events, False, seed)
    mc_events, mc_flow, mc_total, mc_read = read_events(mc_files, args.max_events, True, seed)
    extraction_seconds = time.perf_counter() - extraction_start
    summary = {"configuration": {"mode": args.mode, "metric": args.metric, "regions": regions,
                "excluded_region_abbreviations": sorted(EXCLUDED_ETA_PAIR_ABBREVIATIONS),
                "excluded_from_this_run": [abbreviated_region(region) for region in excluded_regions]
                    if not args.include_excluded_regions else [],
                "seed_name": args.seed_name, "seed": seed, "selection": SELECTION,
                "scan_grid": SCAN_GRID, "histogram_configuration_label": HISTOGRAM_CONFIGURATION_LABEL,
                "scale": args.scale, "resolution": args.resolution},
               "inputs": {"data_files": len(data_files), "mc_files": len(mc_files),
                "data_entries": data_total, "mc_entries": mc_total,
                "data_processed": data_read, "mc_processed": mc_read},
               "cutflows": {"data": data_flow, "mc": mc_flow}, "regions": {}}
    root_file = ROOT.TFile(str(output_root / "histograms.root"), "RECREATE")
    for region in regions:
        region_output = output_root / abbreviated_region(region) if args.region_subdirectories else output_root
        region_output.mkdir(parents=True, exist_ok=True)
        if args.mode == "scan":
            results, best, data_hist, mc_hist = scan_region(region, data_events[region], mc_events[region], args.metric, region_output)
            summary["regions"][region] = {"data_selected": len(data_events[region]),
                "mc_selected": len(mc_events[region]), "best": best, "grid_results": results}
        else:
            data_hist, data_count = fill(data_events[region], "data", region, 0.0, 0.0, False)
            mc_hist, mc_count = fill(mc_events[region], "mc", region, args.scale, args.resolution, True)
            value, ndof, alpha = objective(data_hist["mass"], mc_hist["mass"], args.metric)
            summary["regions"][region] = {"data_selected": data_count, "mc_selected": mc_count,
                "objective": value, "ndof": ndof, "normalization": alpha}
            label = (f"{region}: scale={args.scale:+.3f}, resolution={args.resolution:.3f}, "
                     f"N_{{data}}={data_count}, N_{{MC}}={mc_count}")
            for name in data_hist:
                draw_comparison(data_hist[name], mc_hist[name], f"{region}_{name}", region_output, label)
        root_file.cd()
        for collection in (data_hist, mc_hist):
            for histogram in collection.values():
                histogram.Write()
    root_file.Close()
    if args.mode == "scan":
        draw_parameter_summary(summary, output_root)
    summary["configuration"]["output_directory"] = str(output_root)
    summary["configuration"]["region_subdirectories"] = args.region_subdirectories
    summary["timing_seconds"] = {"event_extraction": extraction_seconds,
                                 "total": time.perf_counter() - start_time}
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    brief = {"configuration": summary["configuration"], "inputs": summary["inputs"],
             "cutflows": summary["cutflows"],
             "regions": {name: {key: value for key, value in result.items() if key != "grid_results"}
                         for name, result in summary["regions"].items()}}
    print(json.dumps(brief, indent=2)); print(f"Wrote full results to {output_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
