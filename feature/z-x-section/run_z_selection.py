#!/usr/bin/env python3
"""Inclusive Z->mumu event selection and charge/isolation ABCD counting.

This is deliberately a reconstruction-level first pass.  Acceptance and
truth-matching are kept for a later step so the data selection and the ABCD
control regions can first be frozen and validated independently.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

try:
    import ROOT
except ImportError as exc:
    raise SystemExit(
        "PyROOT is required. Create the environment documented in environment.yml."
    ) from exc


HERE = Path(__file__).resolve().parent
TREE_NAME = "pfExtractor/pfTree"

# -----------------------------------------------------------------------------
# Analysis configuration. Command-line options can override the numerical
# thresholds without changing the checked-in defaults.
# -----------------------------------------------------------------------------
SELECTION = {
    "require_trigger": True,
    "trigger_branch": "IsoMu24",
    "require_trigger_match": True,
    "muon_min_pt_gev": 25.0,
    "muon_max_abs_eta": 2.4,
    "require_tight_id": True,
    "require_pf_muon": True,
    "max_abs_dxy_cm": 0.05,
    "max_abs_dz_cm": 0.10,
    "mass_window_gev": [91.1876 - 15.0, 91.1876 + 15.0],
    "additional_lepton_veto": {
        "enabled": True,
        "minimum_pt_gev": 10.0,
        "electron_definition": (
            "PV-associated PF candidate with candIsElectron != 0 "
            "and candFromPV >= 2 (loose-electron proxy)"
        ),
    },
}

ISOLATION = {
    "definition": "delta-beta PF relative isolation",
    "isolated_max": 0.15,
    "anti_isolated_min": 0.25,
    # "both" follows the proposed two-isolated versus two-non-isolated ABCD.
    # "at_least_one" is available if the D region is too sparse.
    "anti_isolation_mode": "both",
}

# Region names avoid reusing A, which already denotes detector acceptance.
ABCD_REGIONS = {
    "SR": {"code": 1, "charge": "opposite-sign", "isolation": "both isolated"},
    "B": {"code": 2, "charge": "same-sign", "isolation": "both isolated"},
    "C": {"code": 3, "charge": "opposite-sign", "isolation": "anti-isolated"},
    "D": {"code": 4, "charge": "same-sign", "isolation": "anti-isolated"},
    "unassigned": {"code": 0, "charge": "either", "isolation": "mixed/transition"},
}


CPP_HELPERS = r"""
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace zxs {
using ROOT::VecOps::RVec;

RVec<int> select_kinematic(const RVec<float>& pt, const RVec<float>& eta,
                           double min_pt, double max_abs_eta) {
    RVec<int> result;
    for (std::size_t i = 0; i < pt.size(); ++i) {
        if (pt[i] > min_pt && std::abs(eta[i]) < max_abs_eta) result.push_back(i);
    }
    std::sort(result.begin(), result.end(),
              [&](int a, int b) { return pt[a] > pt[b]; });
    return result;
}

RVec<int> select_quality(const RVec<float>& pt, const RVec<float>& eta,
                         const RVec<int>& tight, const RVec<int>& pf,
                         const RVec<float>& dxy, const RVec<float>& dz,
                         double min_pt, double max_abs_eta,
                         bool require_tight, bool require_pf,
                         double max_abs_dxy, double max_abs_dz) {
    RVec<int> result;
    for (std::size_t i = 0; i < pt.size(); ++i) {
        if (!(pt[i] > min_pt && std::abs(eta[i]) < max_abs_eta)) continue;
        if (require_tight && !tight[i]) continue;
        if (require_pf && !pf[i]) continue;
        if (!(std::abs(dxy[i]) < max_abs_dxy && std::abs(dz[i]) < max_abs_dz)) continue;
        result.push_back(i);
    }
    std::sort(result.begin(), result.end(),
              [&](int a, int b) { return pt[a] > pt[b]; });
    return result;
}

double relative_isolation(const RVec<float>& pt,
                          const RVec<float>& charged,
                          const RVec<float>& neutral,
                          const RVec<float>& photon,
                          const RVec<float>& pileup,
                          int index) {
    if (index < 0 || static_cast<std::size_t>(index) >= pt.size() || pt[index] <= 0.)
        return 1.e9;
    const double corrected_neutral = std::max(
        0.0, static_cast<double>(neutral[index] + photon[index] - 0.5f * pileup[index]));
    return (charged[index] + corrected_neutral) / pt[index];
}

double dimuon_mass(const RVec<float>& pt, const RVec<float>& eta,
                   const RVec<float>& phi, const RVec<float>& mass,
                   const RVec<int>& indices) {
    using P4 = ROOT::Math::PtEtaPhiMVector;
    const int first = indices[0], second = indices[1];
    const P4 one(pt[first], eta[first], phi[first], mass[first]);
    const P4 two(pt[second], eta[second], phi[second], mass[second]);
    return (one + two).M();
}

int isolation_state(double first, double second,
                    double isolated_max, double anti_isolated_min,
                    bool require_both_anti_isolated) {
    if (first < isolated_max && second < isolated_max) return 1;
    const bool anti = require_both_anti_isolated
        ? (first > anti_isolated_min && second > anti_isolated_min)
        : (first > anti_isolated_min || second > anti_isolated_min);
    return anti ? -1 : 0;
}

int abcd_region(int charge_product, int iso_state) {
    if (charge_product < 0 && iso_state == 1) return 1;  // signal region
    if (charge_product > 0 && iso_state == 1) return 2;  // B
    if (charge_product < 0 && iso_state == -1) return 3; // C
    if (charge_product > 0 && iso_state == -1) return 4; // D
    return 0;
}

bool passes_additional_lepton_veto(const RVec<float>& muon_pt,
                                   const RVec<int>& muon_loose,
                                   const RVec<int>& selected_indices,
                                   const RVec<float>& candidate_pt,
                                   const RVec<int>& candidate_is_electron,
                                   const RVec<int>& candidate_from_pv,
                                   double minimum_pt) {
    const int first = selected_indices[0], second = selected_indices[1];
    for (std::size_t i = 0; i < muon_pt.size(); ++i) {
        if (static_cast<int>(i) == first || static_cast<int>(i) == second) continue;
        if (muon_pt[i] > minimum_pt && muon_loose[i] != 0) return false;
    }
    for (std::size_t i = 0; i < candidate_pt.size(); ++i)
        if (candidate_pt[i] > minimum_pt && candidate_is_electron[i] != 0
            && candidate_from_pv[i] >= 2) return false;
    return true;
}
} // namespace zxs
"""


def default_raw_data_directory() -> Path:
    for parent in (HERE, *HERE.parents):
        candidate = parent / "HackathonDataRaw"
        if candidate.is_dir():
            return candidate
    return HERE.parents[2] / "HackathonDataRaw"


def arguments() -> argparse.Namespace:
    raw = default_raw_data_directory()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=raw / "2016H")
    parser.add_argument(
        "--data-files-from", type=Path,
        help="Read data ROOT paths from a manifest (one path per line) instead of --data-dir",
    )
    parser.add_argument(
        "--mc-dir", type=Path,
        default=raw / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8",
    )
    parser.add_argument(
        "--mc-files-from", type=Path,
        help="Read MC ROOT paths from a manifest (one path per line) instead of --mc-dir",
    )
    parser.add_argument("--sample", choices=("data", "mc", "both"), default="both")
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    parser.add_argument("--label", default="local_selection")
    parser.add_argument("--max-files", type=int, default=-1, help="Per sample")
    parser.add_argument("--max-events", type=int, default=-1, help="Per sample")
    parser.add_argument("--threads", type=int, default=1,
                        help="RDataFrame threads; use 1 with --max-events")
    parser.add_argument("--muon-min-pt", type=float, default=SELECTION["muon_min_pt_gev"])
    parser.add_argument("--muon-max-abs-eta", type=float,
                        default=SELECTION["muon_max_abs_eta"])
    parser.add_argument("--mass-min", type=float, default=SELECTION["mass_window_gev"][0])
    parser.add_argument("--mass-max", type=float, default=SELECTION["mass_window_gev"][1])
    parser.add_argument("--isolated-max", type=float, default=ISOLATION["isolated_max"])
    parser.add_argument("--anti-isolated-min", type=float,
                        default=ISOLATION["anti_isolated_min"])
    parser.add_argument("--anti-isolation-mode", choices=("both", "at_least_one"),
                        default=ISOLATION["anti_isolation_mode"])
    parser.add_argument("--no-trigger", action="store_true",
                        help="Disable the IsoMu24 requirement for diagnostics")
    parser.add_argument(
        "--no-trigger-match", action="store_true",
        help="Disable the requirement that at least one selected muon matches the trigger",
    )
    parser.add_argument(
        "--no-additional-lepton-veto", action="store_true",
        help="Disable the extra loose-muon/PF-electron-candidate veto",
    )
    parser.add_argument("--write-skim", action="store_true",
                        help="Write a compact ROOT tree for requested ABCD regions")
    parser.add_argument(
        "--skim-regions", nargs="+", choices=tuple(ABCD_REGIONS),
        default=list(ABCD_REGIONS),
        help="ABCD assignments retained in the skim; default also keeps mixed/transition events",
    )
    args = parser.parse_args()
    if not 0.0 < args.isolated_max < args.anti_isolated_min:
        parser.error("Require 0 < --isolated-max < --anti-isolated-min")
    if not 0.0 < args.mass_min < args.mass_max:
        parser.error("Require 0 < --mass-min < --mass-max")
    if args.max_events >= 0 and args.threads != 1:
        parser.error("ROOT RDataFrame Range is not supported with implicit multithreading; use --threads 1")
    return args


def required_branches(is_mc: bool) -> set[str]:
    required = {
        "run", "lumi", "event", "nMuons", "IsoMu24", "muonPt", "muonEta",
        "muonPhi", "muonMass", "muonCharge", "muonDxy", "muonDz",
        "muonIsoCharged", "muonIsoNeutral", "muonIsoPhoton", "muonIsoPU",
        "muonIsLoose", "muonIsTight", "muonIsPF", "muonIsTrigMatched",
        "candPt", "candIsElectron", "candFromPV",
    }
    if is_mc:
        required.add("genWeight")
    return required


def manifest_files(filename: Path) -> list[Path]:
    """Load deterministic local/CephFS paths from a newline-delimited manifest."""
    if not filename.is_file():
        raise FileNotFoundError(f"Input manifest does not exist: {filename}")
    files: list[Path] = []
    for line_number, raw_line in enumerate(filename.read_text().splitlines(), start=1):
        value = raw_line.strip()
        if not value or value.startswith("#"):
            continue
        expanded = Path(os.path.expandvars(value)).expanduser()
        if not expanded.is_absolute():
            expanded = filename.resolve().parent / expanded
        if expanded.suffix.lower() != ".root":
            raise ValueError(
                f"Manifest {filename}, line {line_number}: expected a .root path, got {value!r}"
            )
        files.append(expanded)
    if not files:
        raise RuntimeError(f"Input manifest contains no ROOT files: {filename}")
    return files


def validate_candidates(candidates: list[Path], maximum: int,
                        is_mc: bool) -> tuple[list[Path], list[dict[str, str]]]:
    if maximum >= 0:
        candidates = candidates[:maximum]
    usable: list[Path] = []
    rejected: list[dict[str, str]] = []
    for filename in candidates:
        source = ROOT.TFile.Open(str(filename), "READ")
        tree = source.Get(TREE_NAME) if source and not source.IsZombie() else None
        if tree and tree.InheritsFrom("TTree"):
            available = {branch.GetName() for branch in tree.GetListOfBranches()}
            missing = sorted(required_branches(is_mc) - available)
            if missing:
                reason = f"missing branches: {', '.join(missing)}"
                print(
                    f"WARNING: skipping schema-incompatible file {filename}; "
                    f"missing {', '.join(missing)}",
                    file=sys.stderr,
                )
                rejected.append({"file": str(filename.resolve()), "reason": reason})
            else:
                usable.append(filename)
        else:
            print(f"WARNING: skipping file without {TREE_NAME}: {filename}", file=sys.stderr)
            rejected.append({
                "file": str(filename.resolve()),
                "reason": f"missing or invalid {TREE_NAME}",
            })
        if source:
            source.Close()
    if not usable:
        raise RuntimeError("No usable ROOT files found in the requested input")
    return usable, rejected


def discover(directory: Path, maximum: int, is_mc: bool) -> tuple[list[Path], list[dict[str, str]]]:
    if not directory.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {directory}")
    return validate_candidates(sorted(directory.glob("*.root")), maximum, is_mc)


def discover_input(directory: Path, manifest: Path | None, maximum: int,
                   is_mc: bool) -> tuple[list[Path], list[dict[str, str]]]:
    """Discover inputs from a directory or an explicit HTCondor-safe manifest."""
    if manifest is None:
        return discover(directory, maximum, is_mc)
    return validate_candidates(manifest_files(manifest), maximum, is_mc)


def make_chain(files: list[Path]) -> Any:
    chain = ROOT.TChain(TREE_NAME)
    for filename in files:
        chain.Add(str(filename))
    return chain


def validate_branches(chain: Any, is_mc: bool) -> None:
    available = {branch.GetName() for branch in chain.GetListOfBranches()}
    missing = sorted(required_branches(is_mc) - available)
    if missing:
        raise RuntimeError(f"Missing required branches: {', '.join(missing)}")


def weighted_actions(node: Any) -> tuple[Any, Any, Any]:
    return node.Count(), node.Sum("event_weight"), node.Sum("event_weight_sq")


def as_number(result: Any, integer: bool = False) -> int | float:
    value = result.GetValue()
    return int(value) if integer else float(value)


def abcd_estimate(regions: dict[str, dict[str, float]]) -> dict[str, float | None]:
    b, c, d = (regions[name]["weighted_yield"] for name in ("B", "C", "D"))
    vb, vc, vd = (regions[name]["sum_weights_squared"] for name in ("B", "C", "D"))
    if b <= 0.0 or c <= 0.0 or d <= 0.0:
        return {"estimate": None, "stat_uncertainty": None, "formula": "B*C/D"}
    estimate = b * c / d
    variance = 0.0
    if b != 0.0:
        variance += vb / (b * b)
    if c != 0.0:
        variance += vc / (c * c)
    variance += vd / (d * d)
    uncertainty = abs(estimate) * math.sqrt(max(0.0, variance))
    return {
        "estimate": estimate,
        "stat_uncertainty": uncertainty,
        "formula": "B*C/D",
    }


def write_cutflow_csv(filename: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "sample", "step", "events", "weighted_yield", "sum_weights_squared",
        "step_efficiency", "cumulative_efficiency", "weighted_step_efficiency",
        "weighted_cumulative_efficiency",
    ]
    with filename.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def process_sample(name: str, files: list[Path], rejected_files: list[dict[str, str]],
                   args: argparse.Namespace, output: Path) -> dict[str, Any]:
    is_mc = name == "mc"
    chain = make_chain(files)
    validate_branches(chain, is_mc)
    entries = int(chain.GetEntries())
    dataframe = ROOT.RDataFrame(chain)
    if args.max_events >= 0:
        dataframe = dataframe.Range(min(args.max_events, entries))

    dataframe = dataframe.Define(
        "event_weight", "static_cast<double>(genWeight)" if is_mc else "1.0"
    ).Define("event_weight_sq", "event_weight * event_weight")

    require_trigger = SELECTION["require_trigger"] and not args.no_trigger
    require_trigger_match = (
        SELECTION["require_trigger_match"]
        and not args.no_trigger
        and not args.no_trigger_match
    )
    trigger_expression = f'{SELECTION["trigger_branch"]} != 0' if require_trigger else "true"
    tight = str(SELECTION["require_tight_id"]).lower()
    pf = str(SELECTION["require_pf_muon"]).lower()

    nodes: list[tuple[str, Any]] = [("all", dataframe)]
    node = dataframe.Filter(trigger_expression, "trigger")
    nodes.append(("trigger", node))
    node = node.Filter("nMuons >= 2", "at least two reconstructed muons")
    nodes.append(("two_reconstructed_muons", node))
    node = node.Define(
        "kinematic_muon_indices",
        f"zxs::select_kinematic(muonPt, muonEta, {args.muon_min_pt}, {args.muon_max_abs_eta})",
    ).Filter("kinematic_muon_indices.size() >= 2", "at least two kinematic muons")
    nodes.append(("two_kinematic_muons", node))
    node = node.Define(
        "selected_muon_indices",
        "zxs::select_quality(muonPt, muonEta, muonIsTight, muonIsPF, muonDxy, muonDz, "
        f"{args.muon_min_pt}, {args.muon_max_abs_eta}, {tight}, {pf}, "
        f'{SELECTION["max_abs_dxy_cm"]}, {SELECTION["max_abs_dz_cm"]})',
    ).Filter("selected_muon_indices.size() >= 2", "at least two quality muons")
    nodes.append(("two_quality_muons", node))

    require_additional_veto = (
        SELECTION["additional_lepton_veto"]["enabled"]
        and not args.no_additional_lepton_veto
    )
    veto_expression = (
        "zxs::passes_additional_lepton_veto(muonPt, muonIsLoose, "
        "selected_muon_indices, candPt, candIsElectron, candFromPV, "
        f'{SELECTION["additional_lepton_veto"]["minimum_pt_gev"]})'
        if require_additional_veto else "true"
    )
    node = node.Filter(veto_expression, "additional loose-lepton veto")
    nodes.append(("additional_lepton_veto", node))

    node = (
        node.Define("lead_index", "selected_muon_indices[0]")
        .Define("sublead_index", "selected_muon_indices[1]")
        .Define("lead_pt", "static_cast<double>(muonPt[lead_index])")
        .Define("sublead_pt", "static_cast<double>(muonPt[sublead_index])")
        .Define("lead_eta", "static_cast<double>(muonEta[lead_index])")
        .Define("sublead_eta", "static_cast<double>(muonEta[sublead_index])")
        .Define("lead_charge", "muonCharge[lead_index]")
        .Define("sublead_charge", "muonCharge[sublead_index]")
        .Define("lead_trigger_matched", "muonIsTrigMatched[lead_index] != 0")
        .Define("sublead_trigger_matched", "muonIsTrigMatched[sublead_index] != 0")
        .Define("charge_product", "lead_charge * sublead_charge")
    )
    trigger_match_expression = (
        "lead_trigger_matched || sublead_trigger_matched"
        if require_trigger_match else "true"
    )
    node = node.Filter(trigger_match_expression, "selected pair trigger match")
    nodes.append(("selected_pair_trigger_match", node))

    node = (
        node.Define(
            "dimuon_mass",
            "zxs::dimuon_mass(muonPt, muonEta, muonPhi, muonMass, selected_muon_indices)",
        )
        .Filter(f"dimuon_mass > {args.mass_min} && dimuon_mass < {args.mass_max}",
                "dimuon mass window")
    )
    nodes.append(("mass_window", node))

    anti_both = str(args.anti_isolation_mode == "both").lower()
    classified = (
        node.Define(
            "lead_rel_iso",
            "zxs::relative_isolation(muonPt, muonIsoCharged, muonIsoNeutral, "
            "muonIsoPhoton, muonIsoPU, lead_index)",
        )
        .Define(
            "sublead_rel_iso",
            "zxs::relative_isolation(muonPt, muonIsoCharged, muonIsoNeutral, "
            "muonIsoPhoton, muonIsoPU, sublead_index)",
        )
        .Define(
            "isolation_state",
            f"zxs::isolation_state(lead_rel_iso, sublead_rel_iso, {args.isolated_max}, "
            f"{args.anti_isolated_min}, {anti_both})",
        )
        .Define("abcd_region_code", "zxs::abcd_region(charge_product, isolation_state)")
    )

    region_nodes = {
        region: classified.Filter(f"abcd_region_code == {config['code']}", region)
        for region, config in ABCD_REGIONS.items()
    }
    nodes.append(("signal_region", region_nodes["SR"]))

    action_sets = {step: weighted_actions(current) for step, current in nodes}
    region_actions = {region: weighted_actions(current) for region, current in region_nodes.items()}
    all_actions = [action for actions in action_sets.values() for action in actions]
    all_actions.extend(action for actions in region_actions.values() for action in actions)

    skim_action = None
    skim_file = None
    if args.write_skim:
        requested_codes = [ABCD_REGIONS[region]["code"] for region in args.skim_regions]
        expression = " || ".join(f"abcd_region_code == {code}" for code in requested_codes)
        skim_node = classified.Filter(expression, "requested skim regions")
        skim_file = output / f"{name}_abcd_skim.root"
        columns = [
            "run", "lumi", "event", "event_weight", "lead_index", "sublead_index",
            "lead_pt", "sublead_pt", "lead_eta", "sublead_eta", "lead_charge",
            "sublead_charge", "lead_trigger_matched", "sublead_trigger_matched",
            "lead_rel_iso", "sublead_rel_iso", "dimuon_mass", "abcd_region_code",
        ]
        options = ROOT.RDF.RSnapshotOptions()
        options.fMode = "RECREATE"
        options.fLazy = True
        skim_action = skim_node.Snapshot("Events", str(skim_file), columns, options)
        all_actions.append(skim_action)

    ROOT.RDF.RunGraphs(all_actions)

    cutflow: list[dict[str, Any]] = []
    previous = None
    previous_weighted = None
    initial = None
    initial_weighted = None
    for step, actions in action_sets.items():
        count = as_number(actions[0], integer=True)
        weighted = as_number(actions[1])
        sumw2 = as_number(actions[2])
        if initial is None:
            initial = count
            initial_weighted = weighted
        cutflow.append({
            "sample": name,
            "step": step,
            "events": count,
            "weighted_yield": weighted,
            "sum_weights_squared": sumw2,
            "step_efficiency": count / previous if previous else (1.0 if count else 0.0),
            "cumulative_efficiency": count / initial if initial else 0.0,
            "weighted_step_efficiency": (
                weighted / previous_weighted if previous_weighted else (1.0 if weighted else 0.0)
            ),
            "weighted_cumulative_efficiency": (
                weighted / initial_weighted if initial_weighted else 0.0
            ),
        })
        previous = count
        previous_weighted = weighted

    regions: dict[str, dict[str, Any]] = {}
    for region, actions in region_actions.items():
        regions[region] = {
            **ABCD_REGIONS[region],
            "events": as_number(actions[0], integer=True),
            "weighted_yield": as_number(actions[1]),
            "sum_weights_squared": as_number(actions[2]),
        }

    return {
        "sample": name,
        "is_mc": is_mc,
        "files": [str(path.resolve()) for path in files],
        "files_used": len(files),
        "rejected_files": rejected_files,
        "tree_entries": entries,
        "entries_requested": args.max_events,
        "cutflow": cutflow,
        "regions": regions,
        "raw_abcd_nonprompt_estimate_in_sr": abcd_estimate(regions),
        "skim": str(skim_file.resolve()) if skim_file else None,
    }


def main() -> int:
    args = arguments()
    ROOT.gROOT.SetBatch(True)
    ROOT.gInterpreter.Declare(CPP_HELPERS)
    if args.threads > 1:
        ROOT.EnableImplicitMT(args.threads)

    output = args.output_dir / args.label
    output.mkdir(parents=True, exist_ok=True)
    configuration = {
        "selection": {
            **SELECTION,
            "require_trigger": SELECTION["require_trigger"] and not args.no_trigger,
            "require_trigger_match": (
                SELECTION["require_trigger_match"]
                and not args.no_trigger
                and not args.no_trigger_match
            ),
            "muon_min_pt_gev": args.muon_min_pt,
            "muon_max_abs_eta": args.muon_max_abs_eta,
            "mass_window_gev": [args.mass_min, args.mass_max],
            "pair_choice": "two highest-pT quality muons, independent of charge and isolation",
            "additional_lepton_veto": {
                **SELECTION["additional_lepton_veto"],
                "enabled": (
                    SELECTION["additional_lepton_veto"]["enabled"]
                    and not args.no_additional_lepton_veto
                ),
            },
        },
        "isolation": {
            **ISOLATION,
            "isolated_max": args.isolated_max,
            "anti_isolated_min": args.anti_isolated_min,
            "anti_isolation_mode": args.anti_isolation_mode,
        },
        "abcd_regions": ABCD_REGIONS,
    }

    sample_inputs: list[tuple[str, Path, Path | None]] = []
    if args.sample in ("data", "both"):
        sample_inputs.append(("data", args.data_dir, args.data_files_from))
    if args.sample in ("mc", "both"):
        sample_inputs.append(("mc", args.mc_dir, args.mc_files_from))

    started = time.monotonic()
    results: dict[str, Any] = {}
    all_cutflow: list[dict[str, Any]] = []
    for name, directory, manifest in sample_inputs:
        files, rejected = discover_input(directory, manifest, args.max_files, name == "mc")
        print(f"Processing {name}: {len(files)} files", flush=True)
        result = process_sample(name, files, rejected, args, output)
        results[name] = result
        all_cutflow.extend(result["cutflow"])
        estimate = result["raw_abcd_nonprompt_estimate_in_sr"]
        print(
            f"  SR={result['regions']['SR']['events']:,}, "
            f"B={result['regions']['B']['events']:,}, "
            f"C={result['regions']['C']['events']:,}, "
            f"D={result['regions']['D']['events']:,}, "
            f"raw ABCD={estimate['estimate']}",
            flush=True,
        )

    summary = {
        "title": "Inclusive Z to dimuon reconstruction and ABCD selection",
        "configuration": configuration,
        "execution": {
            "sample": args.sample,
            "max_files_per_sample": args.max_files,
            "max_events_per_sample": args.max_events,
            "threads": args.threads,
            "wall_time_seconds": time.monotonic() - started,
        },
        "samples": results,
        "caveats": [
            "The ABCD estimate is raw: prompt contamination has not been subtracted.",
            "Charge/isolation factorization and sideband closure have not yet been tested.",
            "No luminosity normalization, truth acceptance, or reconstruction efficiency is applied.",
            "An arbitrary data-file subset cannot be paired with the full Run2016H luminosity.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_cutflow_csv(output / "cutflow.csv", all_cutflow)
    print(f"Wrote {output / 'summary.json'} and {output / 'cutflow.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
