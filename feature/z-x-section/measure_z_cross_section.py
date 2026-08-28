#!/usr/bin/env python3
"""Corrected reconstruction, dressed-Z truth acceptance, and efficiency study."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

import run_z_selection as core

ROOT = core.ROOT
HERE = Path(__file__).resolve().parent


def arguments() -> argparse.Namespace:
    raw = core.default_raw_data_directory()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "measurement_config.json")
    parser.add_argument("--data-dir", type=Path, default=raw / "2016H")
    parser.add_argument("--data-files-from", type=Path,
                        help="Read data ROOT paths from a manifest instead of --data-dir")
    parser.add_argument(
        "--mc-dir", type=Path,
        default=raw / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8",
    )
    parser.add_argument("--mc-files-from", type=Path,
                        help="Read MC ROOT paths from a manifest instead of --mc-dir")
    parser.add_argument("--sample", choices=("data", "mc", "both"), default="both")
    parser.add_argument("--mc-role", choices=("signal", "background"), default="signal",
                        help="Interpret the MC as the Z signal or as a prompt background")
    parser.add_argument("--process-name", default="z_to_mumu",
                        help="Stable process label recorded in the output")
    parser.add_argument("--defer-normalization", action="store_true",
                        help="Keep MC yields in raw genWeight units for distributed merging")
    parser.add_argument("--max-files", type=int, default=-1, help="Per sample")
    parser.add_argument("--max-events", type=int, default=-1, help="Per sample")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    parser.add_argument("--label", default="local_corrected_measurement")
    parser.add_argument("--no-muon-momentum", action="store_true")
    parser.add_argument("--no-muon-efficiency", action="store_true")
    parser.add_argument("--enable-experimental-pileup", action="store_true")
    parser.add_argument("--truth-mass-window", nargs=2, type=float, metavar=("MIN", "MAX"))
    parser.add_argument("--reco-mass-window", nargs=2, type=float, metavar=("MIN", "MAX"))
    parser.add_argument("--write-selected-skim", action="store_true")
    args = parser.parse_args()
    if args.max_events >= 0 and args.threads != 1:
        parser.error("Use --threads 1 with --max-events")
    return args


def resolve(config_file: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_file.resolve().parent / path).resolve()


def load_sf_module(directory: Path) -> Any:
    filename = directory / "MuonPerformance_EfficiencyCorrections.py"
    spec = importlib.util.spec_from_file_location("zxs_muon_efficiency", filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import scale-factor helper {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module._load_sf_maps()
    return module


def cpp_array(values: np.ndarray) -> str:
    rows = []
    for row in np.asarray(values, dtype=float):
        rows.append("{" + ",".join(f"{item:.17g}" for item in row) + "}")
    return "{" + ",".join(rows) + "}"


def pileup_table(filename: Path) -> np.ndarray:
    payload = json.loads(filename.read_text())
    maximum = int(payload["configuration"]["overflow_vertex"])
    result = np.ones((3, maximum + 1), dtype=float)
    for ir, region in enumerate(("BB", "BE", "EE")):
        covered = np.zeros(maximum + 1, dtype=bool)
        for group in payload["regions"][region]:
            lo = max(0, int(group["minimum"]))
            hi = min(maximum, int(group["maximum"]))
            result[ir, lo:hi + 1] = float(group["weight"])
            covered[lo:hi + 1] = True
        if not covered.all():
            raise ValueError(f"Pileup payload does not cover every nVertices value in {region}")
    return result


def build_cpp(config: dict[str, Any], config_file: Path) -> str:
    corrections = config["corrections"]
    momentum = json.loads(resolve(config_file, corrections["muon_momentum"]["payload"]).read_text())
    sf_dir = resolve(config_file, corrections["muon_efficiency"]["payload_directory"])
    sf = load_sf_module(sf_dir)
    pu_file = resolve(config_file, corrections["pileup"]["payload"])
    pu = pileup_table(pu_file)
    golden_json = config["normalization"].get("golden_json")
    if golden_json is None:
        certified_lumi = "bool certified_lumi(unsigned int, unsigned int) { return true; }"
    else:
        lumis = json.loads(resolve(config_file, golden_json).read_text())
        cases = []
        for run, ranges in sorted(lumis.items(), key=lambda item: int(item[0])):
            condition = " || ".join(
                f"(lumi >= {int(bounds[0])}U && lumi <= {int(bounds[1])}U)"
                for bounds in ranges
            ) or "false"
            cases.append(f"case {int(run)}U: return {condition};")
        certified_lumi = (
            "bool certified_lumi(unsigned int run, unsigned int lumi) { "
            "switch (run) {" + "".join(cases) + "default: return false;} }"
        )
    maps = {
        "reco_nom": sf._MuonReco_nominal_sf,
        "reco_sys": sf._MuonReco_systematic_sf,
        "reco_stat": sf._MuonReco_statistic_sf,
        "trig_nom": sf._MuonTrigger_nominal_sf,
        "trig_sys": sf._MuonTrigger_systematic_sf,
        "trig_stat": sf._MuonTrigger_statistic_sf,
        "iso_nom": sf._MuonIso_nominal_sf,
        "iso_sys": sf._MuonIso_systematic_sf,
        "iso_stat": sf._MuonIso_statistic_sf,
    }
    declarations = "\n".join(
        f"static const double {name}[10][15] = {cpp_array(values)};"
        for name, values in maps.items()
    )
    pu_rows = "{" + ",".join(
        "{" + ",".join(f"{x:.17g}" for x in row) + "}" for row in pu
    ) + "}"
    truth = config["physics_definition"]
    fid = truth["truth_fiducial"]
    tmass = truth["truth_mass_window_gev"]
    return f"""
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace zxm {{
using ROOT::VecOps::RVec;
using P4 = ROOT::Math::PxPyPzEVector;
{declarations}
static const double pu_weights[3][{pu.shape[1]}] = {pu_rows};
{certified_lumi}

struct TruthZ {{
  bool exists = false;
  bool fiducial = false;
  double mass = -1.;
  double lead_pt = -1., sublead_pt = -1.;
  double lead_eta = 0., sublead_eta = 0.;
  double lead_phi = 0., sublead_phi = 0.;
  int lead_pdgid = 0, sublead_pdgid = 0;
}};

uint64_t mix(uint64_t x) {{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}}

double gaussian(uint64_t run, uint64_t lumi, uint64_t event, uint64_t index, uint64_t seed) {{
  uint64_t key = mix(seed ^ mix(run) ^ mix(lumi << 16) ^ mix(event) ^ mix(index + 1));
  uint64_t key2 = mix(key);
  const double u1 = (static_cast<double>((key >> 11) + 1)) / 9007199254740993.0;
  const double u2 = (static_cast<double>((key2 >> 11) + 1)) / 9007199254740993.0;
  return std::sqrt(-2. * std::log(u1)) * std::cos(2. * M_PI * u2);
}}

RVec<float> corrected_pts(const RVec<float>& pt, const RVec<float>& eta,
                          unsigned int run, unsigned int lumi, unsigned long long event,
                          bool apply, unsigned long long seed) {{
  RVec<float> result(pt);
  if (!apply) return result;
  for (std::size_t i = 0; i < result.size(); ++i) {{
    const bool barrel = std::abs(eta[i]) < 1.4;
    const double scale = barrel ? {momentum['barrel_scale']:.17g} : {momentum['endcap_scale']:.17g};
    const double resolution = barrel ? {momentum['barrel_resolution']:.17g} : {momentum['endcap_resolution']:.17g};
    result[i] *= 1. + scale + resolution * gaussian(run, lumi, event, i, seed);
  }}
  return result;
}}

int eta_bin(double eta) {{ return std::clamp(static_cast<int>(std::floor((eta + 2.5) / 0.5)), 0, 9); }}
int pt_bin(double pt) {{ return std::clamp(static_cast<int>(std::floor((pt - 25.) / 5.)), 0, 14); }}
double lookup(const double table[10][15], double pt, double eta) {{ return table[eta_bin(eta)][pt_bin(pt)]; }}
double varied(const double nom[10][15], const double unc[10][15], double pt, double eta, int direction) {{
  return lookup(nom, pt, eta) + direction * lookup(unc, pt, eta);
}}
int detector_region(double eta1, double eta2) {{
  const bool b1 = std::abs(eta1) < 1.4, b2 = std::abs(eta2) < 1.4;
  return b1 && b2 ? 0 : ((!b1 && !b2) ? 2 : 1);
}}
double pileup(int region, int nvertices, bool apply) {{
  if (!apply) return 1.;
  return pu_weights[std::clamp(region, 0, 2)][std::clamp(nvertices, 0, {pu.shape[1]-1})];
}}

bool has_ancestor(int index, int target, const RVec<int>& pdgid, const RVec<int>& mother) {{
  int current = index;
  for (int depth = 0; depth < 100; ++depth) {{
    if (current < 0 || current >= static_cast<int>(mother.size())) return false;
    current = mother[current];
    if (current < 0 || current >= static_cast<int>(pdgid.size())) return false;
    if (std::abs(pdgid[current]) == std::abs(target)) return true;
  }}
  return false;
}}
double delta_phi(double a, double b) {{ return std::remainder(a - b, 2. * M_PI); }}
double delta_r(double e1, double p1, double e2, double p2) {{
  return std::hypot(e1 - e2, delta_phi(p1, p2));
}}
P4 p4(double pt, double eta, double phi, double mass) {{
  const double px = pt * std::cos(phi), py = pt * std::sin(phi), pz = pt * std::sinh(eta);
  return P4(px, py, pz, std::sqrt(px*px + py*py + pz*pz + mass*mass));
}}

TruthZ dressed_truth_z(const RVec<float>& pt, const RVec<float>& eta,
                       const RVec<float>& phi, const RVec<float>& mass,
                       const RVec<int>& pdgid, const RVec<int>& status,
                       const RVec<int>& mother) {{
  std::vector<int> muons;
  for (std::size_t i = 0; i < pdgid.size(); ++i)
    if (std::abs(pdgid[i]) == 13 && status[i] == 1 && has_ancestor(i, 23, pdgid, mother)) muons.push_back(i);
  std::vector<P4> dressed;
  for (int index : muons) dressed.push_back(p4(pt[index], eta[index], phi[index], mass[index]));
  for (std::size_t ip = 0; ip < pdgid.size(); ++ip) {{
    if (pdgid[ip] != 22 || status[ip] != 1 || !has_ancestor(ip, 23, pdgid, mother)) continue;
    int best = -1; double best_dr = {truth['dressing_cone_dr']:.17g};
    for (std::size_t im = 0; im < muons.size(); ++im) {{
      const int index = muons[im];
      const double dr = delta_r(eta[ip], phi[ip], eta[index], phi[index]);
      if (dr < best_dr) {{ best = static_cast<int>(im); best_dr = dr; }}
    }}
    if (best >= 0) dressed[best] += p4(pt[ip], eta[ip], phi[ip], 0.);
  }}
  TruthZ result; double closest = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < muons.size(); ++i) for (std::size_t j = i + 1; j < muons.size(); ++j) {{
    if (pdgid[muons[i]] * pdgid[muons[j]] >= 0) continue;
    const double candidate_mass = (dressed[i] + dressed[j]).M();
    const double distance = std::abs(candidate_mass - 91.1876);
    if (distance >= closest) continue;
    closest = distance; result.exists = candidate_mass >= {tmass[0]:.17g} && candidate_mass < {tmass[1]:.17g};
    result.mass = candidate_mass;
    const bool first_leads = dressed[i].Pt() >= dressed[j].Pt();
    const std::size_t lead = first_leads ? i : j, sub = first_leads ? j : i;
    result.lead_pt = dressed[lead].Pt(); result.sublead_pt = dressed[sub].Pt();
    result.lead_eta = dressed[lead].Eta(); result.sublead_eta = dressed[sub].Eta();
    result.lead_phi = dressed[lead].Phi(); result.sublead_phi = dressed[sub].Phi();
    result.lead_pdgid = pdgid[muons[lead]]; result.sublead_pdgid = pdgid[muons[sub]];
    result.fiducial = result.exists && result.lead_pt > {fid['muon_min_pt_gev']:.17g}
      && result.sublead_pt > {fid['muon_min_pt_gev']:.17g}
      && std::abs(result.lead_eta) < {fid['muon_max_abs_eta']:.17g}
      && std::abs(result.sublead_eta) < {fid['muon_max_abs_eta']:.17g};
  }}
  return result;
}}

bool truth_matched(const TruthZ& truth, double eta1, double phi1, int charge1,
                   double eta2, double phi2, int charge2, double maximum_dr) {{
  if (!truth.fiducial) return false;
  const int truth_charge_lead = truth.lead_pdgid > 0 ? -1 : 1;
  const int truth_charge_sub = truth.sublead_pdgid > 0 ? -1 : 1;
  const bool direct = charge1 == truth_charge_lead && charge2 == truth_charge_sub
    && delta_r(eta1, phi1, truth.lead_eta, truth.lead_phi) < maximum_dr
    && delta_r(eta2, phi2, truth.sublead_eta, truth.sublead_phi) < maximum_dr;
  const bool swapped = charge1 == truth_charge_sub && charge2 == truth_charge_lead
    && delta_r(eta1, phi1, truth.sublead_eta, truth.sublead_phi) < maximum_dr
    && delta_r(eta2, phi2, truth.lead_eta, truth.lead_phi) < maximum_dr;
  return direct || swapped;
}}
}}
"""


def action(node: Any, weight: str) -> tuple[Any, Any, Any]:
    return node.Count(), node.Sum(weight), node.Define(f"{weight}_sq_tmp", f"{weight}*{weight}").Sum(f"{weight}_sq_tmp")


def values(actions: tuple[Any, Any, Any]) -> dict[str, float | int]:
    return {
        "events": int(actions[0].GetValue()),
        "sum_weights": float(actions[1].GetValue()),
        "sum_weights_squared": float(actions[2].GetValue()),
    }


def safe_ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator != 0.0 else None


def report_rows(report: Any) -> list[dict[str, float | int | str]]:
    result = []
    for cut in report.GetValue():
        result.append({
            "step": str(cut.GetName()),
            "events": int(cut.GetPass()),
            "events_before": int(cut.GetAll()),
            "step_efficiency_percent": float(cut.GetEff()),
        })
    return result


def normalization(config: dict[str, Any], processed_sumw: float,
                  defer: bool = False) -> tuple[float, dict[str, Any]]:
    norm = config["normalization"]
    denominator = norm["sum_generator_weights"]
    denominator = processed_sumw if denominator is None else float(denominator)
    required = (norm["luminosity_pb_inverse"], norm["dy_cross_section_pb"])
    ready = all(value is not None for value in required) and denominator != 0.0 and not defer
    factor = 1.0
    if ready:
        factor = (float(norm["luminosity_pb_inverse"]) * float(norm["dy_cross_section_pb"])
                  * float(norm["filter_efficiency"]) * float(norm["k_factor"]) / denominator)
    return factor, {
        "absolute_ready": ready,
        "deferred_to_final_stage": defer,
        "factor": factor,
        "processed_sum_generator_weights": processed_sumw,
        "normalization_sum_generator_weights": denominator,
        **norm,
    }


def process_mc(files: list[Path], rejected: list[dict[str, str]], args: argparse.Namespace,
               config: dict[str, Any], output: Path) -> dict[str, Any]:
    chain = core.make_chain(files)
    required_truth = {"genPt", "genEta", "genPhi", "genMass", "genPdgId", "genStatus", "genMotherIndex", "nVertices"}
    available = {branch.GetName() for branch in chain.GetListOfBranches()}
    missing = sorted(required_truth - available)
    if missing:
        raise RuntimeError(f"MC truth/weight branches missing: {', '.join(missing)}")
    rdf = ROOT.RDataFrame(chain)
    if args.max_events >= 0:
        rdf = rdf.Range(min(args.max_events, int(chain.GetEntries())))
    processed_sumw = float(rdf.Sum("genWeight").GetValue())
    norm_factor, norm_info = normalization(config, processed_sumw, args.defer_normalization)
    corr = config["corrections"]
    apply_momentum = corr["muon_momentum"]["enabled"] and not args.no_muon_momentum
    apply_sf = corr["muon_efficiency"]["enabled"] and not args.no_muon_efficiency
    apply_pu = corr["pileup"]["enabled"] or args.enable_experimental_pileup
    seed = int(corr["muon_momentum"]["seed"])
    reco = config["physics_definition"]["reconstruction"]
    mass = reco["mass_window_gev"]

    base = (rdf.Define("gen_base_weight", f"static_cast<double>(genWeight)*{norm_factor:.17g}")
            .Define("truth_z", "zxm::dressed_truth_z(genPt,genEta,genPhi,genMass,genPdgId,genStatus,genMotherIndex)")
            .Define("truth_region", "zxm::detector_region(truth_z.lead_eta,truth_z.sublead_eta)")
            .Define("truth_pu_weight", f"zxm::pileup(truth_region,nVertices,{str(apply_pu).lower()})")
            .Define("truth_eff_den_weight", "gen_base_weight*truth_pu_weight"))
    truth_total = base.Filter("truth_z.exists", "truth Z mass definition")
    truth_fid = truth_total.Filter("truth_z.fiducial", "dressed truth fiducial")

    node = base.Define("analysisMuonPt", f"zxm::corrected_pts(muonPt,muonEta,run,lumi,event,{str(apply_momentum).lower()},{seed}ULL)")
    reconstruction_nodes = {"all": node}
    node = node.Filter("IsoMu24 != 0", "IsoMu24")
    reconstruction_nodes["trigger"] = node
    node = node.Filter("nMuons >= 2", "two reconstructed muons")
    reconstruction_nodes["two_reconstructed_muons"] = node
    node = (node.Define("selected_muon_indices",
                        f"zxs::select_quality(analysisMuonPt,muonEta,muonIsTight,muonIsPF,muonDxy,muonDz,{reco['muon_min_pt_gev']},{reco['muon_max_abs_eta']},true,true,0.05,0.10)")
            .Filter("selected_muon_indices.size() >= 2", "two quality muons"))
    reconstruction_nodes["two_quality_muons"] = node
    node = (node.Define("lead_index", "selected_muon_indices[0]")
            .Define("sublead_index", "selected_muon_indices[1]")
            .Define("lead_pt", "static_cast<double>(analysisMuonPt[lead_index])")
            .Define("sublead_pt", "static_cast<double>(analysisMuonPt[sublead_index])")
            .Define("lead_eta", "static_cast<double>(muonEta[lead_index])")
            .Define("sublead_eta", "static_cast<double>(muonEta[sublead_index])")
            .Define("lead_phi", "static_cast<double>(muonPhi[lead_index])")
            .Define("sublead_phi", "static_cast<double>(muonPhi[sublead_index])")
            .Define("lead_charge", "muonCharge[lead_index]")
            .Define("sublead_charge", "muonCharge[sublead_index]")
            .Define("lead_trigger_matched", "muonIsTrigMatched[lead_index] != 0")
            .Define("sublead_trigger_matched", "muonIsTrigMatched[sublead_index] != 0")
            .Filter("lead_trigger_matched || sublead_trigger_matched", "selected pair trigger match"))
    reconstruction_nodes["selected_pair_trigger_match"] = node
    node = (node.Define("dimuon_mass", "zxs::dimuon_mass(analysisMuonPt,muonEta,muonPhi,muonMass,selected_muon_indices)")
            .Filter(f"dimuon_mass>{mass[0]} && dimuon_mass<{mass[1]}", "reconstructed mass window"))
    reconstruction_nodes["mass_window"] = node
    node = (node
            .Define("lead_rel_iso", "zxs::relative_isolation(analysisMuonPt,muonIsoCharged,muonIsoNeutral,muonIsoPhoton,muonIsoPU,lead_index)")
            .Define("sublead_rel_iso", "zxs::relative_isolation(analysisMuonPt,muonIsoCharged,muonIsoNeutral,muonIsoPhoton,muonIsoPU,sublead_index)")
            .Define("charge_product", "lead_charge*sublead_charge")
            .Define("iso_state_both", "zxs::isolation_state(lead_rel_iso,sublead_rel_iso,0.15,0.25,true)")
            .Define("iso_state_at_least_one", "zxs::isolation_state(lead_rel_iso,sublead_rel_iso,0.15,0.25,false)")
            .Define("selected_signal", "lead_charge*sublead_charge<0 && lead_rel_iso<0.15 && sublead_rel_iso<0.15")
            .Define("truth_reco_matched", f"zxm::truth_matched(truth_z,lead_eta,lead_phi,lead_charge,sublead_eta,sublead_phi,sublead_charge,{reco['truth_match_dr']})"))

    enabled = str(apply_sf).lower()
    node = (node.Define("reco_sf_nom", f"{enabled} ? zxm::lookup(zxm::reco_nom,lead_pt,lead_eta)*zxm::lookup(zxm::reco_nom,sublead_pt,sublead_eta) : 1.")
            .Define("reco_sf_sys_up", f"{enabled} ? zxm::varied(zxm::reco_nom,zxm::reco_sys,lead_pt,lead_eta,1)*zxm::varied(zxm::reco_nom,zxm::reco_sys,sublead_pt,sublead_eta,1) : 1.")
            .Define("reco_sf_sys_down", f"{enabled} ? zxm::varied(zxm::reco_nom,zxm::reco_sys,lead_pt,lead_eta,-1)*zxm::varied(zxm::reco_nom,zxm::reco_sys,sublead_pt,sublead_eta,-1) : 1.")
            .Define("reco_sf_stat_up", f"{enabled} ? zxm::varied(zxm::reco_nom,zxm::reco_stat,lead_pt,lead_eta,1)*zxm::varied(zxm::reco_nom,zxm::reco_stat,sublead_pt,sublead_eta,1) : 1.")
            .Define("reco_sf_stat_down", f"{enabled} ? zxm::varied(zxm::reco_nom,zxm::reco_stat,lead_pt,lead_eta,-1)*zxm::varied(zxm::reco_nom,zxm::reco_stat,sublead_pt,sublead_eta,-1) : 1.")
            .Define("iso_sf_nom", f"{enabled} ? zxm::lookup(zxm::iso_nom,lead_pt,lead_eta)*zxm::lookup(zxm::iso_nom,sublead_pt,sublead_eta) : 1.")
            .Define("iso_sf_sys_up", f"{enabled} ? zxm::varied(zxm::iso_nom,zxm::iso_sys,lead_pt,lead_eta,1)*zxm::varied(zxm::iso_nom,zxm::iso_sys,sublead_pt,sublead_eta,1) : 1.")
            .Define("iso_sf_sys_down", f"{enabled} ? zxm::varied(zxm::iso_nom,zxm::iso_sys,lead_pt,lead_eta,-1)*zxm::varied(zxm::iso_nom,zxm::iso_sys,sublead_pt,sublead_eta,-1) : 1.")
            .Define("iso_sf_stat_up", f"{enabled} ? zxm::varied(zxm::iso_nom,zxm::iso_stat,lead_pt,lead_eta,1)*zxm::varied(zxm::iso_nom,zxm::iso_stat,sublead_pt,sublead_eta,1) : 1.")
            .Define("iso_sf_stat_down", f"{enabled} ? zxm::varied(zxm::iso_nom,zxm::iso_stat,lead_pt,lead_eta,-1)*zxm::varied(zxm::iso_nom,zxm::iso_stat,sublead_pt,sublead_eta,-1) : 1.")
            .Define("trigger_pt", "lead_trigger_matched ? lead_pt : sublead_pt")
            .Define("trigger_eta", "lead_trigger_matched ? lead_eta : sublead_eta")
            .Define("trigger_sf_nom", f"{enabled} ? zxm::lookup(zxm::trig_nom,trigger_pt,trigger_eta) : 1.")
            .Define("trigger_sf_sys_up", f"{enabled} ? zxm::varied(zxm::trig_nom,zxm::trig_sys,trigger_pt,trigger_eta,1) : 1.")
            .Define("trigger_sf_sys_down", f"{enabled} ? zxm::varied(zxm::trig_nom,zxm::trig_sys,trigger_pt,trigger_eta,-1) : 1.")
            .Define("trigger_sf_stat_up", f"{enabled} ? zxm::varied(zxm::trig_nom,zxm::trig_stat,trigger_pt,trigger_eta,1) : 1.")
            .Define("trigger_sf_stat_down", f"{enabled} ? zxm::varied(zxm::trig_nom,zxm::trig_stat,trigger_pt,trigger_eta,-1) : 1.")
            .Define("reco_region", "zxm::detector_region(lead_eta,sublead_eta)")
            .Define("pileup_weight", f"zxm::pileup(reco_region,nVertices,{str(apply_pu).lower()})")
            .Define("event_weight_nominal", "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_nom*trigger_sf_nom"))
    variations = {
        "reco_sys_up": "gen_base_weight*pileup_weight*reco_sf_sys_up*iso_sf_nom*trigger_sf_nom",
        "reco_sys_down": "gen_base_weight*pileup_weight*reco_sf_sys_down*iso_sf_nom*trigger_sf_nom",
        "reco_stat_up": "gen_base_weight*pileup_weight*reco_sf_stat_up*iso_sf_nom*trigger_sf_nom",
        "reco_stat_down": "gen_base_weight*pileup_weight*reco_sf_stat_down*iso_sf_nom*trigger_sf_nom",
        "iso_sys_up": "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_sys_up*trigger_sf_nom",
        "iso_sys_down": "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_sys_down*trigger_sf_nom",
        "iso_stat_up": "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_stat_up*trigger_sf_nom",
        "iso_stat_down": "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_stat_down*trigger_sf_nom",
        "trigger_sys_up": "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_nom*trigger_sf_sys_up",
        "trigger_sys_down": "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_nom*trigger_sf_sys_down",
        "trigger_stat_up": "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_nom*trigger_sf_stat_up",
        "trigger_stat_down": "gen_base_weight*pileup_weight*reco_sf_nom*iso_sf_nom*trigger_sf_stat_down",
    }
    for name, expression in variations.items():
        node = node.Define(f"event_weight_{name}", expression)
    selected = node.Filter("selected_signal", "OS isolated signal")
    matched = selected.Filter("truth_z.fiducial && truth_reco_matched", "truth fiducial matched selected")
    region_nodes = {
        mode: {
            "A": node.Filter(f"charge_product<0 && {state}==1"),
            "B": node.Filter(f"charge_product>0 && {state}==1"),
            "C": node.Filter(f"charge_product<0 && {state}==-1"),
            "D": node.Filter(f"charge_product>0 && {state}==-1"),
        }
        for mode, state in (("both", "iso_state_both"),
                            ("at_least_one", "iso_state_at_least_one"))
    }
    reconstruction_report = matched.Report()

    truth_actions = {
        "truth_total": action(truth_total, "gen_base_weight"),
        "truth_fiducial": action(truth_fid, "gen_base_weight"),
        "truth_fiducial_pileup_weighted": action(truth_fid, "truth_eff_den_weight"),
        "selected": action(selected, "event_weight_nominal"),
        "selected_truth_matched": action(matched, "event_weight_nominal"),
    }
    variation_actions = {name: action(matched, f"event_weight_{name}") for name in variations}
    reconstruction_actions = {
        name: action(current, "gen_base_weight") for name, current in reconstruction_nodes.items()
    }
    region_actions = {
        mode: {name: action(current, "event_weight_nominal")
               for name, current in regions.items()}
        for mode, regions in region_nodes.items()
    }
    all_actions = [
        x for group in (*truth_actions.values(), *variation_actions.values(),
                        *reconstruction_actions.values(),
                        *(value for group in region_actions.values() for value in group.values()))
        for x in group
    ]
    ROOT.RDF.RunGraphs(all_actions)
    yields = {name: values(group) for name, group in truth_actions.items()}
    variation_yields = {name: values(group) for name, group in variation_actions.items()}
    weighted_reconstruction_cutflow = {
        name: values(group) for name, group in reconstruction_actions.items()
    }
    region_yields = {
        mode: {name: values(group) for name, group in regions.items()}
        for mode, regions in region_actions.items()
    }
    total = yields["truth_total"]["sum_weights"]
    fiducial = yields["truth_fiducial"]["sum_weights"]
    fiducial_eff = yields["truth_fiducial_pileup_weighted"]["sum_weights"]
    selected_matched = yields["selected_truth_matched"]["sum_weights"]
    acceptance = safe_ratio(fiducial, total)
    efficiency = safe_ratio(selected_matched, fiducial_eff)
    variation_efficiencies = {
        name: {
            "efficiency": safe_ratio(value["sum_weights"], fiducial_eff),
            "relative_to_nominal": (
                None if efficiency in (None, 0.0)
                else safe_ratio(value["sum_weights"], fiducial_eff) / efficiency - 1.0
            ),
        }
        for name, value in variation_yields.items()
    }
    result = {
        "schema_version": 2,
        "sample_kind": "mc",
        "process_name": args.process_name,
        "mc_role": args.mc_role,
        "files": [str(path.resolve()) for path in files],
        "rejected_files": rejected,
        "normalization": norm_info,
        "corrections": {"muon_momentum": apply_momentum, "muon_efficiency": apply_sf,
                        "experimental_pileup": apply_pu, "momentum_seed": seed},
        "yields": yields,
        "variation_yields": variation_yields,
        "regions_by_anti_isolation": region_yields,
        "variation_efficiencies": variation_efficiencies,
        "reconstruction_cutflow": report_rows(reconstruction_report),
        "weighted_reconstruction_cutflow": weighted_reconstruction_cutflow,
        "acceptance": acceptance,
        "efficiency": efficiency,
        "acceptance_times_efficiency": None if acceptance is None or efficiency is None else acceptance * efficiency,
    }
    if args.write_selected_skim:
        columns = ["run", "lumi", "event", "dimuon_mass", "lead_pt", "sublead_pt",
                   "lead_eta", "sublead_eta", "lead_rel_iso", "sublead_rel_iso",
                   "truth_reco_matched", "event_weight_nominal", "pileup_weight",
                   "reco_sf_nom", "iso_sf_nom", "trigger_sf_nom"]
        selected.Snapshot("Events", str(output / "mc_selected.root"), columns)
    return result


def process_data(files: list[Path], rejected: list[dict[str, str]], args: argparse.Namespace,
                 config: dict[str, Any], output: Path) -> dict[str, Any]:
    chain = core.make_chain(files)
    rdf = ROOT.RDataFrame(chain)
    if args.max_events >= 0:
        rdf = rdf.Range(min(args.max_events, int(chain.GetEntries())))
    lumi_keys_action = rdf.Define(
        "zxm_processed_lumi_key",
        "(static_cast<unsigned long long>(run) << 32) | static_cast<unsigned long long>(lumi)",
    ).Take["unsigned long long"]("zxm_processed_lumi_key")
    reco = config["physics_definition"]["reconstruction"]
    mass = reco["mass_window_gev"]
    node = (rdf.Filter("zxm::certified_lumi(run,lumi)", "certified luminosity section")
            .Filter("IsoMu24 != 0", "IsoMu24")
            .Filter("nMuons >= 2", "two reconstructed muons")
            .Define("selected_muon_indices",
                    f"zxs::select_quality(muonPt,muonEta,muonIsTight,muonIsPF,muonDxy,muonDz,{reco['muon_min_pt_gev']},{reco['muon_max_abs_eta']},true,true,0.05,0.10)")
            .Filter("selected_muon_indices.size() >= 2", "two quality muons")
            .Define("lead_index", "selected_muon_indices[0]").Define("sublead_index", "selected_muon_indices[1]")
            .Define("lead_trigger_matched", "muonIsTrigMatched[lead_index] != 0")
            .Define("sublead_trigger_matched", "muonIsTrigMatched[sublead_index] != 0")
            .Filter("lead_trigger_matched || sublead_trigger_matched", "selected pair trigger match")
            .Define("lead_pt", "static_cast<double>(muonPt[lead_index])").Define("sublead_pt", "static_cast<double>(muonPt[sublead_index])")
            .Define("lead_eta", "static_cast<double>(muonEta[lead_index])").Define("sublead_eta", "static_cast<double>(muonEta[sublead_index])")
            .Define("lead_charge", "muonCharge[lead_index]").Define("sublead_charge", "muonCharge[sublead_index]")
            .Define("dimuon_mass", "zxs::dimuon_mass(muonPt,muonEta,muonPhi,muonMass,selected_muon_indices)")
            .Filter(f"dimuon_mass>{mass[0]} && dimuon_mass<{mass[1]}", "reconstructed mass window")
            .Define("lead_rel_iso", "zxs::relative_isolation(muonPt,muonIsoCharged,muonIsoNeutral,muonIsoPhoton,muonIsoPU,lead_index)")
            .Define("sublead_rel_iso", "zxs::relative_isolation(muonPt,muonIsoCharged,muonIsoNeutral,muonIsoPhoton,muonIsoPU,sublead_index)")
            .Define("charge_product", "lead_charge*sublead_charge")
            .Define("iso_state_both", "zxs::isolation_state(lead_rel_iso,sublead_rel_iso,0.15,0.25,true)")
            .Define("iso_state_at_least_one", "zxs::isolation_state(lead_rel_iso,sublead_rel_iso,0.15,0.25,false)"))
    regions_by_mode = {
        mode: {
            "A": node.Filter(f"charge_product<0 && {state}==1"),
            "B": node.Filter(f"charge_product>0 && {state}==1"),
            "C": node.Filter(f"charge_product<0 && {state}==-1"),
            "D": node.Filter(f"charge_product>0 && {state}==-1"),
        }
        for mode, state in (("both", "iso_state_both"),
                            ("at_least_one", "iso_state_at_least_one"))
    }
    regions = regions_by_mode["both"]
    common_report = node.Report()
    actions = {
        mode: {name: current.Count() for name, current in region_nodes.items()}
        for mode, region_nodes in regions_by_mode.items()
    }
    ROOT.RDF.RunGraphs([lumi_keys_action] + [value for group in actions.values() for value in group.values()])
    counts_by_mode = {
        mode: {name: int(value.GetValue()) for name, value in group.items()}
        for mode, group in actions.items()
    }
    counts = counts_by_mode["both"]
    background = None if counts["D"] == 0 else counts["B"] * counts["C"] / counts["D"]
    lumi_pairs = sorted({(int(key) >> 32, int(key) & 0xFFFFFFFF)
                         for key in lumi_keys_action.GetValue()})
    result = {"schema_version": 2, "sample_kind": "data", "process_name": "data",
              "files": [str(path.resolve()) for path in files], "rejected_files": rejected,
              "processed_luminosity_sections": [[run, lumi] for run, lumi in lumi_pairs],
              "regions": counts, "regions_by_anti_isolation": counts_by_mode,
              "raw_abcd_background": background,
              "reconstruction_cutflow": report_rows(common_report)}
    if args.write_selected_skim:
        regions["A"].Snapshot("Events", str(output / "data_selected.root"),
                              ["run", "lumi", "event", "dimuon_mass", "lead_pt", "sublead_pt", "lead_eta", "sublead_eta"])
    return result


def main() -> int:
    args = arguments()
    config = json.loads(args.config.read_text())
    if args.truth_mass_window is not None:
        config["physics_definition"]["truth_mass_window_gev"] = args.truth_mass_window
    if args.reco_mass_window is not None:
        config["physics_definition"]["reconstruction"]["mass_window_gev"] = args.reco_mass_window
    output = args.output_dir / args.label
    output.mkdir(parents=True, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    ROOT.gInterpreter.Declare(core.CPP_HELPERS)
    ROOT.gInterpreter.Declare(build_cpp(config, args.config))
    if args.threads > 1:
        ROOT.EnableImplicitMT(args.threads)
    started = time.monotonic()
    results: dict[str, Any] = {}
    if args.sample in ("data", "both"):
        files, rejected = core.discover_input(
            args.data_dir, args.data_files_from, args.max_files, False
        )
        results["data"] = process_data(files, rejected, args, config, output)
    if args.sample in ("mc", "both"):
        files, rejected = core.discover_input(
            args.mc_dir, args.mc_files_from, args.max_files, True
        )
        results["mc"] = process_mc(files, rejected, args, config, output)

    cross_sections: dict[str, Any] = {"fiducial_pb": None, "full_pb": None}
    if "data" in results and "mc" in results:
        ndata = results["data"]["regions"]["A"]
        background = results["data"]["raw_abcd_background"]
        signal = None if background is None else ndata - background
        efficiency = results["mc"]["efficiency"]
        acceptance = results["mc"]["acceptance"]
        lumi = config["normalization"]["luminosity_pb_inverse"]
        absolute_ready = results["mc"]["normalization"]["absolute_ready"]
        if signal is not None and efficiency and acceptance and lumi and absolute_ready:
            cross_sections["fiducial_pb"] = signal / (efficiency * float(lumi))
            cross_sections["full_pb"] = signal / (acceptance * efficiency * float(lumi))
        cross_sections.update({"data_A": ndata, "provisional_abcd_background": background,
                               "provisional_signal_yield": signal})
    summary = {
        "title": "Corrected inclusive dressed-Z acceptance and efficiency study",
        "configuration": config,
        "execution": {"wall_time_seconds": time.monotonic() - started, "arguments": vars(args)},
        "results": results,
        "cross_sections": cross_sections,
        "caveats": [
            "PDG-23 ancestry is an operational MC truth definition; Z/gamma* interference is not experimentally separable.",
            "Absolute cross sections remain null until luminosity and DY normalization metadata are supplied.",
            "The raw ABCD background has no prompt subtraction or validated closure and is provisional.",
            "The highest-pT matched-muon trigger SF and coherent uncertainty variations are provisional prescriptions.",
            "The reconstructed-vertex pileup correction is experimental and disabled unless explicitly enabled."
        ],
    }
    serializable = json.loads(json.dumps(summary, default=lambda value: str(value)))
    (output / "measurement_summary.json").write_text(json.dumps(serializable, indent=2) + "\n")
    print(json.dumps({"output": str(output), "cross_sections": cross_sections,
                      "mc_acceptance": results.get("mc", {}).get("acceptance"),
                      "mc_efficiency": results.get("mc", {}).get("efficiency")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
