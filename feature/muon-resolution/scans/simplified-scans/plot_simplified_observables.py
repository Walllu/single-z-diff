#!/usr/bin/env python3
"""Build and render full-sample BB/BE/EE muon, event, and Z observables."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import simplified_scan as scan


OBSERVABLES = {
    "mass": (None, "m_{#mu#mu} [GeV]", "Z kinematics"),
    "lead_pt": ([50, 26.0, 126.0], "leading muon p_{T} [GeV]", "muons"),
    "sublead_pt": ([50, 15.0, 115.0], "subleading muon p_{T} [GeV]", "muons"),
    "lead_eta": ([50, -2.5, 2.5], "leading muon #eta", "muons"),
    "sublead_eta": ([50, -2.5, 2.5], "subleading muon #eta", "muons"),
    "lead_phi": ([36, -math.pi, math.pi], "leading muon #phi", "muons"),
    "sublead_phi": ([36, -math.pi, math.pi], "subleading muon #phi", "muons"),
    "nvertices": ([50, 0.0, 50.0], "number of reconstructed vertices", "event"),
    "z_pt": ([50, 0.0, 100.0], "Z candidate p_{T} [GeV]", "Z kinematics"),
    "z_rapidity": ([48, -2.4, 2.4], "Z candidate rapidity", "Z kinematics"),
    "z_phi": ([36, -math.pi, math.pi], "Z candidate #phi", "Z kinematics"),
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="production simplified-scan summary.json")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--mc-dir", type=Path, default=None)
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--regions", nargs="+", choices=scan.REGIONS, default=list(scan.REGIONS))
    parser.add_argument("--observables", nargs="+", choices=tuple(OBSERVABLES),
                        default=list(OBSERVABLES))
    parser.add_argument("--render-only", action="store_true")
    return parser.parse_args()


def edges(payload: dict, region: str, observable: str) -> np.ndarray:
    specification = OBSERVABLES[observable][0]
    if specification is None:
        return np.asarray(payload["binning"][region]["edges_gev"], dtype=float)
    bins, low, high = specification
    return np.linspace(low, high, int(bins) + 1)


def values(events: list[dict], seeds: list[int], parameters: dict[str, float] | None) -> dict:
    arrays = scan.prepare(events, seeds)
    replicas = len(seeds)
    if parameters is None:
        lead_factor = np.ones((len(events), replicas))
        sub_factor = np.ones((len(events), replicas))
    else:
        lead_scale = np.where(arrays["lead_barrel"], parameters["barrel_scale"],
                              parameters["endcap_scale"])
        sub_scale = np.where(arrays["sublead_barrel"], parameters["barrel_scale"],
                             parameters["endcap_scale"])
        lead_resolution = np.where(arrays["lead_barrel"], parameters["barrel_resolution"],
                                   parameters["endcap_resolution"])
        sub_resolution = np.where(arrays["sublead_barrel"], parameters["barrel_resolution"],
                                  parameters["endcap_resolution"])
        lead_factor = 1 + lead_scale[:, None] + lead_resolution[:, None] * arrays["gaussians"][:, :, 0]
        sub_factor = 1 + sub_scale[:, None] + sub_resolution[:, None] * arrays["gaussians"][:, :, 1]

    lead_pt = arrays["lead_pt"][:, None] * lead_factor
    sublead_pt = arrays["sublead_pt"][:, None] * sub_factor
    lead_eta = np.repeat(arrays["lead_eta"][:, None], replicas, axis=1)
    sublead_eta = np.repeat(arrays["sublead_eta"][:, None], replicas, axis=1)
    lead_phi = np.repeat(arrays["lead_phi"][:, None], replicas, axis=1)
    sublead_phi = np.repeat(arrays["sublead_phi"][:, None], replicas, axis=1)
    lead_mass = arrays["lead_mass"][:, None]
    sublead_mass = arrays["sublead_mass"][:, None]
    px = lead_pt * np.cos(lead_phi) + sublead_pt * np.cos(sublead_phi)
    py = lead_pt * np.sin(lead_phi) + sublead_pt * np.sin(sublead_phi)
    pz = lead_pt * np.sinh(lead_eta) + sublead_pt * np.sinh(sublead_eta)
    energy = np.sqrt((lead_pt * np.cosh(lead_eta)) ** 2 + lead_mass**2)
    energy += np.sqrt((sublead_pt * np.cosh(sublead_eta)) ** 2 + sublead_mass**2)
    mass = np.sqrt(np.maximum(0.0, energy**2 - px**2 - py**2 - pz**2))
    z_pt = np.hypot(px, py)
    z_phi = np.arctan2(py, px)
    numerator = np.maximum(energy + pz, 1e-12)
    denominator = np.maximum(energy - pz, 1e-12)
    z_rapidity = 0.5 * np.log(numerator / denominator)
    nvertices = np.repeat(arrays["nvertices"][:, None], replicas, axis=1)
    return {
        "mass": mass, "lead_pt": lead_pt, "sublead_pt": sublead_pt,
        "lead_eta": lead_eta, "sublead_eta": sublead_eta,
        "lead_phi": lead_phi, "sublead_phi": sublead_phi,
        "nvertices": nvertices, "z_pt": z_pt,
        "z_rapidity": z_rapidity, "z_phi": z_phi,
    }


def template(observable_values: np.ndarray, mass_values: np.ndarray,
             bin_edges: np.ndarray) -> scan.hf.Template:
    replicas = observable_values.shape[1]
    accepted_mass = (mass_values > 70.0) & (mass_values < 110.0)
    bin_index = np.searchsorted(bin_edges, observable_values, side="right") - 1
    valid = accepted_mass & (observable_values >= bin_edges[0]) & (observable_values < bin_edges[-1])
    bin_index[~valid] = -1
    content = np.zeros(len(bin_edges) - 1)
    variance = np.zeros_like(content)
    for b in range(len(content)):
        fraction = np.count_nonzero(bin_index == b, axis=1) / replicas
        content[b] = fraction.sum()
        variance[b] = np.square(fraction).sum()
    return scan.hf.Template(bin_edges, content, variance, float(content.sum()))


def from_root(histogram) -> scan.hf.Template:
    bins = histogram.GetNbinsX()
    bin_edges = np.asarray([histogram.GetXaxis().GetBinLowEdge(i) for i in range(1, bins + 2)])
    content = np.asarray([histogram.GetBinContent(i) for i in range(1, bins + 1)])
    variance = np.asarray([histogram.GetBinError(i) ** 2 for i in range(1, bins + 1)])
    return scan.hf.Template(bin_edges, content, variance, float(content.sum()))


def render(args: argparse.Namespace, payload: dict, output: Path) -> int:
    source = output / "observable_histograms.root"
    root_file = scan.ROOT.TFile.Open(str(source), "READ")
    if not root_file or root_file.IsZombie():
        raise SystemExit(f"Could not open {source}")
    for region in args.regions:
        for observable in args.observables:
            histograms = [root_file.Get(f"{kind}_{region}_{observable}")
                          for kind in ("data", "nominal_mc", "corrected_mc")]
            if any(not item for item in histograms):
                raise SystemExit(f"Missing {region}/{observable} histograms in {source}")
            scan.draw_validation(region, *(from_root(item) for item in histograms),
                                 payload["best_fit"], output / region,
                                 observable, OBSERVABLES[observable][1])
    root_file.Close()
    return 0


def main() -> int:
    start = time.perf_counter(); args = arguments()
    payload = json.loads(args.summary.read_text())
    output = args.output_dir or args.summary.parent
    if args.render_only:
        return render(args, payload, output)

    config = payload["configuration"]
    data_dir = args.data_dir or Path(config["data_directory"])
    mc_dir = args.mc_dir or Path(config["mc_directory"])
    max_files = config.get("max_files", -1) if args.max_files is None else args.max_files
    max_events = config.get("max_events", -1) if args.max_events is None else args.max_events
    seeds = [int(seed) for seed in config["seeds"]]
    best = payload["best_fit"]
    data_files = scan.baseline.discover(data_dir, max_files)
    mc_files = scan.baseline.discover(mc_dir, max_files)
    raw_data, data_flow, data_entries, data_processed = scan.baseline.read_events(
        data_files, max_events, False, seeds[0])
    raw_mc, mc_flow, mc_entries, mc_processed = scan.baseline.read_events(
        mc_files, max_events, True, seeds[0])
    data, mc = scan.collapse(raw_data), scan.collapse(raw_mc)
    del raw_data, raw_mc

    output.mkdir(parents=True, exist_ok=True)
    root_file = scan.ROOT.TFile(str(output / "observable_histograms.root"), "RECREATE")
    numerical = {
        "source_summary": str(args.summary.resolve()), "seeds": seeds,
        "inputs": {"data_files": len(data_files), "mc_files": len(mc_files),
                   "data_entries": data_entries, "mc_entries": mc_entries,
                   "data_processed": data_processed, "mc_processed": mc_processed},
        "cutflows": {"data": data_flow, "mc": mc_flow}, "regions": {},
    }
    for region in args.regions:
        data_values = values(data[region], [seeds[0]], None)
        nominal_values = values(mc[region], [seeds[0]], None)
        corrected_values = values(mc[region], seeds, best)
        numerical["regions"][region] = {}
        for observable in args.observables:
            bin_edges = edges(payload, region, observable)
            data_hist = template(data_values[observable], data_values["mass"], bin_edges)
            nominal_hist = template(nominal_values[observable], nominal_values["mass"], bin_edges)
            corrected_hist = template(corrected_values[observable], corrected_values["mass"], bin_edges)
            nominal_metric = scan.hf.barlow_beeston(data_hist, nominal_hist)
            corrected_metric = scan.hf.barlow_beeston(data_hist, corrected_hist)
            numerical["regions"][region][observable] = {
                "category": OBSERVABLES[observable][2],
                "data_entries": data_hist.accepted,
                "nominal_mc_entries": nominal_hist.accepted,
                "corrected_mc_entries": corrected_hist.accepted,
                "nominal": scan.hf.slim_metrics(nominal_metric),
                "corrected": scan.hf.slim_metrics(corrected_metric),
                "delta_vs_nominal": nominal_metric["objective"] - corrected_metric["objective"],
            }
            root_file.cd()
            for kind, histogram in (("data", data_hist), ("nominal_mc", nominal_hist),
                                    ("corrected_mc", corrected_hist)):
                scan.hf.root_histogram(f"{kind}_{region}_{observable}", histogram,
                                       OBSERVABLES[observable][1]).Write()
        del data_values, nominal_values, corrected_values
    root_file.Close()
    numerical["timing_seconds"] = {"histogram_production": time.perf_counter() - start}
    (output / "observable_validation_summary.json").write_text(json.dumps(numerical, indent=2) + "\n")

    # Isolate every canvas in a fresh ROOT process to avoid Cocoa backing-store
    # corruption during long batch rendering jobs.
    for region in args.regions:
        for observable in args.observables:
            command = [sys.executable, str(Path(__file__).resolve()), str(args.summary.resolve()),
                       "--output-dir", str(output), "--render-only",
                       "--regions", region, "--observables", observable]
            subprocess.run(command, check=True)
    print(json.dumps({"output": str(output), "regions": args.regions,
                      "observables": args.observables,
                      "timing_seconds": time.perf_counter() - start}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
