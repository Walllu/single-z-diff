#!/usr/bin/env python3
"""High-fidelity signed-eta muon scale/resolution template scan.

This scan uses event-keyed Gaussian replicas, an analytic single-source
Barlow--Beeston profile likelihood, adaptive resonance-aware mass bins, and a
deterministic train/evaluation split. It intentionally leaves the earlier scan
products unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import sys
import time
from array import array
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCANS = HERE.parent
sys.path.insert(0, str(SCANS))

import ROOT  # noqa: E402
import scan_zmumu_scale_resolution as baseline  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)


# -----------------------------------------------------------------------------
# User-facing configuration
# -----------------------------------------------------------------------------
RNG_SEEDS = [314159, 271828, 161803, 141421, 173205, 223607, 244949, 264575]
DEFAULT_REPLICAS = 3

GRID_PROFILES = {
    "production": {
        "scale_shifts": [round(-0.005 + 0.0005 * i, 7) for i in range(21)],
        "resolution_smearings": [round(0.0025 * i, 7) for i in range(21)],
    },
    "test": {
        "scale_shifts": [round(-0.005 + 0.001 * i, 7) for i in range(11)],
        "resolution_smearings": [round(0.005 * i, 7) for i in range(11)],
    },
    "smoke": {
        "scale_shifts": [-0.004, -0.002, 0.0, 0.002, 0.004],
        "resolution_smearings": [0.0, 0.01, 0.02, 0.03, 0.04],
    },
}

SPLIT_CONFIGURATION = {
    "training_fraction": 0.80,
    "split_seed": 8675309,
}

DYNAMIC_BINNING = {
    "range_gev": [70.0, 110.0],
    "segments": [
        # low, high, finest permitted width; low-stat adjacent bins are merged
        [70.0, 82.0, 2.0],
        [82.0, 84.0, 1.0],
        [84.0, 98.0, 0.5],
        [98.0, 100.0, 1.0],
        [100.0, 110.0, 2.0],
    ],
    "minimum_training_data_per_bin": 10,
    "minimum_training_mc_per_bin": 20,
}


@dataclass
class Template:
    edges: np.ndarray
    content: np.ndarray
    variance: np.ndarray
    accepted: float


def arguments() -> argparse.Namespace:
    raw = baseline.default_raw_data_directory()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=raw / "2016H")
    parser.add_argument("--mc-dir", type=Path,
                        default=raw / "ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8")
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--label", default=None)
    parser.add_argument("--regions", nargs="+", default=["barrel__barrel"])
    parser.add_argument("--include-excluded-regions", action="store_true")
    parser.add_argument("--grid-profile", choices=tuple(GRID_PROFILES), default="production")
    parser.add_argument("--replicas", type=int, default=DEFAULT_REPLICAS,
                        help="Number of event-keyed seed replicas averaged per grid point")
    parser.add_argument("--seeds", nargs="+", type=int, default=None,
                        help="Explicit seed list; overrides the first --replicas seeds")
    parser.add_argument("--train-fraction", type=float,
                        default=SPLIT_CONFIGURATION["training_fraction"])
    parser.add_argument("--split-seed", type=int, default=SPLIT_CONFIGURATION["split_seed"])
    parser.add_argument("--minimum-data-per-bin", type=int,
                        default=DYNAMIC_BINNING["minimum_training_data_per_bin"])
    parser.add_argument("--minimum-mc-per-bin", type=int,
                        default=DYNAMIC_BINNING["minimum_training_mc_per_bin"])
    parser.add_argument("--max-files", type=int, default=-1)
    parser.add_argument("--max-events", type=int, default=-1)
    return parser.parse_args()


def stable_uniform(seed: int, sample: str, event_key: str) -> float:
    payload = f"{seed}:{sample}:{event_key}".encode()
    value = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
    return value / 2**64


def keyed_gaussians(seed: int, event_key: str) -> tuple[float, float]:
    rng = random.Random(f"{seed}:{event_key}")
    return rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)


def split_events(events: dict[str, list[dict]], fraction: float, split_seed: int,
                 sample: str) -> tuple[dict[str, list[dict]], dict[str, list[dict]]]:
    train = {region: [] for region in baseline.ETA_PAIR_REGIONS}
    evaluation = {region: [] for region in baseline.ETA_PAIR_REGIONS}
    for region, values in events.items():
        for event in values:
            target = train if stable_uniform(split_seed, sample, event["event_key"]) < fraction else evaluation
            target[region].append(event)
    return train, evaluation


def base_binning_edges() -> np.ndarray:
    edges: list[float] = []
    for low, high, width in DYNAMIC_BINNING["segments"]:
        segment = [round(low + width * i, 8) for i in range(round((high - low) / width) + 1)]
        if edges and math.isclose(edges[-1], segment[0]):
            segment = segment[1:]
        edges.extend(segment)
    return np.asarray(edges, dtype=float)


def nominal_masses(events: list[dict]) -> np.ndarray:
    return np.asarray([baseline.dimuon_mass(e["lead"], e["sublead"], 1.0, 1.0)
                       for e in events], dtype=float)


def adaptive_edges(data_events: list[dict], mc_events: list[dict], min_data: int,
                   min_mc: int) -> tuple[np.ndarray, dict]:
    """Merge resonance-aware base bins until both training samples meet minima."""
    base = base_binning_edges()
    data_counts, _ = np.histogram(nominal_masses(data_events), bins=base)
    mc_counts, _ = np.histogram(nominal_masses(mc_events), bins=base)
    output = [float(base[0])]
    accumulated_data = accumulated_mc = 0
    closed_at_high = False
    for index, (data_count, mc_count) in enumerate(zip(data_counts, mc_counts)):
        accumulated_data += int(data_count)
        accumulated_mc += int(mc_count)
        upper = float(base[index + 1])
        if accumulated_data >= min_data and accumulated_mc >= min_mc:
            output.append(upper)
            accumulated_data = accumulated_mc = 0
            closed_at_high = index == len(data_counts) - 1
    if not closed_at_high:
        # Merge a deficient final tail with the preceding accepted bin.
        if len(output) > 1:
            output.pop()
        output.append(float(base[-1]))
    edges = np.asarray(output, dtype=float)
    final_data, _ = np.histogram(nominal_masses(data_events), bins=edges)
    final_mc, _ = np.histogram(nominal_masses(mc_events), bins=edges)
    details = {
        "base_edges_gev": base.tolist(), "edges_gev": edges.tolist(),
        "minimum_training_data_per_bin": min_data,
        "minimum_training_mc_per_bin": min_mc,
        "training_data_counts": final_data.astype(int).tolist(),
        "training_mc_counts": final_mc.astype(int).tolist(),
    }
    return edges, details


def prepare_arrays(events: list[dict], seeds: list[int]) -> dict[str, np.ndarray]:
    n = len(events)
    result: dict[str, np.ndarray] = {}
    for position, name in enumerate(("lead", "sublead")):
        for component, index in (("pt", 0), ("eta", 1), ("phi", 2), ("mass", 3)):
            result[f"{name}_{component}"] = np.fromiter(
                (event[name][index] for event in events), dtype=float, count=n)
    result["nvertices"] = np.fromiter((event["nvertices"] for event in events), dtype=float, count=n)
    gaussian = np.empty((n, len(seeds), 2), dtype=float)
    for i, event in enumerate(events):
        for j, seed in enumerate(seeds):
            gaussian[i, j] = keyed_gaussians(seed, event["event_key"])
    result["gaussians"] = gaussian
    return result


def mass_matrix(values: dict[str, np.ndarray], scale: float, resolution: float,
                replica_indices: list[int] | None = None) -> np.ndarray:
    gaussians = values["gaussians"]
    indices = replica_indices if replica_indices is not None else list(range(gaussians.shape[1]))
    factors1 = 1.0 + scale + resolution * gaussians[:, indices, 0]
    factors2 = 1.0 + scale + resolution * gaussians[:, indices, 1]
    pt1 = values["lead_pt"][:, None] * factors1
    pt2 = values["sublead_pt"][:, None] * factors2
    eta1, eta2 = values["lead_eta"][:, None], values["sublead_eta"][:, None]
    phi1, phi2 = values["lead_phi"][:, None], values["sublead_phi"][:, None]
    m1, m2 = values["lead_mass"][:, None], values["sublead_mass"][:, None]
    px = pt1 * np.cos(phi1) + pt2 * np.cos(phi2)
    py = pt1 * np.sin(phi1) + pt2 * np.sin(phi2)
    pz = pt1 * np.sinh(eta1) + pt2 * np.sinh(eta2)
    energy = np.sqrt((pt1 * np.cosh(eta1)) ** 2 + m1 ** 2)
    energy += np.sqrt((pt2 * np.cosh(eta2)) ** 2 + m2 ** 2)
    return np.sqrt(np.maximum(0.0, energy**2 - px**2 - py**2 - pz**2))


def data_template(events: list[dict], edges: np.ndarray) -> Template:
    content, _ = np.histogram(nominal_masses(events), bins=edges)
    values = content.astype(float)
    return Template(edges, values, values.copy(), float(values.sum()))


def mc_template(values: dict[str, np.ndarray], edges: np.ndarray, scale: float,
                resolution: float, replica_indices: list[int] | None = None) -> Template:
    masses = mass_matrix(values, scale, resolution, replica_indices)
    replicas = masses.shape[1]
    bin_index = np.searchsorted(edges, masses, side="right") - 1
    valid = (masses > edges[0]) & (masses < edges[-1])
    bin_index[~valid] = -1
    content = np.zeros(len(edges) - 1, dtype=float)
    variance = np.zeros_like(content)
    for b in range(len(content)):
        fraction = np.count_nonzero(bin_index == b, axis=1).astype(float) / replicas
        content[b] = fraction.sum()
        # Replicas integrate the smearing kernel; physical MC events remain the
        # independent statistical units for Sumw2/Barlow--Beeston purposes.
        variance[b] = np.square(fraction).sum()
    return Template(edges, content, variance, float(content.sum()))


def poisson_deviance(observed: float, expected: float) -> float:
    if observed <= 0:
        return 2.0 * expected
    expected = max(expected, 1e-12)
    return 2.0 * (expected - observed + observed * math.log(observed / expected))


def barlow_beeston(data: Template, mc: Template) -> dict:
    """Profile one MC-stat nuisance per bin for a single simulated source.

    Fractional averaged-replica templates are represented by effective counts
    m_eff=m^2/sumw2. With one source this analytic profile is equivalent to the
    full per-bin Barlow--Beeston construction.
    """
    total_data, total_mc = data.content.sum(), mc.content.sum()
    alpha = total_data / total_mc if total_mc > 0 else 0.0
    total = 0.0
    contributions: list[float] = []
    pulls: list[float] = []
    empty_expected_bins = 0
    for observed, simulated, variance in zip(data.content, mc.content, mc.variance):
        expected = alpha * simulated
        if simulated <= 0 or variance <= 0 or alpha <= 0:
            contribution = poisson_deviance(float(observed), max(float(expected), 1e-12))
            if observed > 0 and expected <= 0:
                empty_expected_bins += 1
        else:
            effective_mc = simulated * simulated / variance
            tau = effective_mc / expected
            profiled_data_mean = (observed + effective_mc) / (1.0 + tau)
            contribution = poisson_deviance(float(observed), float(profiled_data_mean))
            contribution += poisson_deviance(float(effective_mc), float(tau * profiled_data_mean))
        total += contribution
        contributions.append(float(contribution))
        sign = 1.0 if observed >= expected else -1.0
        pulls.append(sign * math.sqrt(max(0.0, contribution)))
    return {
        "objective": float(total), "ndof": max(0, len(data.content) - 1),
        "normalization": float(alpha), "empty_expected_bins": empty_expected_bins,
        "bin_contributions": contributions, "signed_profile_pulls": pulls,
    }


def root_histogram(name: str, template: Template, axis_title: str = "m_{#mu#mu} [GeV]"):
    hist = ROOT.TH1D(name, f";{axis_title};Events", len(template.edges) - 1,
                     array("d", template.edges.tolist()))
    hist.SetDirectory(0)
    for index, (content, variance) in enumerate(zip(template.content, template.variance), 1):
        hist.SetBinContent(index, float(content))
        hist.SetBinError(index, math.sqrt(max(0.0, float(variance))))
    return hist


def draw_scan(region: str, results: list[dict], best: dict, grid: dict, output: Path) -> None:
    scales, resolutions = grid["scale_shifts"], grid["resolution_smearings"]
    dx = scales[1] - scales[0] if len(scales) > 1 else 0.001
    dy = resolutions[1] - resolutions[0] if len(resolutions) > 1 else 0.0025
    surface = ROOT.TH2D("hf_surface_" + region,
                        f"{baseline.abbreviated_region(region)};scale shift;added resolution;#Delta BB deviance",
                        len(scales), min(scales) - dx / 2, max(scales) + dx / 2,
                        len(resolutions), min(resolutions) - dy / 2, max(resolutions) + dy / 2)
    surface.SetDirectory(0)
    for item in results:
        surface.Fill(item["scale"], item["resolution"], item["objective"] - best["objective"])
    surface.SetMinimum(0); surface.SetMaximum(30)
    canvas = ROOT.TCanvas("hf_canvas_" + region, "", 850, 700)
    canvas.SetRightMargin(0.17); surface.Draw("COLZ")
    contours = surface.Clone("hf_contours_" + region)
    contours.SetContour(3, array("d", [2.30, 6.18, 11.83]))
    contours.SetLineColor(ROOT.kBlack); contours.SetLineWidth(2); contours.Draw("CONT3 SAME")
    marker = ROOT.TMarker(best["scale"], best["resolution"], 20)
    marker.SetMarkerColor(ROOT.kRed + 1); marker.SetMarkerSize(1.5); marker.Draw("SAME")
    canvas.SaveAs(str(output / f"{region}_hf_scan.png"))
    canvas.SaveAs(str(output / f"{region}_hf_scan.pdf"))


def slim_metrics(metrics: dict) -> dict:
    return {key: value for key, value in metrics.items()
            if key not in {"bin_contributions", "signed_profile_pulls"}}


def scan_region(region: str, train_data: list[dict], train_mc: list[dict],
                eval_data: list[dict], eval_mc: list[dict], seeds: list[int],
                grid: dict, args: argparse.Namespace, output: Path) -> tuple[dict, list]:
    edges, binning = adaptive_edges(train_data, train_mc, args.minimum_data_per_bin,
                                    args.minimum_mc_per_bin)
    train_data_hist = data_template(train_data, edges)
    eval_data_hist = data_template(eval_data, edges)
    train_arrays = prepare_arrays(train_mc, seeds)
    eval_arrays = prepare_arrays(eval_mc, seeds)

    zero_train_mc = mc_template(train_arrays, edges, 0.0, 0.0)
    zero_eval_mc = mc_template(eval_arrays, edges, 0.0, 0.0)
    zero_train = barlow_beeston(train_data_hist, zero_train_mc)
    zero_eval = barlow_beeston(eval_data_hist, zero_eval_mc)

    results: list[dict] = []
    best: dict | None = None
    for scale in grid["scale_shifts"]:
        for resolution in grid["resolution_smearings"]:
            template = mc_template(train_arrays, edges, scale, resolution)
            metrics = barlow_beeston(train_data_hist, template)
            item = {"scale": scale, "resolution": resolution,
                    **slim_metrics(metrics), "mc_accepted": template.accepted}
            results.append(item)
            if best is None or item["objective"] < best["objective"]:
                best = item
    assert best is not None

    best_train_mc = mc_template(train_arrays, edges, best["scale"], best["resolution"])
    best_eval_mc = mc_template(eval_arrays, edges, best["scale"], best["resolution"])
    evaluation = barlow_beeston(eval_data_hist, best_eval_mc)
    evaluation_by_seed = {}
    for index, seed in enumerate(seeds):
        seed_mc = mc_template(eval_arrays, edges, best["scale"], best["resolution"], [index])
        metrics = barlow_beeston(eval_data_hist, seed_mc)
        evaluation_by_seed[str(seed)] = {
            **slim_metrics(metrics), "delta_vs_zero": zero_eval["objective"] - metrics["objective"],
            "mc_accepted": seed_mc.accepted,
        }

    result = {
        "counts": {"training_data": len(train_data), "training_mc": len(train_mc),
                   "evaluation_data": len(eval_data), "evaluation_mc": len(eval_mc)},
        "binning": binning,
        "zero_correction": {"training": slim_metrics(zero_train),
                            "evaluation": slim_metrics(zero_eval)},
        "best": best,
        "training_delta_vs_zero": zero_train["objective"] - best["objective"],
        "evaluation_at_best": {**slim_metrics(evaluation),
                               "delta_vs_zero": zero_eval["objective"] - evaluation["objective"]},
        "evaluation_by_seed": evaluation_by_seed,
        "grid_results": results,
    }
    draw_scan(region, results, best, grid, output)
    histograms = [
        root_histogram(f"train_data_{region}_mass", train_data_hist),
        root_histogram(f"train_zero_mc_{region}_mass", zero_train_mc),
        root_histogram(f"train_best_mc_{region}_mass", best_train_mc),
        root_histogram(f"eval_data_{region}_mass", eval_data_hist),
        root_histogram(f"eval_zero_mc_{region}_mass", zero_eval_mc),
        root_histogram(f"eval_best_mc_{region}_mass", best_eval_mc),
    ]
    return result, histograms


def configuration_label(args: argparse.Namespace, seeds: list[int]) -> str:
    if args.label:
        return args.label
    fraction = round(100 * args.train_fraction)
    return f"hf_{args.grid_profile}_replicas{len(seeds)}_train{fraction}"


def main() -> int:
    start = time.perf_counter()
    args = arguments()
    if not 0 < args.train_fraction < 1:
        raise SystemExit("--train-fraction must lie strictly between zero and one")
    if args.replicas < 1:
        raise SystemExit("--replicas must be positive")
    seeds = args.seeds if args.seeds else RNG_SEEDS[:args.replicas]
    if not seeds:
        raise SystemExit("At least one RNG seed is required")
    requested = baseline.ETA_PAIR_REGIONS if args.regions == ["all"] else args.regions
    unknown = sorted(set(requested) - set(baseline.ETA_PAIR_REGIONS))
    if unknown:
        raise SystemExit(f"Unknown regions: {unknown}")
    excluded = [r for r in requested if baseline.abbreviated_region(r)
                in baseline.EXCLUDED_ETA_PAIR_ABBREVIATIONS]
    regions = requested if args.include_excluded_regions else [r for r in requested if r not in excluded]
    grid = GRID_PROFILES[args.grid_profile]
    output = args.output_dir / configuration_label(args, seeds)
    output.mkdir(parents=True, exist_ok=True)

    extraction_start = time.perf_counter()
    data_files = baseline.discover(args.data_dir, args.max_files)
    mc_files = baseline.discover(args.mc_dir, args.max_files)
    data_events, data_flow, data_total, data_processed = baseline.read_events(
        data_files, args.max_events, False, seeds[0])
    # The baseline reader categorizes all 25 regions. Release unrequested event
    # records before loading MC so directed high-fidelity runs have a bounded
    # memory footprint.
    data_events = {region: data_events[region] if region in regions else []
                   for region in baseline.ETA_PAIR_REGIONS}
    mc_events, mc_flow, mc_total, mc_processed = baseline.read_events(
        mc_files, args.max_events, True, seeds[0])
    mc_events = {region: mc_events[region] if region in regions else []
                 for region in baseline.ETA_PAIR_REGIONS}
    train_data, eval_data = split_events(data_events, args.train_fraction, args.split_seed, "data")
    train_mc, eval_mc = split_events(mc_events, args.train_fraction, args.split_seed, "mc")
    extraction_seconds = time.perf_counter() - extraction_start

    summary = {
        "schema_version": 1,
        "configuration": {
            "method": "single-source analytic Barlow-Beeston profile",
            "event_keyed_rng": True, "seeds": seeds, "replicas_averaged": len(seeds),
            "grid_profile": args.grid_profile, "scan_grid": grid,
            "train_fraction": args.train_fraction, "evaluation_fraction": 1.0 - args.train_fraction,
            "split_seed": args.split_seed, "split_method": "SHA256 threshold by sample and event key",
            "selection": baseline.SELECTION, "dynamic_binning": DYNAMIC_BINNING,
            "minimum_data_per_bin": args.minimum_data_per_bin,
            "minimum_mc_per_bin": args.minimum_mc_per_bin,
            "regions": regions, "excluded_from_run": [baseline.abbreviated_region(r) for r in excluded],
            "data_directory": str(args.data_dir.resolve()), "mc_directory": str(args.mc_dir.resolve()),
            "max_files": args.max_files, "max_events": args.max_events,
            "output_directory": str(output.resolve()),
        },
        "inputs": {"data_files": len(data_files), "mc_files": len(mc_files),
                   "data_entries": data_total, "mc_entries": mc_total,
                   "data_processed": data_processed, "mc_processed": mc_processed},
        "cutflows": {"data": data_flow, "mc": mc_flow},
        "regions": {},
    }
    root_file = ROOT.TFile(str(output / "high_fidelity_histograms.root"), "RECREATE")
    scan_start = time.perf_counter()
    for region in regions:
        region_output = output / baseline.abbreviated_region(region)
        region_output.mkdir(parents=True, exist_ok=True)
        result, histograms = scan_region(region, train_data[region], train_mc[region],
                                         eval_data[region], eval_mc[region], seeds, grid,
                                         args, region_output)
        summary["regions"][region] = result
        root_file.cd()
        for histogram in histograms:
            histogram.Write()
    root_file.Close()
    summary["timing_seconds"] = {
        "event_extraction_and_split": extraction_seconds,
        "regional_scans": time.perf_counter() - scan_start,
        "total": time.perf_counter() - start,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    brief = {
        "output": str(output), "seeds": seeds, "regions": len(regions),
        "timing_seconds": summary["timing_seconds"],
        "results": {baseline.abbreviated_region(r): {
            "best": v["best"], "train_delta": v["training_delta_vs_zero"],
            "evaluation_delta": v["evaluation_at_best"]["delta_vs_zero"],
        } for r, v in summary["regions"].items()},
    }
    print(json.dumps(brief, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
