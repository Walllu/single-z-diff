#!/usr/bin/env python3
"""Merge additive Z cross-section inputs and form only global MC ratios."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any


TRIPLE = ("events", "sum_weights", "sum_weights_squared")


def add_triple(target: dict[str, float | int], source: dict[str, Any]) -> None:
    target["events"] = int(target.get("events", 0)) + int(source["events"])
    target["sum_weights"] = float(target.get("sum_weights", 0.0)) + float(source["sum_weights"])
    target["sum_weights_squared"] = float(target.get("sum_weights_squared", 0.0)) + float(source["sum_weights_squared"])


def merge_triple_maps(payloads: list[dict[str, Any]], key: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for payload in payloads:
        for name, value in payload.get(key, {}).items():
            add_triple(result.setdefault(name, {}), value)
    return result


def merge_nested_triples(payloads: list[dict[str, Any]], key: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for payload in payloads:
        for outer, values in payload.get(key, {}).items():
            for name, value in values.items():
                add_triple(result.setdefault(outer, {}).setdefault(name, {}), value)
    return result


def ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def compress_lumis(pairs: set[tuple[int, int]]) -> dict[str, list[list[int]]]:
    by_run: dict[int, list[int]] = {}
    for run, lumi in sorted(pairs):
        by_run.setdefault(run, []).append(lumi)
    result: dict[str, list[list[int]]] = {}
    for run, lumis in by_run.items():
        ranges: list[list[int]] = []
        start = previous = lumis[0]
        for current in lumis[1:]:
            if current != previous + 1:
                ranges.append([start, previous])
                start = current
            previous = current
        ranges.append([start, previous])
        result[str(run)] = ranges
    return result


def merge_cutflow(payloads: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ordered: list[str] = []
    totals: dict[str, list[int]] = {}
    for payload in payloads:
        for row in payload.get("reconstruction_cutflow", []):
            step = str(row["step"])
            if step not in totals:
                ordered.append(step)
                totals[step] = [0, 0]
            totals[step][0] += int(row["events"])
            totals[step][1] += int(row["events_before"])
    return [{"step": step, "events": totals[step][0], "events_before": totals[step][1],
             "step_efficiency_percent": (100.0 * totals[step][0] / totals[step][1]
                                           if totals[step][1] else 0.0)}
            for step in ordered]


def merge_mc(payloads: list[dict[str, Any]], role: str, name: str) -> dict[str, Any]:
    yields = merge_triple_maps(payloads, "yields")
    variations = merge_triple_maps(payloads, "variation_yields")
    regions = merge_nested_triples(payloads, "regions_by_anti_isolation")
    weighted_cutflow = merge_triple_maps(payloads, "weighted_reconstruction_cutflow")
    processed_sumw = sum(float(p["normalization"]["processed_sum_generator_weights"])
                         for p in payloads)
    total = float(yields.get("truth_total", {}).get("sum_weights", 0.0))
    fiducial = float(yields.get("truth_fiducial", {}).get("sum_weights", 0.0))
    fiducial_eff = float(yields.get("truth_fiducial_pileup_weighted", {}).get("sum_weights", 0.0))
    selected_matched = float(yields.get("selected_truth_matched", {}).get("sum_weights", 0.0))
    acceptance = ratio(fiducial, total)
    efficiency = ratio(selected_matched, fiducial_eff)
    variation_efficiencies = {
        variation: {
            "efficiency": ratio(float(value["sum_weights"]), fiducial_eff),
            "relative_to_nominal": (
                None if not efficiency else float(value["sum_weights"]) / fiducial_eff / efficiency - 1.0
            ),
        }
        for variation, value in variations.items()
    }
    return {
        "process_name": name,
        "mc_role": role,
        "files": sorted(path for payload in payloads for path in payload["files"]),
        "processed_sum_generator_weights": processed_sumw,
        "yields": yields,
        "variation_yields": variations,
        "variation_efficiencies": variation_efficiencies,
        "regions_by_anti_isolation": regions,
        "reconstruction_cutflow": merge_cutflow(payloads),
        "weighted_reconstruction_cutflow": weighted_cutflow,
        "acceptance": acceptance,
        "efficiency": efficiency,
        "acceptance_times_efficiency": (
            None if acceptance is None or efficiency is None else acceptance * efficiency
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parts_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    discovered = sorted(args.parts_dir.glob("**/part_*/measurement_summary.json"))
    if not discovered:
        raise RuntimeError(f"No measurement summaries below {args.parts_dir}")
    incomplete = [path for path in discovered if not (path.parent / "SUCCESS").is_file()]
    if incomplete:
        raise RuntimeError(f"Measurement summaries without SUCCESS marker: {incomplete[:5]}")
    summaries = discovered
    documents = [json.loads(path.read_text()) for path in summaries]
    configuration = documents[0]["configuration"]
    if any(doc["configuration"] != configuration for doc in documents[1:]):
        raise RuntimeError("Refusing to merge measurement parts with different configurations")

    groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    all_files: list[str] = []
    for document in documents:
        for result_key, payload in document["results"].items():
            kind = payload.get("sample_kind", result_key)
            role = payload.get("mc_role", "data" if kind == "data" else "signal")
            name = payload.get("process_name", "data" if kind == "data" else "z_to_mumu")
            groups.setdefault((kind, role, name), []).append(payload)
            all_files.extend(payload["files"])
            if kind == "mc" and not payload["normalization"].get("deferred_to_final_stage", False):
                raise RuntimeError("MC part was normalized before merging; rerun with --defer-normalization")
    duplicates = sorted(path for path, count in Counter(all_files).items() if count > 1)
    if duplicates:
        raise RuntimeError(f"Duplicate input files across measurement parts: {duplicates[:5]}")

    data_groups = [payloads for (kind, _, _), payloads in groups.items() if kind == "data"]
    if len(data_groups) > 1:
        raise RuntimeError("Found multiple independently named data groups")
    data = None
    processed_lumis: set[tuple[int, int]] = set()
    if data_groups:
        payloads = data_groups[0]
        regions: dict[str, dict[str, int]] = {}
        for payload in payloads:
            for mode, counts in payload.get("regions_by_anti_isolation", {"both": payload["regions"]}).items():
                for region, count in counts.items():
                    regions.setdefault(mode, {}).setdefault(region, 0)
                    regions[mode][region] += int(count)
            processed_lumis.update((int(run), int(lumi))
                                   for run, lumi in payload.get("processed_luminosity_sections", []))
        data = {
            "files": sorted(path for payload in payloads for path in payload["files"]),
            "processed_luminosity_sections": [[run, lumi] for run, lumi in sorted(processed_lumis)],
            "regions_by_anti_isolation": regions,
            "reconstruction_cutflow": merge_cutflow(payloads),
        }

    signal = None
    backgrounds: dict[str, Any] = {}
    for (kind, role, name), payloads in groups.items():
        if kind != "mc":
            continue
        merged = merge_mc(payloads, role, name)
        if role == "signal":
            if signal is not None:
                raise RuntimeError("Only one signal process is supported")
            signal = merged
        else:
            if name in backgrounds:
                raise RuntimeError(f"Duplicate background process name: {name}")
            backgrounds[name] = merged

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    lumi_json = compress_lumis(processed_lumis) if processed_lumis else {}
    (output / "processed_lumis.json").write_text(json.dumps(lumi_json, indent=2) + "\n")
    result = {
        "schema_version": 1,
        "title": "Merged additive inputs for the inclusive Z cross-section",
        "configuration": configuration,
        "provenance": {"part_summaries": [str(path.resolve()) for path in summaries],
                       "input_files": len(all_files)},
        "data": data,
        "signal": signal,
        "backgrounds": backgrounds,
        "luminosity": {
            "observed_unique_sections": len(processed_lumis),
            "processed_lumis_json": str((output / "processed_lumis.json").resolve()),
            "integrated_luminosity_pb_inverse": None,
            "note": "An observed run/lumisection inventory is not a luminosity measurement; intersect with certification and evaluate with approved luminosity data.",
        },
    }
    (output / "measurement_inputs.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(output), "parts": len(summaries),
                      "data_files": len(data["files"]) if data else 0,
                      "signal_files": len(signal["files"]) if signal else 0,
                      "backgrounds": sorted(backgrounds),
                      "acceptance": signal.get("acceptance") if signal else None,
                      "efficiency": signal.get("efficiency") if signal else None}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
