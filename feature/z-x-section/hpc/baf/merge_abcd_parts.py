#!/usr/bin/env python3
"""Merge BAF per-job ABCD summaries and ROOT skims without double counting."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parents[1]
sys.path.insert(0, str(ANALYSIS))
import run_z_selection as core  # noqa: E402


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parts_dir", type=Path,
                        help="Directory containing data/part_*/ and mc/part_*/")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--label", default="full_sideband_input")
    parser.add_argument("--no-hadd", action="store_true",
                        help="Merge JSON/cutflows only; do not merge ROOT skims")
    return parser.parse_args()


def load_parts(root: Path, sample: str) -> list[tuple[Path, dict[str, Any]]]:
    summaries = sorted((root / sample).glob("part_*/summary.json"))
    if not summaries:
        raise RuntimeError(f"No {sample} summaries found below {root}")
    result = []
    seen_files: dict[str, Path] = {}
    for filename in summaries:
        success = filename.parent / "SUCCESS"
        if not success.is_file():
            raise RuntimeError(f"Incomplete job output (no SUCCESS marker): {filename.parent}")
        payload = json.loads(filename.read_text())
        current = payload["samples"][sample]
        for input_file in current["files"]:
            if input_file in seen_files:
                raise RuntimeError(
                    f"Duplicate input {input_file} in {seen_files[input_file]} and {filename}; "
                    "use a fresh run label or remove the duplicate retry output"
                )
            seen_files[input_file] = filename
        result.append((filename, payload))
    return result


def summed_cutflow(parts: list[tuple[Path, dict[str, Any]]], sample: str) -> list[dict[str, Any]]:
    order = [row["step"] for row in parts[0][1]["samples"][sample]["cutflow"]]
    totals = {step: {"events": 0, "weighted_yield": 0.0, "sum_weights_squared": 0.0}
              for step in order}
    for _, payload in parts:
        rows = payload["samples"][sample]["cutflow"]
        if [row["step"] for row in rows] != order:
            raise RuntimeError(f"Incompatible cutflow steps in {sample} parts")
        for row in rows:
            target = totals[row["step"]]
            for field in target:
                target[field] += row[field]
    initial = totals[order[0]]
    previous: dict[str, Any] | None = None
    result = []
    for step in order:
        value = totals[step]
        result.append({
            "sample": sample, "step": step, **value,
            "step_efficiency": value["events"] / previous["events"] if previous and previous["events"] else (1.0 if value["events"] else 0.0),
            "cumulative_efficiency": value["events"] / initial["events"] if initial["events"] else 0.0,
            "weighted_step_efficiency": value["weighted_yield"] / previous["weighted_yield"] if previous and previous["weighted_yield"] else (1.0 if value["weighted_yield"] else 0.0),
            "weighted_cumulative_efficiency": value["weighted_yield"] / initial["weighted_yield"] if initial["weighted_yield"] else 0.0,
        })
        previous = value
    return result


def aggregate(parts: list[tuple[Path, dict[str, Any]]], sample: str) -> dict[str, Any]:
    first = parts[0][1]["samples"][sample]
    cutflow = summed_cutflow(parts, sample)
    regions: dict[str, dict[str, Any]] = {}
    for name, definition in core.ABCD_REGIONS.items():
        values = [payload["samples"][sample]["regions"][name] for _, payload in parts]
        regions[name] = {
            **definition,
            "events": sum(item["events"] for item in values),
            "weighted_yield": sum(item["weighted_yield"] for item in values),
            "sum_weights_squared": sum(item["sum_weights_squared"] for item in values),
        }
    files = [item for _, payload in parts for item in payload["samples"][sample]["files"]]
    rejected = [item for _, payload in parts
                for item in payload["samples"][sample]["rejected_files"]]
    return {
        "sample": sample,
        "is_mc": sample == "mc",
        "files": files,
        "files_used": len(files),
        "rejected_files": rejected,
        "tree_entries": sum(payload["samples"][sample]["tree_entries"] for _, payload in parts),
        "entries_requested": -1,
        "cutflow": cutflow,
        "regions": regions,
        "raw_abcd_nonprompt_estimate_in_sr": core.abcd_estimate(regions),
        "skim": None,
        "merged_parts": len(parts),
    }


def main() -> int:
    args = arguments()
    output = args.output_dir / args.label
    output.mkdir(parents=True, exist_ok=True)
    part_sets = {sample: load_parts(args.parts_dir, sample) for sample in ("data", "mc")}
    configurations = [
        payload["configuration"]
        for parts in part_sets.values()
        for _, payload in parts
    ]
    reference = configurations[0]
    if any(configuration != reference for configuration in configurations[1:]):
        raise RuntimeError("Per-job selections differ; refusing to merge incompatible outputs")

    all_cutflow: list[dict[str, Any]] = []
    samples: dict[str, Any] = {}
    for sample in ("data", "mc"):
        parts = part_sets[sample]
        samples[sample] = aggregate(parts, sample)
        all_cutflow.extend(samples[sample]["cutflow"])
        if not args.no_hadd:
            inputs = [payload["samples"][sample]["skim"] for _, payload in parts]
            merged = output / f"{sample}_abcd_skim.root"
            subprocess.run(["hadd", "-f", str(merged), *inputs], check=True)
            samples[sample]["skim"] = str(merged.resolve())
    summary = {
        "title": "Merged BAF inclusive Z dimuon reconstruction and ABCD selection",
        "configuration": reference,
        "execution": {"source_parts_directory": str(args.parts_dir.resolve())},
        "samples": samples,
        "caveats": [
            "The merge rejects duplicate input paths to prevent retry double counting.",
            "Prompt subtraction and closure are evaluated in the downstream diagnostic step.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    core.write_cutflow_csv(output / "cutflow.csv", all_cutflow)
    print(f"Merged outputs written to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
