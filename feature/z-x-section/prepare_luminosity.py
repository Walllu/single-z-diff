#!/usr/bin/env python3
"""Intersect processed/certified lumisections and optionally parse BRIL output."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any


def expand(mask: dict[str, list[list[int]]]) -> set[tuple[int, int]]:
    return {(int(run), lumi) for run, ranges in mask.items()
            for start, stop in ranges for lumi in range(int(start), int(stop) + 1)}


def compress(pairs: set[tuple[int, int]]) -> dict[str, list[list[int]]]:
    result: dict[str, list[list[int]]] = {}
    for run in sorted({run for run, _ in pairs}):
        lumis = sorted(lumi for candidate, lumi in pairs if candidate == run)
        ranges: list[list[int]] = []
        if lumis:
            start = previous = lumis[0]
            for current in lumis[1:]:
                if current != previous + 1:
                    ranges.append([start, previous])
                    start = current
                previous = current
            ranges.append([start, previous])
        result[str(run)] = ranges
    return result


def parse_brilcalc(filename: Path) -> tuple[float, dict[str, Any]]:
    lines = [line.strip().lstrip("#").strip() for line in filename.read_text().splitlines()
             if line.strip()]
    for index, line in enumerate(lines):
        fields = [field.strip() for field in next(csv.reader([line]))]
        recorded = next((position for position, field in enumerate(fields)
                         if field.lower().startswith("totrecorded(")), None)
        if recorded is None:
            continue
        match = re.search(r"\((/[^)]+)\)", fields[recorded])
        if match is None:
            raise ValueError(f"Cannot determine BRIL luminosity unit from {fields[recorded]!r}")
        for candidate in lines[index + 1:]:
            values = [value.strip() for value in next(csv.reader([candidate]))]
            if len(values) <= recorded:
                continue
            try:
                raw = float(values[recorded])
            except ValueError:
                continue
            unit = match.group(1).lower()
            to_pb = {"/ub": 1.0e-6, "/nb": 1.0e-3, "/pb": 1.0, "/fb": 1.0e3}
            if unit not in to_pb:
                raise ValueError(f"Unsupported BRIL luminosity unit {unit}")
            return raw * to_pb[unit], {"raw_recorded": raw, "raw_unit": unit,
                                       "source": str(filename.resolve())}
    raise ValueError(f"No totrecorded summary row found in {filename}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("processed_lumis", type=Path)
    parser.add_argument("--golden-json", type=Path, required=True)
    parser.add_argument("--brilcalc-csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    processed = expand(json.loads(args.processed_lumis.read_text()))
    certified = expand(json.loads(args.golden_json.read_text()))
    selected = processed & certified
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    selected_file = output / "processed_certified_lumis.json"
    selected_file.write_text(json.dumps(compress(selected), indent=2) + "\n")
    integrated = None
    bril = None
    if args.brilcalc_csv:
        integrated, bril = parse_brilcalc(args.brilcalc_csv)
    summary = {
        "schema_version": 1,
        "processed_sections": len(processed),
        "certified_sections": len(certified),
        "processed_and_certified_sections": len(selected),
        "processed_not_certified_sections": len(processed - certified),
        "certified_not_observed_sections": len(certified - processed),
        "processed_certified_lumis_json": str(selected_file),
        "integrated_luminosity_pb_inverse": integrated,
        "brilcalc": bril,
        "next_command": (
            f"brilcalc lumi -u /pb -i {selected_file} --normtag PATH_TO_APPROVED_NORMTAG -o luminosity.csv"
            if args.brilcalc_csv is None else None
        ),
        "warning": "The event selection must use the same golden JSON; intersection after processing cannot remove uncertified events already counted.",
    }
    (output / "luminosity_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
