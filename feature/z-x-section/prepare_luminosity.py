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
    candidates: list[tuple[float, str, bool]] = []
    summary_seen = False
    for index, line in enumerate(lines):
        if "summary" in line.lower():
            summary_seen = True
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
            if any(value.lower().startswith("totrecorded(") for value in values):
                break
            if len(values) <= recorded:
                continue
            try:
                raw = float(values[recorded])
            except ValueError:
                continue
            unit = match.group(1).lower()
            candidates.append((raw, unit, summary_seen))
    if candidates:
        summary_candidates = [candidate for candidate in candidates if candidate[2]]
        raw, unit, _ = (summary_candidates or candidates)[-1]
        to_pb = {"/ub": 1.0e-6, "/nb": 1.0e-3, "/pb": 1.0, "/fb": 1.0e3}
        if unit not in to_pb:
            raise ValueError(f"Unsupported BRIL luminosity unit {unit}")
        return raw * to_pb[unit], {"raw_recorded": raw, "raw_unit": unit,
                                   "source": str(filename.resolve())}
    raise ValueError(f"No totrecorded summary row found in {filename}")


def parse_lumibyls(
    filename: Path, selected: set[tuple[int, int]]
) -> tuple[float, dict[str, Any]]:
    """Sum an official BRIL by-lumisection CSV over an exact run/LS selection."""
    rows: dict[tuple[int, int], float] = {}
    unit: str | None = None
    header: list[str] | None = None
    positions: tuple[int, int, int] | None = None
    for raw_line in filename.read_text().splitlines():
        line = raw_line.strip().lstrip("#").strip()
        if not line:
            continue
        values = [value.strip() for value in next(csv.reader([line]))]
        run_position = next(
            (index for index, value in enumerate(values)
             if value.lower().startswith("run:fill")),
            None,
        )
        ls_position = next(
            (index for index, value in enumerate(values) if value.lower() == "ls"),
            None,
        )
        recorded_position = next(
            (index for index, value in enumerate(values)
             if value.lower().startswith("recorded(")),
            None,
        )
        if None not in (run_position, ls_position, recorded_position):
            header = values
            positions = (int(run_position), int(ls_position), int(recorded_position))
            match = re.search(r"\((/[^)]+)\)", values[int(recorded_position)])
            if match is None:
                raise ValueError(
                    f"Cannot determine luminosity unit from {values[int(recorded_position)]!r}"
                )
            unit = match.group(1).lower()
            continue
        if positions is None:
            continue
        run_position, ls_position, recorded_position = positions
        if len(values) <= max(positions):
            continue
        try:
            run = int(values[run_position].split(":", 1)[0])
            lumi = int(values[ls_position].split(":", 1)[0])
            recorded = float(values[recorded_position])
        except ValueError:
            continue
        key = (run, lumi)
        if key in rows and rows[key] != recorded:
            raise ValueError(f"Conflicting luminosity values for run/LS {key}")
        rows[key] = recorded
    if header is None or unit is None or not rows:
        raise ValueError(f"No run/LS luminosity rows found in {filename}")
    to_pb = {"/ub": 1.0e-6, "/nb": 1.0e-3, "/pb": 1.0, "/fb": 1.0e3}
    if unit not in to_pb:
        raise ValueError(f"Unsupported luminosity unit {unit}")
    matched = selected & set(rows)
    missing = selected - set(rows)
    integrated = sum(rows[key] for key in matched) * to_pb[unit]
    return integrated, {
        "kind": "official_luminosity_by_lumisection_csv",
        "source": str(filename.resolve()),
        "raw_unit": unit,
        "available_rows": len(rows),
        "matched_selected_sections": len(matched),
        "selected_sections_without_luminosity": len(missing),
        "missing_examples": [[run, lumi] for run, lumi in sorted(missing)[:10]],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("processed_lumis", type=Path)
    parser.add_argument("--golden-json", type=Path, required=True)
    luminosity_source = parser.add_mutually_exclusive_group()
    luminosity_source.add_argument("--brilcalc-csv", type=Path)
    luminosity_source.add_argument(
        "--lumibyls-csv", type=Path,
        help="Official luminosity-by-lumisection CSV, such as pp_2016lumibyls.csv",
    )
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
    source_details = None
    if args.brilcalc_csv:
        integrated, source_details = parse_brilcalc(args.brilcalc_csv)
    elif args.lumibyls_csv:
        integrated, source_details = parse_lumibyls(args.lumibyls_csv, selected)
    summary = {
        "schema_version": 2,
        "processed_sections": len(processed),
        "certified_sections": len(certified),
        "processed_and_certified_sections": len(selected),
        "processed_not_certified_sections": len(processed - certified),
        "certified_not_observed_sections": len(certified - processed),
        "processed_certified_lumis_json": str(selected_file),
        "integrated_luminosity_pb_inverse": integrated,
        "luminosity_source": source_details,
        "brilcalc": source_details if args.brilcalc_csv else None,
        "next_command": (
            f"brilcalc lumi -c web -u /pb -i {selected_file} "
            "--normtag PATH_TO_APPROVED_NORMTAG --output-style csv -o luminosity.csv"
            if args.brilcalc_csv is None and args.lumibyls_csv is None else None
        ),
        "warning": "The event selection must use the same golden JSON; intersection after processing cannot remove uncertified events already counted.",
    }
    (output / "luminosity_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
