#!/usr/bin/env python3
"""Compare a pasted CERNBox listing with pasted ``ls -lh`` output."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any


UUID = r"[0-9A-F]{8}-[0-9A-F]{4}-[0-9A-F]{4}-[0-9A-F]{4}-[0-9A-F]{12}"


def decimal_size(value: str, unit: str) -> tuple[float, float]:
    multipliers = {"KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}
    multiplier = multipliers[unit]
    decimals = len(value.partition(".")[2]) if "." in value else 0
    estimate = float(value) * multiplier
    rounding = 0.5 * 10.0 ** (-decimals) * multiplier
    return estimate, rounding


def binary_size(value: str, unit: str) -> tuple[float, float]:
    multipliers = {"K": 1024.0, "M": 1024.0**2, "G": 1024.0**3, "T": 1024.0**4}
    multiplier = multipliers[unit]
    decimals = len(value.partition(".")[2]) if "." in value else 0
    estimate = float(value) * multiplier
    # GNU ls human-readable sizes are aggressively rounded; allow one displayed
    # unit at integer precision plus a small fixed margin in the comparison.
    rounding = 10.0 ** (-decimals) * multiplier
    return estimate, rounding


def parse_cernbox(filename: Path) -> dict[str, dict[str, Any]]:
    text = filename.read_text()
    pattern = re.compile(
        rf"(?m)^({UUID})\s*\n\.root\s*\n\s*\n?\s*([0-9]+(?:\.[0-9]+)?)\s*(KB|MB|GB|TB)\b"
    )
    result = {}
    for name, value, unit in pattern.findall(text):
        estimate, rounding = decimal_size(value, unit)
        result[f"{name}.root"] = {
            "display": f"{value} {unit}", "bytes_estimate": estimate,
            "rounding_bytes": rounding,
        }
    return result


def parse_ls(filename: Path) -> tuple[dict[str, dict[str, Any]], str | None]:
    result = {}
    total = None
    line_pattern = re.compile(
        rf"\s([0-9]+(?:\.[0-9]+)?)([KMGT])\s+\S+\s+\d+\s+\d{{2}}:\d{{2}}\s+({UUID}\.root)$"
    )
    for line in filename.read_text().splitlines():
        if line.startswith("total "):
            total = line.split(maxsplit=1)[1]
        match = line_pattern.search(line)
        if not match:
            continue
        value, unit, name = match.groups()
        estimate, rounding = binary_size(value, unit)
        result[name] = {
            "display": f"{value}{unit}", "bytes_estimate": estimate,
            "rounding_bytes": rounding,
        }
    return result, total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cernbox", type=Path)
    parser.add_argument("cephfs", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cernbox = parse_cernbox(args.cernbox)
    cephfs, ls_total = parse_ls(args.cephfs)
    if not cernbox or not cephfs:
        raise RuntimeError(
            f"Failed to parse listings: CERNBox={len(cernbox)}, CephFS={len(cephfs)} files"
        )
    shared = sorted(cernbox.keys() & cephfs.keys())
    missing = sorted(cernbox.keys() - cephfs.keys())
    extra = sorted(cephfs.keys() - cernbox.keys())
    consistent = []
    inconsistent = []
    for name in shared:
        left, right = cernbox[name], cephfs[name]
        difference = abs(left["bytes_estimate"] - right["bytes_estimate"])
        tolerance = left["rounding_bytes"] + right["rounding_bytes"] + 2.0 * 1024.0**2
        row = {
            "basename": name,
            "cernbox_display": left["display"],
            "cephfs_display": right["display"],
            "estimated_difference_bytes": difference,
            "rounding_tolerance_bytes": tolerance,
        }
        (consistent if difference <= tolerance else inconsistent).append(row)

    cernbox_sum = sum(item["bytes_estimate"] for item in cernbox.values())
    cephfs_sum = sum(item["bytes_estimate"] for item in cephfs.values())
    missing_sum = sum(cernbox[name]["bytes_estimate"] for name in missing)
    extra_sum = sum(cephfs[name]["bytes_estimate"] for name in extra)
    report = {
        "title": "Comparison of pasted CERNBox and CephFS human-readable listings",
        "inputs": {"cernbox": str(args.cernbox), "cephfs": str(args.cephfs)},
        "counts": {
            "cernbox": len(cernbox), "cephfs": len(cephfs), "shared": len(shared),
            "shared_size_consistent": len(consistent),
            "shared_size_inconsistent": len(inconsistent),
            "missing_from_cephfs": len(missing), "extra_on_cephfs": len(extra),
        },
        "display_size_estimates_bytes": {
            "cernbox_sum": cernbox_sum,
            "cephfs_sum": cephfs_sum,
            "cernbox_missing_from_cephfs_sum": missing_sum,
            "cephfs_extra_sum": extra_sum,
            "cephfs_ls_reported_total": ls_total,
            "note": "CERNBox units are decimal; ls -h units are binary and rounded.",
        },
        "missing_from_cephfs": missing,
        "extra_on_cephfs": extra,
        "shared_size_inconsistent": inconsistent,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    missing_file = args.output.with_name("missing_from_cephfs.txt")
    missing_file.write_text("".join(f"{name}\n" for name in missing))
    replace_file = args.output.with_name("replace_on_cephfs.txt")
    replace_file.write_text(
        "".join(f"{item['basename']}\n" for item in inconsistent)
    )
    print(json.dumps({"counts": report["counts"], "sizes": report["display_size_estimates_bytes"]}, indent=2))
    print(f"Wrote {args.output}, {missing_file}, and {replace_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
