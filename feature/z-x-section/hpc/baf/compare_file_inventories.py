#!/usr/bin/env python3
"""Inventory ROOT files and compare a complete source with a partial copy."""

from __future__ import annotations

import argparse
import json
import socket
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    scan = subparsers.add_parser("scan", help="Inventory ROOT basenames and exact byte sizes")
    scan.add_argument("directory", type=Path)
    scan.add_argument("--output", type=Path, required=True)
    scan.add_argument("--max-files", type=int, default=-1,
                      help="Diagnostic limit after sorting; omit for a production inventory")

    compare = subparsers.add_parser(
        "compare", help="Compare a known-complete reference inventory with a candidate copy"
    )
    compare.add_argument("reference", type=Path,
                         help="Inventory of the known-complete/source directory")
    compare.add_argument("candidate", type=Path,
                         help="Inventory of the possibly partial destination directory")
    compare.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def human_bytes(value: int) -> str:
    amount = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if amount < 1024.0 or unit == "TiB":
            return f"{amount:.2f} {unit}"
        amount /= 1024.0
    raise AssertionError("unreachable")


def scan(directory: Path, output: Path, maximum: int) -> None:
    root = directory.expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Directory does not exist: {root}")
    paths = sorted(path for path in root.rglob("*.root") if path.is_file())
    if maximum >= 0:
        paths = paths[:maximum]
    if not paths:
        raise RuntimeError(f"No ROOT files found below {root}")
    files = [
        {
            "basename": path.name,
            "relative_path": path.relative_to(root).as_posix(),
            "size_bytes": path.stat().st_size,
        }
        for path in paths
    ]
    payload = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "root_directory": str(root),
        "selection": {"pattern": "*.root", "recursive": True, "max_files": maximum},
        "file_count": len(files),
        "total_size_bytes": sum(item["size_bytes"] for item in files),
        "files": files,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"Inventoried {payload['file_count']} files, "
        f"{human_bytes(payload['total_size_bytes'])}, below {root}"
    )
    print(f"Wrote {output}")


def load(filename: Path) -> dict[str, Any]:
    payload = json.loads(filename.read_text())
    if payload.get("schema_version") != 1 or not isinstance(payload.get("files"), list):
        raise ValueError(f"Unsupported or malformed inventory: {filename}")
    return payload


def unique_by_basename(payload: dict[str, Any], label: str) -> dict[str, dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in payload["files"]:
        grouped[item["basename"]].append(item)
    duplicates = {name: items for name, items in grouped.items() if len(items) > 1}
    if duplicates:
        examples = ", ".join(sorted(duplicates)[:10])
        raise RuntimeError(
            f"{label} inventory has {len(duplicates)} duplicate ROOT basenames "
            f"({examples}); basename matching would be ambiguous"
        )
    return {name: items[0] for name, items in grouped.items()}


def write_lines(filename: Path, values: list[str]) -> None:
    filename.write_text("".join(f"{value}\n" for value in values))


def compare(reference_file: Path, candidate_file: Path, output: Path) -> int:
    reference_payload = load(reference_file)
    candidate_payload = load(candidate_file)
    reference = unique_by_basename(reference_payload, "reference")
    candidate = unique_by_basename(candidate_payload, "candidate")

    exact: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    mismatched: list[dict[str, Any]] = []
    for name, source in sorted(reference.items()):
        destination = candidate.get(name)
        if destination is None:
            missing.append(source)
        elif int(source["size_bytes"]) == int(destination["size_bytes"]):
            exact.append(source)
        else:
            mismatched.append({
                "basename": name,
                "reference_relative_path": source["relative_path"],
                "reference_size_bytes": int(source["size_bytes"]),
                "candidate_relative_path": destination["relative_path"],
                "candidate_size_bytes": int(destination["size_bytes"]),
            })
    extra = [candidate[name] for name in sorted(candidate.keys() - reference.keys())]
    bytes_to_copy = (
        sum(int(item["size_bytes"]) for item in missing)
        + sum(int(item["reference_size_bytes"]) for item in mismatched)
    )
    summary = {
        "title": "ROOT file inventory comparison by basename and exact byte size",
        "reference_inventory": str(reference_file.resolve()),
        "candidate_inventory": str(candidate_file.resolve()),
        "reference_root": reference_payload.get("root_directory"),
        "candidate_root": candidate_payload.get("root_directory"),
        "counts": {
            "reference": len(reference),
            "candidate": len(candidate),
            "exact_matches": len(exact),
            "missing_from_candidate": len(missing),
            "size_mismatches": len(mismatched),
            "extra_in_candidate": len(extra),
        },
        "bytes_to_copy_or_replace": bytes_to_copy,
        "missing": missing,
        "size_mismatches": mismatched,
        "extra": extra,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_lines(output / "missing_from_candidate.txt",
                [item["relative_path"] for item in missing])
    write_lines(output / "replace_size_mismatches.txt",
                [item["reference_relative_path"] for item in mismatched])
    write_lines(
        output / "copy_or_replace_from_reference.txt",
        [item["relative_path"] for item in missing]
        + [item["reference_relative_path"] for item in mismatched],
    )
    write_lines(output / "extra_in_candidate.txt",
                [item["relative_path"] for item in extra])
    print(f"Reference: {len(reference)} files; candidate: {len(candidate)} files")
    print(
        f"Exact: {len(exact)}; missing: {len(missing)}; "
        f"wrong size: {len(mismatched)}; extra: {len(extra)}"
    )
    print(f"Copy/replace volume: {human_bytes(bytes_to_copy)}")
    print(f"Wrote comparison products to {output}")
    return 0 if not missing and not mismatched else 1


def main() -> int:
    args = arguments()
    if args.command == "scan":
        scan(args.directory, args.output, args.max_files)
        return 0
    return compare(args.reference, args.candidate, args.output_dir)


if __name__ == "__main__":
    raise SystemExit(main())
