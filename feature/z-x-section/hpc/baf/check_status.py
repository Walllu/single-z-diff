#!/usr/bin/env python3
"""Compare a prepared BAF submission's expected inputs with successful parts."""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from pathlib import Path


def manifest_paths(directory: Path, sample: str) -> set[str]:
    result: set[str] = set()
    for filename in sorted((directory / "manifests").glob(f"{sample}_*.txt")):
        result.update(line.strip() for line in filename.read_text().splitlines()
                      if line.strip() and not line.lstrip().startswith("#"))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("submit_dir", type=Path)
    args = parser.parse_args()
    submit = args.submit_dir.expanduser().resolve()
    metadata = json.loads((submit / "submission.json").read_text())
    buddy = os.environ.get("BUDDY")
    if not buddy:
        raise RuntimeError("BUDDY is not set")
    parts = Path(buddy) / "z-xsec/results" / metadata["run_label"] / "parts"
    complete = True
    for sample in ("data", "mc"):
        expected = manifest_paths(submit, sample)
        observed: list[str] = []
        successful_parts = 0
        for summary_file in sorted((parts / sample).glob("part_*/summary.json")):
            if not (summary_file.parent / "SUCCESS").is_file():
                continue
            successful_parts += 1
            payload = json.loads(summary_file.read_text())
            observed.extend(payload["samples"][sample]["files"])
        counts = Counter(observed)
        observed_set = set(observed)
        missing = sorted(expected - observed_set)
        unexpected = sorted(observed_set - expected)
        duplicates = sorted(path for path, count in counts.items() if count > 1)
        expected_jobs = metadata["samples"][sample]["jobs"]
        okay = not missing and not unexpected and not duplicates and successful_parts == expected_jobs
        complete &= okay
        print(
            f"{sample}: successful parts {successful_parts}/{expected_jobs}; "
            f"covered files {len(observed_set)}/{len(expected)}; "
            f"missing={len(missing)}, unexpected={len(unexpected)}, duplicates={len(duplicates)}"
        )
        for label, paths in (("missing", missing), ("unexpected", unexpected),
                             ("duplicate", duplicates)):
            for path in paths[:10]:
                print(f"  {label}: {path}")
            if len(paths) > 10:
                print(f"  ... and {len(paths) - 10} more {label} paths")
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
