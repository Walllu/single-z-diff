#!/usr/bin/env python3
"""Render SVG summaries of correction impacts and efficiency-SF variations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np


HERE = Path(__file__).resolve().parent
RUNS = {
    "Uncorrected": "local_uncorrected_measurement",
    "Momentum only": "local_momentum_only_measurement",
    "Efficiency SF only": "local_sf_only_measurement",
    "Nominal corrected": "local_corrected_measurement",
    "Nominal + toy PU": "local_corrected_measurement_toy_pileup",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=HERE / "outputs")
    parser.add_argument("--output-dir", type=Path, default=HERE / "plots")
    parser.add_argument("--label", default="local_measurement_validation")
    return parser.parse_args()


def main() -> int:
    args = arguments()
    payloads = {}
    for label, directory in RUNS.items():
        filename = args.input_dir / directory / "measurement_summary.json"
        if not filename.is_file():
            raise FileNotFoundError(f"Missing measurement result {filename}")
        payloads[label] = json.loads(filename.read_text())["results"]["mc"]
    output = args.output_dir / args.label
    output.mkdir(parents=True, exist_ok=True)
    hep.style.use("CMS")

    labels = list(payloads)
    efficiency = np.asarray([payloads[name]["efficiency"] for name in labels])
    correction = np.asarray([payloads[name]["acceptance_times_efficiency"] for name in labels])
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.bar(x - 0.2, efficiency, 0.4, color="#3b82f6", label=r"$\epsilon$")
    ax.bar(x + 0.2, correction, 0.4, color="#f59e0b", label=r"$A\epsilon$")
    for xpos, value in zip(x - 0.2, efficiency):
        ax.text(xpos, value + 0.006, f"{value:.3f}", ha="center", fontsize=9)
    for xpos, value in zip(x + 0.2, correction):
        ax.text(xpos, value + 0.006, f"{value:.3f}", ha="center", fontsize=9)
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylim(0.0, 0.95)
    ax.set_ylabel("Correction factor")
    ax.legend()
    hep.cms.label("Simulation", data=False, com=13, ax=ax)
    fig.savefig(output / "acceptance_efficiency_correction_impact.svg", bbox_inches="tight")
    plt.close(fig)

    nominal = payloads["Nominal corrected"]
    variations = nominal["variation_efficiencies"]
    names = list(variations)
    shifts = 100.0 * np.asarray([variations[name]["relative_to_nominal"] for name in names])
    colors = ["#ef4444" if value >= 0 else "#3b82f6" for value in shifts]
    fig, ax = plt.subplots(figsize=(10, 7))
    ypos = np.arange(len(names))
    ax.barh(ypos, shifts, color=colors)
    ax.axvline(0.0, color="black", linewidth=1.0)
    ax.set_yticks(ypos, [name.replace("_", " ") for name in names])
    ax.set_xlabel("Relative efficiency shift [%]")
    ax.grid(axis="x", alpha=0.2)
    hep.cms.label("Simulation", data=False, com=13, ax=ax)
    fig.savefig(output / "muon_efficiency_sf_variations.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote 2 SVG plots to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
