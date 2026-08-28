#!/usr/bin/env python3
"""Derive a conservative unresolved-background envelope from ABCD closure output."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def safe_ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def mode_components(payload: dict[str, Any], mode: str, confidence_z: float) -> dict[str, Any]:
    result = payload["results"][mode]
    closure = result.get("closure_factor") or {}
    if closure.get("value") is None or closure.get("variance") is None:
        return {
            "available": False,
            "reason": "ABCD closure factor is undefined, usually because a required sideband control region is empty",
            "closure_factor": closure.get("value"),
            "closure_factor_stat_uncertainty": None,
        }
    windows = result["windows"]
    signal = windows["signal"]
    kappa = float(result["closure_factor"]["value"])
    kappa_variance = float(result["closure_factor"]["variance"])
    base_signal = float(signal["base_bc_over_d"])
    nominal_background = kappa * base_signal
    data_a = float(signal["components"]["A"]["data"])
    sideband_variations: dict[str, Any] = {}
    residual_bounds: dict[str, Any] = {}
    for name in ("low_sideband", "high_sideband"):
        window = windows[name]
        base = float(window["base_bc_over_d"])
        observed = float(window["observed_A_nonprompt"])
        sideband_kappa = safe_ratio(observed, base)
        varied_background = None if sideband_kappa is None else sideband_kappa * base_signal
        sideband_variations[name] = {
            "closure_factor": sideband_kappa,
            "signal_window_background": varied_background,
            "absolute_signal_yield_shift": (
                None if varied_background is None else abs(varied_background - nominal_background)
            ),
        }

        prediction = kappa * base
        residual = observed - prediction
        observed_variance = float(window["observed_A_variance"])
        base_variance = float(window["base_variance"])
        # The closure factor and sideband counts are not fully independent. This
        # deliberately neglects their negative covariance, making the bound conservative.
        residual_variance = max(
            0.0, observed_variance + kappa * kappa * base_variance
            + base * base * kappa_variance
        )
        width = float(window["bounds_gev"][1]) - float(window["bounds_gev"][0])
        signal_width = (float(signal["bounds_gev"][1])
                        - float(signal["bounds_gev"][0]))
        upper = max(0.0, residual + confidence_z * math.sqrt(residual_variance))
        residual_bounds[name] = {
            "residual_events": residual,
            "residual_stat_uncertainty": math.sqrt(residual_variance),
            "one_sided_upper_events_in_sideband": upper,
            "width_scaled_upper_events_in_signal_window": upper * signal_width / width,
        }
    return {
        "available": True,
        "closure_factor": kappa,
        "closure_factor_stat_uncertainty": math.sqrt(max(0.0, kappa_variance)),
        "nominal_signal_window_background": nominal_background,
        "candidate_signal_yield": data_a - nominal_background,
        "sideband_closure_variations": sideband_variations,
        "sideband_residual_bounds": residual_bounds,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("closure_summary", type=Path)
    parser.add_argument("--nominal-mode", choices=("both", "at_least_one"), default="both")
    parser.add_argument("--confidence-z", type=float, default=1.64,
                        help="One-sided Gaussian quantile; 1.64 is approximately 95%%")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.closure_summary.read_text())
    modes = {mode: mode_components(source, mode, args.confidence_z)
             for mode in ("both", "at_least_one")}
    nominal = modes[args.nominal_mode]
    if not nominal["available"]:
        payload = {
            "schema_version": 1,
            "title": "Conservative unresolved-background envelope from sideband closure",
            "status": "insufficient_statistics",
            "available": False,
            "source_closure_summary": str(args.closure_summary.resolve()),
            "nominal_mode": args.nominal_mode,
            "confidence_z": args.confidence_z,
            "modes": modes,
            "components_events": {},
            "recommended_absolute_events": None,
            "recommended_fraction_of_candidate_yield": None,
            "reason": nominal["reason"],
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
        print(json.dumps({"output": str(args.output.resolve()),
                          "status": payload["status"],
                          "reason": payload["reason"]}, indent=2))
        return 0
    alternate_mode = "at_least_one" if args.nominal_mode == "both" else "both"
    components: dict[str, float] = {}
    if modes[alternate_mode]["available"]:
        components["anti_isolation_definition_shift"] = abs(
            modes[alternate_mode]["nominal_signal_window_background"]
            - nominal["nominal_signal_window_background"]
        )
    for sideband, value in nominal["sideband_closure_variations"].items():
        shift = value["absolute_signal_yield_shift"]
        if shift is not None:
            components[f"{sideband}_closure_shift"] = float(shift)
    for sideband, value in nominal["sideband_residual_bounds"].items():
        components[f"{sideband}_residual_95_upper"] = float(
            value["width_scaled_upper_events_in_signal_window"]
        )
    recommended = max(components.values(), default=0.0)
    candidate = float(nominal["candidate_signal_yield"])
    payload = {
        "schema_version": 1,
        "title": "Conservative unresolved-background envelope from sideband closure",
        "status": "available",
        "available": True,
        "source_closure_summary": str(args.closure_summary.resolve()),
        "nominal_mode": args.nominal_mode,
        "confidence_z": args.confidence_z,
        "modes": modes,
        "components_events": components,
        "recommended_absolute_events": recommended,
        "recommended_fraction_of_candidate_yield": (
            recommended / abs(candidate) if candidate else None
        ),
        "interpretation": (
            "Use zero central yield and this symmetric event envelope while prompt samples "
            "are absent. It covers the largest anti-isolation, sideband-transfer, or "
            "one-sided residual excursion and is therefore a total unresolved-background "
            "model uncertainty, not a process-specific cross-section uncertainty."
        ),
        "double_counting_warning": (
            "Do not add the same anti-isolation/low-high closure shifts again as independent "
            "systematics. Statistical propagation of the disjoint A/B/C/D counts remains separate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"output": str(args.output.resolve()),
                      "recommended_absolute_events": recommended,
                      "recommended_fraction_of_candidate_yield": payload["recommended_fraction_of_candidate_yield"],
                      "dominant_component": max(components, key=components.get) if components else None}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
