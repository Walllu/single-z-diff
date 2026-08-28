#!/usr/bin/env python3
"""Apply global normalization and a configurable uncertainty model to merged inputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def quadrature(*values: float) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in values))


def ratio_uncertainty(value: float | None, denominator: dict[str, Any]) -> float | None:
    if value is None or not 0.0 <= value <= 1.0:
        return None
    sumw = float(denominator["sum_weights"])
    sumw2 = float(denominator["sum_weights_squared"])
    if sumw2 <= 0.0 or sumw == 0.0:
        return None
    effective_entries = sumw * sumw / sumw2
    return math.sqrt(max(0.0, value * (1.0 - value) / effective_entries))


def sample_factor(sample: dict[str, Any], metadata: dict[str, Any], luminosity: float) -> float:
    cross_section = metadata.get("cross_section_pb")
    if cross_section is None:
        raise ValueError(f"Missing cross_section_pb for {sample['process_name']}")
    sumw = metadata.get("sum_generator_weights")
    if sumw is None:
        sumw = sample["processed_sum_generator_weights"]
    if float(sumw) == 0.0:
        raise ValueError(f"Zero generator-weight denominator for {sample['process_name']}")
    return (luminosity * float(cross_section) * float(metadata.get("filter_efficiency", 1.0))
            * float(metadata.get("matching_efficiency", 1.0))
            * float(metadata.get("k_factor", 1.0)) / float(sumw))


def scaled_regions(sample: dict[str, Any], metadata: dict[str, Any], luminosity: float,
                   anti_mode: str) -> tuple[dict[str, float], dict[str, float], float]:
    factor = sample_factor(sample, metadata, luminosity)
    source = sample["regions_by_anti_isolation"][anti_mode]
    yields = {region: factor * float(value["sum_weights"]) for region, value in source.items()}
    variances = {region: factor * factor * float(value["sum_weights_squared"])
                 for region, value in source.items()}
    return yields, variances, factor


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("measurement_inputs", type=Path)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    merged = json.loads(args.measurement_inputs.read_text())
    config = json.loads(args.config.read_text())
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    data = merged.get("data")
    signal = merged.get("signal")
    if data is None or signal is None:
        raise RuntimeError("Finalization requires merged data and signal inputs")
    anti_mode = config["nonprompt_background"]["anti_isolation_mode"]
    data_regions = data["regions_by_anti_isolation"][anti_mode]
    observed_a = float(data_regions["A"])
    luminosity_value = config["luminosity"].get("integrated_pb_inverse")
    luminosity_summary = config["luminosity"].get("summary_json")
    if luminosity_value is None and luminosity_summary:
        summary_path = Path(luminosity_summary)
        if not summary_path.is_absolute():
            summary_path = (args.config.resolve().parent / summary_path).resolve()
        luminosity_value = json.loads(summary_path.read_text()).get(
            "integrated_luminosity_pb_inverse"
        )
    luminosity = None if luminosity_value is None else float(luminosity_value)

    acceptance = signal["acceptance"]
    efficiency = signal["efficiency"]
    acceptance_stat = ratio_uncertainty(acceptance, signal["yields"]["truth_total"])
    efficiency_stat = ratio_uncertainty(
        efficiency, signal["yields"]["truth_fiducial_pileup_weighted"]
    )
    efficiency_variations = {
        name: value["relative_to_nominal"]
        for name, value in signal.get("variation_efficiencies", {}).items()
    }
    sf_groups: dict[str, float] = {}
    for group in ("reco_sys", "reco_stat", "iso_sys", "iso_stat", "trigger_sys", "trigger_stat"):
        shifts = [abs(float(efficiency_variations[name])) for name in (f"{group}_up", f"{group}_down")
                  if efficiency_variations.get(name) is not None]
        sf_groups[group] = max(shifts, default=0.0)

    prompt_mode = config["prompt_backgrounds"]["mode"]
    prompt_a = 0.0
    prompt_a_mc_variance = 0.0
    prompt_norm_uncertainty = 0.0
    prompt_bcd = {region: 0.0 for region in "BCD"}
    prompt_bcd_variance = {region: 0.0 for region in "BCD"}
    prompt_details: dict[str, Any] = {}
    nonprompt_cfg = config["nonprompt_background"]
    if nonprompt_cfg.get("subtract_signal_leakage", False):
        if luminosity is None:
            raise ValueError("Signal leakage subtraction in B/C/D requires integrated luminosity")
        signal_metadata = config.get("signal_sample", {})
        signal_regions, signal_variances, signal_factor = scaled_regions(
            signal, signal_metadata, luminosity, anti_mode
        )
        for region in "BCD":
            prompt_bcd[region] += signal_regions[region]
            prompt_bcd_variance[region] += signal_variances[region]
        leakage_yields = dict(signal_regions)
        leakage_yields["A"] = 0.0
        leakage_variances = dict(signal_variances)
        leakage_variances["A"] = 0.0
        prompt_details["z_to_mumu_control_region_leakage"] = {
            "normalization_factor": signal_factor,
            "region_yields": leakage_yields,
            "region_mc_variances": leakage_variances,
            "metadata": signal_metadata,
            "note": "Only B/C/D are subtracted; A is the measured signal.",
        }
    if prompt_mode == "explicit":
        if luminosity is None:
            raise ValueError("Explicit prompt-background normalization requires integrated luminosity")
        metadata_by_name = {item["process_name"]: item
                            for item in config["prompt_backgrounds"].get("samples", [])}
        for name, sample in merged.get("backgrounds", {}).items():
            if name not in metadata_by_name:
                raise ValueError(f"No normalization metadata for merged background {name}")
            metadata = metadata_by_name[name]
            yields, variances, factor = scaled_regions(sample, metadata, luminosity, anti_mode)
            prompt_a += yields["A"]
            prompt_a_mc_variance += variances["A"]
            prompt_norm_uncertainty = quadrature(
                prompt_norm_uncertainty,
                yields["A"] * float(metadata.get("normalization_uncertainty_fraction", 0.0)),
            )
            for region in "BCD":
                prompt_bcd[region] += yields[region]
                prompt_bcd_variance[region] += variances[region]
            prompt_details[name] = {"normalization_factor": factor, "region_yields": yields,
                                    "region_mc_variances": variances, "metadata": metadata}
        missing_metadata = sorted(set(metadata_by_name) - set(merged.get("backgrounds", {})))
        if missing_metadata:
            raise ValueError(f"Configured prompt samples have no merged input: {missing_metadata}")
    elif prompt_mode != "missing_uncertainty":
        raise ValueError("prompt_backgrounds.mode must be explicit or missing_uncertainty")

    closure = float(nonprompt_cfg.get("closure_factor", 1.0))
    corrected = {region: float(data_regions[region]) - prompt_bcd[region] for region in "BCD"}
    if corrected["D"] <= 0.0:
        nonprompt = None
        nonprompt_stat_variance = None
    else:
        nonprompt = closure * corrected["B"] * corrected["C"] / corrected["D"]
        derivatives = {
            "B": closure * corrected["C"] / corrected["D"],
            "C": closure * corrected["B"] / corrected["D"],
            "D": -closure * corrected["B"] * corrected["C"] / corrected["D"] ** 2,
        }
        nonprompt_stat_variance = sum(
            derivatives[region] ** 2
            * (float(data_regions[region]) + prompt_bcd_variance[region])
            for region in "BCD"
        )
    nonprompt_uncertainty = None
    prompt_bcd_norm_uncertainty = 0.0
    prompt_combined_norm_uncertainty = prompt_norm_uncertainty
    if nonprompt is not None:
        if prompt_details:
            combined_shifts = []
            for detail in prompt_details.values():
                correlated_shift = sum(
                    derivatives[region] * detail["region_yields"][region]
                    for region in "BCD"
                )
                fraction = float(detail["metadata"].get(
                    "normalization_uncertainty_fraction", 0.0
                ))
                prompt_bcd_norm_uncertainty = quadrature(
                    prompt_bcd_norm_uncertainty, correlated_shift * fraction
                )
                combined_shifts.append(
                    (detail["region_yields"]["A"] + correlated_shift) * fraction
                )
            prompt_combined_norm_uncertainty = quadrature(*combined_shifts)
        nonprompt_uncertainty = quadrature(
            math.sqrt(nonprompt_stat_variance or 0.0),
            abs(nonprompt) * float(nonprompt_cfg.get("closure_factor_uncertainty_fraction", 0.0)),
            abs(nonprompt) * float(nonprompt_cfg.get("additional_uncertainty_fraction", 0.0)),
            prompt_bcd_norm_uncertainty,
        )

    candidate_before_prompt = None if nonprompt is None else observed_a - nonprompt
    signal_yield = None if candidate_before_prompt is None else candidate_before_prompt - prompt_a
    missing_prompt_uncertainty = 0.0
    if prompt_mode == "missing_uncertainty" and candidate_before_prompt is not None:
        missing = config["prompt_backgrounds"].get("missing_component", {})
        envelope_events = 0.0
        envelope_file = missing.get("envelope_json")
        if envelope_file:
            envelope_path = Path(envelope_file)
            if not envelope_path.is_absolute():
                envelope_path = (args.config.resolve().parent / envelope_path).resolve()
            envelope_payload = json.loads(envelope_path.read_text())
            if not envelope_payload.get("available", True):
                raise ValueError(
                    "Configured background envelope is unavailable: "
                    + str(envelope_payload.get("reason", "insufficient control-region statistics"))
                )
            envelope_events = float(envelope_payload["recommended_absolute_events"])
        missing_prompt_uncertainty = quadrature(
            float(missing.get("absolute_events", 0.0)),
            abs(candidate_before_prompt) * float(missing.get("fraction_of_candidate_yield", 0.0)),
            envelope_events,
        )

    ready = (luminosity is not None and luminosity > 0.0 and acceptance not in (None, 0.0)
             and efficiency not in (None, 0.0) and signal_yield is not None and signal_yield > 0.0)
    fiducial = full = None
    uncertainties: dict[str, Any] = {}
    if ready:
        fiducial = signal_yield / (float(efficiency) * luminosity)
        full = fiducial / float(acceptance)
        yield_stat = quadrature(math.sqrt(observed_a), math.sqrt(nonprompt_stat_variance or 0.0),
                                math.sqrt(prompt_a_mc_variance))
        yield_background = quadrature(
            abs(nonprompt or 0.0) * float(nonprompt_cfg.get("closure_factor_uncertainty_fraction", 0.0)),
            abs(nonprompt or 0.0) * float(nonprompt_cfg.get("additional_uncertainty_fraction", 0.0)),
            prompt_combined_norm_uncertainty,
            missing_prompt_uncertainty,
        )
        external = config.get("systematics", {})
        efficiency_relative = quadrature(
            *(sf_groups.values()),
            float(external.get("efficiency_additional_fraction", 0.0)),
            0.0 if efficiency_stat is None else efficiency_stat / float(efficiency),
        )
        acceptance_relative = quadrature(
            float(external.get("acceptance_modeling_fraction", 0.0)),
            0.0 if acceptance_stat is None else acceptance_stat / float(acceptance),
        )
        lumi_relative = float(config["luminosity"].get("uncertainty_fraction", 0.0))
        common_relative = quadrature(
            yield_stat / signal_yield, yield_background / signal_yield,
            efficiency_relative, lumi_relative,
            float(external.get("muon_momentum_fraction", 0.0)),
            float(external.get("pileup_fraction", 0.0)),
        )
        full_relative = quadrature(common_relative, acceptance_relative)
        uncertainties = {
            "yield_stat_events": yield_stat,
            "background_systematic_events": yield_background,
            "missing_prompt_events": missing_prompt_uncertainty,
            "efficiency_relative": efficiency_relative,
            "acceptance_relative": acceptance_relative,
            "luminosity_relative": lumi_relative,
            "fiducial_total_pb": abs(fiducial) * common_relative,
            "full_total_pb": abs(full) * full_relative,
            "muon_sf_relative_components": sf_groups,
        }

    signal_diagnostic = None
    signal_metadata = config.get("signal_sample", {})
    if luminosity is not None and signal_metadata.get("cross_section_pb") is not None:
        factor = sample_factor(signal, signal_metadata, luminosity)
        signal_diagnostic = {
            "normalization_factor": factor,
            "expected_selected_yield": factor * float(signal["yields"]["selected"]["sum_weights"]),
            "cross_section_pb": signal_metadata["cross_section_pb"],
            "filter_efficiency": signal_metadata.get("filter_efficiency", 1.0),
            "matching_efficiency": signal_metadata.get("matching_efficiency", 1.0),
            "k_factor": signal_metadata.get("k_factor", 1.0),
        }

    missing_cfg = config["prompt_backgrounds"].get("missing_component", {})
    has_missing_envelope = bool(
        missing_cfg.get("envelope_json")
        or float(missing_cfg.get("absolute_events", 0.0)) > 0.0
        or float(missing_cfg.get("fraction_of_candidate_yield", 0.0)) > 0.0
    )
    result = {
        "schema_version": 1,
        "title": "Inclusive Z to dimuon cross-section finalization",
        "ready_for_absolute_result": ready,
        "formulae": {
            "fiducial": "(N_A - N_nonprompt - N_prompt) / (epsilon * luminosity)",
            "full": "(N_A - N_nonprompt - N_prompt) / (acceptance * epsilon * luminosity)",
            "mc_normalization": "luminosity * cross_section * filter_efficiency * matching_efficiency * k_factor / sum_generator_weights",
        },
        "inputs": {"measurement_inputs": str(args.measurement_inputs.resolve()),
                   "finalization_config": str(args.config.resolve()),
                   "configuration": config},
        "luminosity_pb_inverse": luminosity,
        "acceptance": acceptance,
        "acceptance_statistical_estimate": acceptance_stat,
        "efficiency": efficiency,
        "efficiency_statistical_estimate": efficiency_stat,
        "observed_A": observed_a,
        "nonprompt_background": {"value": nonprompt, "total_uncertainty": nonprompt_uncertainty,
                                 "prompt_subtracted_BCD": corrected,
                                 "anti_isolation_mode": anti_mode, "closure_factor": closure},
        "prompt_background": {"mode": prompt_mode, "value": prompt_a,
                              "mc_statistical_variance": prompt_a_mc_variance,
                              "normalization_uncertainty": prompt_norm_uncertainty,
                              "normalization_effect_on_abcd_uncertainty": prompt_bcd_norm_uncertainty,
                              "combined_normalization_effect_on_signal_yield": prompt_combined_norm_uncertainty,
                              "missing_component_uncertainty": missing_prompt_uncertainty,
                              "ignored_merged_backgrounds": (
                                  sorted(merged.get("backgrounds", {}))
                                  if prompt_mode == "missing_uncertainty" else []
                              ),
                              "processes": prompt_details},
        "signal_yield": signal_yield,
        "fiducial_cross_section_pb": fiducial,
        "full_cross_section_pb": full,
        "uncertainties": uncertainties,
        "signal_mc_normalization_diagnostic": signal_diagnostic,
        "external_inputs_still_required": [
            name for name, condition in (
                ("approved integrated luminosity and uncertainty", luminosity is None),
                ("prompt-background samples and normalization metadata, or a justified missing-component uncertainty",
                 prompt_mode == "missing_uncertainty" and not has_missing_envelope),
                ("signal cross section/filter efficiency/k-factor for the MC normalization diagnostic",
                 signal_metadata.get("cross_section_pb") is None),
            ) if condition
        ],
    }
    (output / "cross_section_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in
                      ("ready_for_absolute_result", "signal_yield", "fiducial_cross_section_pb",
                       "full_cross_section_pb", "external_inputs_still_required")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
