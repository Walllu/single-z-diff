#!/usr/bin/env python3
"""Apply fitted Z->mumu pair-region scale and resolution corrections to MC muons."""

from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Mapping


ETA_REGIONS = {
    "neg_endcap": (-2.5, -2.1), "neg_transition": (-2.1, -1.4),
    "barrel": (-1.4, 1.4), "pos_transition": (1.4, 2.1),
    "pos_endcap": (2.1, 2.5),
}
DEFAULT_PARAMETERS_FILE = (
    Path(__file__).resolve().parent / "plots" /
    "scan_poisson_seed_nominal_s21_m0p005_p0p005_r21_p0p000_p0p050_bins_tiered_v1" /
    "summary.json"
)


@dataclass(frozen=True)
class Muon:
    pt: float
    eta: float
    phi: float = 0.0
    mass: float = 0.105658


def _muon(value: Muon | Mapping[str, float]) -> Muon:
    return value if isinstance(value, Muon) else Muon(**value)


def _eta_region(eta: float) -> str | None:
    names = list(ETA_REGIONS)
    for position, name in enumerate(names):
        low, high = ETA_REGIONS[name]
        if low <= eta < high or (position == len(names) - 1 and eta == high):
            return name
    return None


@lru_cache(maxsize=8)
def _load_parameters(filename: str) -> dict:
    return json.loads(Path(filename).read_text())


def correct_muon_pair(
    muon1: Muon | Mapping[str, float],
    muon2: Muon | Mapping[str, float],
    parameters_file: str | Path = DEFAULT_PARAMETERS_FILE,
    *,
    apply_correction: bool = True,
    rng_seed: int = 314159,
    event_key: str | int | None = None,
) -> dict:
    """Return corrected MC muons and factors using the nominal leading/subleading pair region.

    Pass a stable event_key (for example ``f"{run}:{lumi}:{event}"``) so each
    event receives different but reproducible Gaussian draws. Missing/excluded
    eta-pair regions return identity factors. These fitted corrections are for
    simulation; do not smear collision data.
    """
    original = [_muon(muon1), _muon(muon2)]
    order = sorted(range(2), key=lambda i: original[i].pt, reverse=True)
    lead_region, sublead_region = _eta_region(original[order[0]].eta), _eta_region(original[order[1]].eta)
    region = f"{lead_region}__{sublead_region}" if lead_region and sublead_region else None
    payload = _load_parameters(str(Path(parameters_file).resolve()))
    fit = payload.get("regions", {}).get(region, {}).get("best") if region else None
    scale, resolution = (float(fit["scale"]), float(fit["resolution"])) if fit else (0.0, 0.0)
    rng = random.Random(f"{rng_seed}:{event_key if event_key is not None else 'default'}")
    gaussians = [rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)]
    proposed_factors = [1.0 + scale + resolution * gaussians[i] for i in range(2)]
    applied_factors = proposed_factors if apply_correction and fit else [1.0, 1.0]
    corrected = []
    for muon, factor in zip(original, applied_factors):
        values = asdict(muon); values["pt"] *= factor; corrected.append(values)
    return {"muons": corrected, "factors": applied_factors,
            "proposed_factors": proposed_factors, "gaussians": gaussians,
            "region": region, "scale": scale, "resolution": resolution,
            "applied": bool(apply_correction and fit), "covered": bool(fit)}


__all__ = ["DEFAULT_PARAMETERS_FILE", "Muon", "correct_muon_pair"]
