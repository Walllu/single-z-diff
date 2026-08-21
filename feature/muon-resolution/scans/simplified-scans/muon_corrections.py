#!/usr/bin/env python3
"""Portable per-muon barrel/endcap scale and resolution corrections.

This module deliberately has no ROOT dependency.  It can therefore be copied or
imported by a downstream Python analysis.  Random draws are keyed by the event
and muon identity, making a chosen seed reproducible and independent of event
processing order.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Calibration:
    barrel_scale: float = 0.0
    barrel_resolution: float = 0.0
    endcap_scale: float = 0.0
    endcap_resolution: float = 0.0
    barrel_max_abs_eta: float = 1.4
    muon_max_abs_eta: float = 2.5

    @classmethod
    def from_summary(cls, filename: str | Path) -> "Calibration":
        payload = json.loads(Path(filename).read_text())
        values = payload.get("best_fit", payload)
        return cls(
            barrel_scale=float(values["barrel_scale"]),
            barrel_resolution=float(values["barrel_resolution"]),
            endcap_scale=float(values["endcap_scale"]),
            endcap_resolution=float(values["endcap_resolution"]),
        )

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class CorrectedMuon:
    pt: float
    eta: float
    phi: float
    mass: float
    factor: float
    region: str
    gaussian: float


def eta_region(eta: float, calibration: Calibration) -> str:
    absolute = abs(float(eta))
    if absolute < calibration.barrel_max_abs_eta:
        return "barrel"
    if absolute < calibration.muon_max_abs_eta:
        return "endcap"
    raise ValueError(f"Muon |eta|={absolute:g} lies outside the calibrated range")


def keyed_gaussian(seed: int, event_key: Any, muon_key: Any) -> float:
    """Return a stable N(0,1) draw without maintaining a sequential RNG stream."""
    if event_key is None:
        raise ValueError("A unique event_key is required for reproducible event-wise smearing")
    token = f"{int(seed)}:{event_key}:{muon_key}".encode("utf-8")
    digest = hashlib.sha256(token).hexdigest()
    return random.Random(digest).gauss(0.0, 1.0)


def correction_factor(
    eta: float,
    calibration: Calibration,
    *,
    seed: int = 314159,
    event_key: Any = None,
    muon_key: Any = 0,
    apply: bool = True,
) -> tuple[float, float, str]:
    """Return ``(factor, gaussian, region)`` for one simulated muon."""
    region = eta_region(eta, calibration)
    if not apply:
        return 1.0, 0.0, region
    gaussian = keyed_gaussian(seed, event_key, muon_key)
    if region == "barrel":
        factor = 1.0 + calibration.barrel_scale + calibration.barrel_resolution * gaussian
    else:
        factor = 1.0 + calibration.endcap_scale + calibration.endcap_resolution * gaussian
    if factor <= 0:
        raise ValueError(f"Non-positive momentum correction factor {factor:g}")
    return factor, gaussian, region


def _component(muon: Any, name: str, index: int) -> float:
    if isinstance(muon, dict):
        return float(muon[name])
    if hasattr(muon, name):
        return float(getattr(muon, name))
    return float(muon[index])


def correct_muon(
    muon: Any,
    calibration: Calibration,
    *,
    seed: int = 314159,
    event_key: Any = None,
    muon_key: Any = 0,
    apply: bool = True,
) -> CorrectedMuon:
    """Correct a ``(pt,eta,phi,mass)`` tuple, mapping, or attribute object."""
    pt = _component(muon, "pt", 0)
    eta = _component(muon, "eta", 1)
    phi = _component(muon, "phi", 2)
    mass = _component(muon, "mass", 3)
    factor, gaussian, region = correction_factor(
        eta, calibration, seed=seed, event_key=event_key,
        muon_key=muon_key, apply=apply,
    )
    # The fitted model changes curvature/momentum while preserving eta, phi and
    # the muon rest mass used in the four-vector reconstruction.
    return CorrectedMuon(pt * factor, eta, phi, mass, factor, region, gaussian)


def invariant_mass(first: CorrectedMuon, second: CorrectedMuon) -> float:
    def vector(muon: CorrectedMuon) -> tuple[float, float, float, float]:
        px = muon.pt * math.cos(muon.phi)
        py = muon.pt * math.sin(muon.phi)
        pz = muon.pt * math.sinh(muon.eta)
        energy = math.sqrt((muon.pt * math.cosh(muon.eta)) ** 2 + muon.mass**2)
        return px, py, pz, energy
    a, b = vector(first), vector(second)
    mass2 = (a[3] + b[3]) ** 2 - (a[0] + b[0]) ** 2 - (a[1] + b[1]) ** 2 - (a[2] + b[2]) ** 2
    return math.sqrt(max(0.0, mass2))


def correct_dimuon(
    first: Any,
    second: Any,
    calibration: Calibration,
    *,
    seed: int = 314159,
    event_key: Any = None,
    apply: bool = True,
) -> tuple[CorrectedMuon, CorrectedMuon, float]:
    """Correct two muons independently and return them with their invariant mass."""
    one = correct_muon(first, calibration, seed=seed, event_key=event_key,
                       muon_key=0, apply=apply)
    two = correct_muon(second, calibration, seed=seed, event_key=event_key,
                       muon_key=1, apply=apply)
    return one, two, invariant_mass(one, two)
