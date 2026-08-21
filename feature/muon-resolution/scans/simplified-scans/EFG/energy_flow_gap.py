#!/usr/bin/env python3
"""Energy Gap Flow observables from Mikael Mieskolainen's 2026-07-14 gist.

The numerical algorithm is preserved; this module adds imports, type handling,
and an event-level convenience wrapper suitable for the local analysis.
"""

from __future__ import annotations

import numpy as np

try:
    import numba
    njit = numba.njit
except ImportError:  # pragma: no cover - portable downstream fallback
    def njit(function):
        return function


@njit
def energy_gap_flow_score(eta, pt, side, eta_min, eta_max, pt_min, beta, q0):
    if side != -1 and side != 1:
        raise ValueError("energy_gap_flow_score: side must be -1 or +1")
    if eta_max <= eta_min:
        raise ValueError("energy_gap_flow_score: eta_max must exceed eta_min")
    if pt_min < 0.0:
        raise ValueError("energy_gap_flow_score: pt_min must be non-negative")
    if beta <= 0.0:
        raise ValueError("energy_gap_flow_score: beta must be positive")
    if q0 <= 0.0:
        raise ValueError("energy_gap_flow_score: q0 (GeV) must be positive")
    oriented = side * eta
    mask = (eta_min < oriented) & (oriented < eta_max) & (pt > pt_min)
    u = oriented[mask]; side_pt = pt[mask]
    if u.size == 0:
        return 1.0
    order = np.argsort(u); u = u[order]; side_pt = side_pt[order]
    tail_pt = np.sum(side_pt); previous_u = eta_min; numerator = 0.0
    for i in range(u.size):
        upper_weight = np.exp(beta * (u[i] - eta_max))
        lower_weight = np.exp(beta * (previous_u - eta_max))
        numerator += np.exp(-tail_pt / q0) * (upper_weight - lower_weight)
        tail_pt -= side_pt[i]; previous_u = u[i]
    numerator += 1.0 - np.exp(beta * (previous_u - eta_max))
    normalization = -np.expm1(-beta * (eta_max - eta_min))
    return min(1.0, max(0.0, numerator / normalization))


@njit
def energy_gap_flow_effective_gap(score, eta_min, eta_max, beta):
    if not np.isfinite(score):
        raise ValueError("energy_gap_flow_effective_gap: score must be finite")
    if score < 0.0 or score > 1.0:
        raise ValueError("energy_gap_flow_effective_gap: score must be in [0, 1]")
    if eta_max <= eta_min:
        raise ValueError("energy_gap_flow_effective_gap: eta_max must exceed eta_min")
    if beta <= 0.0:
        raise ValueError("energy_gap_flow_effective_gap: beta must be positive")
    length = eta_max - eta_min
    if score == 0.0:
        return 0.0
    if score == 1.0:
        return length
    normalization = -np.expm1(-beta * length)
    effective_gap = -np.log1p(-normalization * score) / beta
    return min(length, max(0.0, effective_gap))


@njit
def forward_gap_size(eta, pt, side, eta_min, eta_max, pt_min):
    if side != -1 and side != 1:
        raise ValueError("forward_gap_size: side must be -1 or +1")
    if eta_max <= eta_min:
        raise ValueError("forward_gap_size: eta_max must exceed eta_min")
    if pt_min < 0.0:
        raise ValueError("forward_gap_size: pt_min must be non-negative")
    if eta.size != pt.size:
        raise ValueError("forward_gap_size: eta and pt must have equal size")
    outermost_u = eta_min
    for i, value in enumerate(eta):
        oriented = side * value
        if eta_min < oriented < eta_max and pt[i] > pt_min:
            outermost_u = max(outermost_u, oriented)
    return eta_max - outermost_u


def event_observables(eta, pt, *, eta_min=-2.5, eta_max=2.5,
                      pt_min=0.2, beta=0.1, q0=2.0):
    eta = np.asarray(eta, dtype=np.float64); pt = np.asarray(pt, dtype=np.float64)
    backward = energy_gap_flow_score(eta, pt, -1, eta_min, eta_max, pt_min, beta, q0)
    forward = energy_gap_flow_score(eta, pt, +1, eta_min, eta_max, pt_min, beta, q0)
    score = max(backward, forward)
    traditional = max(forward_gap_size(eta, pt, -1, eta_min, eta_max, pt_min),
                      forward_gap_size(eta, pt, +1, eta_min, eta_max, pt_min))
    return {
        "efg_backward": float(backward), "efg_forward": float(forward),
        "efg_score": float(score),
        "efg_effective_gap": float(energy_gap_flow_effective_gap(score, eta_min, eta_max, beta)),
        "forward_gap": float(traditional),
        "gap_side": -1 if backward > forward else 1,
    }
