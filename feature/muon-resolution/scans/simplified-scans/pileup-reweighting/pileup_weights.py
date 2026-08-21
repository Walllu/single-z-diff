#!/usr/bin/env python3
"""Reusable reconstructed-vertex event weights for BB/BE/EE simulation."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class VertexWeightGroup:
    minimum: int
    maximum: int
    weight: float
    data_count: float = 0.0
    mc_count: float = 0.0


class PileupCalibration:
    """Regional lookup for approximate pileup weights derived from nVertices."""

    def __init__(self, groups: dict[str, list[VertexWeightGroup]], overflow_vertex: int = 50):
        self.groups = groups
        self.overflow_vertex = int(overflow_vertex)

    @classmethod
    def from_json(cls, filename: str | Path) -> "PileupCalibration":
        payload = json.loads(Path(filename).read_text())
        groups = {
            region: [VertexWeightGroup(**item) for item in values]
            for region, values in payload["regions"].items()
        }
        return cls(groups, payload["configuration"]["overflow_vertex"])

    def weight(self, region: str, nvertices: int, *, apply: bool = True) -> float:
        if not apply:
            return 1.0
        if region not in self.groups:
            raise KeyError(f"Unknown pileup-weight region {region!r}")
        value = min(max(0, int(nvertices)), self.overflow_vertex)
        for group in self.groups[region]:
            if group.minimum <= value <= group.maximum:
                return group.weight
        raise RuntimeError(f"No {region} pileup-weight group covers nVertices={value}")

