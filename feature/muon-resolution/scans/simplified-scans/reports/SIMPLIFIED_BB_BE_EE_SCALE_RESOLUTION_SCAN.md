# Simplified BB/BE/EE Muon Scale and Resolution Grid Scan

## Context and lineage

This is the next iteration of the signed 25-region and high-fidelity studies in
the parent `scans` directory. Those studies showed that several ordered signed
eta-pair regions were statistically fragile and that an identity point could be
preferred even when non-mass control distributions were poorly modelled.

This iteration reduces the calibration to detector-motivated per-muon regions:

- barrel: `|eta| < 1.4`;
- endcap: `1.4 <= |eta| < 2.5`.

The fit event categories are unordered BB, BE, and EE. In particular, BE and EB
from a leading/subleading split are combined. Corrections are still applied to
each muon independently according to its own absolute eta.

## Model and fit strategy

Four parameters are shared across all three mass histograms:

| Parameter | Meaning |
|---|---|
| `barrel_scale` | fractional barrel-muon momentum scale shift |
| `barrel_resolution` | added fractional Gaussian barrel width |
| `endcap_scale` | fractional endcap-muon momentum scale shift |
| `endcap_resolution` | added fractional Gaussian endcap width |

For simulated muon `i` in detector region `R`, the fitted model is

```text
pT_corrected(i) = pT(i) * [1 + scale(R) + resolution(R) * G(i)]
```

with independent standard-normal `G(i)` values. Three event-keyed seed replicas
are averaged at every hypothesis. The dimuon four-vector is reconstructed after
correcting both muons.

The simultaneous objective is the sum of the BB, BE, and EE single-source
analytic Barlow--Beeston profile deviances. Each category is independently
shape-normalized, so this fit does not calibrate the data/MC yield.

A brute-force 21-by-21-by-21-by-21 grid would be unnecessarily expensive.
Instead, `simplified_scan.py` uses an explicit staged grid search:

1. scan the two scale coordinates while holding the resolutions fixed;
2. scan the two resolution coordinates at the selected scales;
3. repeat the coarse block scans;
4. scan finer two-dimensional scale and resolution grids;
5. finish with a local `3^4` Cartesian grid allowing all coordinates to move.

This remains a grid search, but profiles correlated coordinate blocks instead of
evaluating every distant four-dimensional combination.

The established selection was retained: IsoMu24 trigger, two tight PF isolated
muons, `pT > 15 GeV`, leading `pT > 26 GeV`, `|eta| < 2.5`, opposite charge,
impact-parameter requirements, and nominal `70 < m_mumu < 110 GeV`. Adaptive
mass bins are finer across the resonance and merged in the tails. The fit uses
the full selected samples; no train/evaluation split was requested for this run.

## Inputs and execution

- Data: 31 usable ROOT files, 4,389,502 input entries.
- Simulation: 37 usable ROOT files, 1,914,000 input entries.
- One additional downloaded MC file,
  `14FF8F13-51E6-784E-BD14-B3C4417B3565.root`, was skipped because it did not
  contain `pfExtractor/pfTree`.
- Selected data: BB 68,008; BE 57,396; EE 18,436.
- Selected MC: BB 306,304; BE 259,701; EE 81,091.
- Event extraction: 432.1 seconds.
- Staged scan: 229.2 seconds.
- Total: 663.7 seconds (11.1 minutes).

The full output is in `plots/full_production_r3`. The exact configuration,
bin edges, cutflows, every evaluated grid point, and regional metrics are stored
in `summary.json`; the compact downstream constants are in `calibration.json`.

## Best grid point

| Detector region | Scale shift | Added resolution |
|---|---:|---:|
| Barrel | -0.125% | 0.375% |
| Endcap | -0.150% | 0.875% |

The minimum is not on any configured parameter boundary. The combined objective
decreased from 405.30 at the identity point to 186.21, an improvement of 219.09.

| Event category | Identity deviance / ndof | Corrected deviance / ndof | Improvement |
|---|---:|---:|---:|
| BB | 154.52 / 42 | 55.95 / 42 | 98.57 |
| BE | 195.15 / 42 | 91.37 / 42 | 103.78 |
| EE | 55.63 / 42 | 38.90 / 42 | 16.73 |

The correction clearly improves all three mass templates. EE has a satisfactory
aggregate corrected objective and BB is much improved. BE remains imperfect at
about 2.18 deviance units per nominal degree of freedom, with coherent residuals
visible in the ratio and pull panels. This is therefore a useful first shared
calibration, not evidence that all data/MC modelling issues have disappeared.

## Outputs and plotting

Each `BB`, `BE`, and `EE` subdirectory contains PNG, SVG, and PDF mass plots with:

1. data, uncorrected shape-normalized MC, and best corrected MC;
2. both data/MC ratios;
3. corrected signed Barlow--Beeston profile pulls.

The top-level fine scale and resolution landscape plots mark the selected joint
minimum. Histograms are preserved in `simplified_histograms.root`. The separate
`plot_simplified_validation.py` utility regenerates these plots without the
expensive ROOT event loop.

## Expanded muon, event, and Z validation

An additional full-sample validation pass was produced with
`plot_simplified_observables.py`. It reread all 31 data files and 37 usable MC
files once, built multi-replica corrected templates, and rendered each canvas in
an isolated ROOT process. Histogram production plus 33 plot renders took 481.7
seconds. The following observables are now available in each BB, BE, and EE
subdirectory:

- dimuon invariant mass;
- leading and subleading muon pT, eta, and phi;
- reconstructed vertex multiplicity;
- Z-candidate pT, rapidity, and phi.

Each three-panel PNG/SVG/PDF plot shows shape-normalized data, uncorrected MC,
best corrected MC, both data/MC ratios, and the corrected signed
Barlow--Beeston profile pulls. Pull panels now include labeled green +/-1 sigma,
orange +/-2 sigma, and red +/-3 sigma reference lines.

The exact metrics are stored in `observable_validation_summary.json`, with all
templates in `observable_histograms.root`. Selected results illustrate the
scope of the correction:

| Region/observable | Nominal deviance | Corrected deviance | Change |
|---|---:|---:|---:|
| BB mass | 154.52 | 55.95 | +98.57 improvement |
| BE mass | 195.15 | 91.37 | +103.78 improvement |
| EE mass | 55.63 | 38.90 | +16.73 improvement |
| EE subleading-muon pT | 263.66 | 100.97 | +162.69 improvement |
| BB leading-muon pT | 75.94 | 90.94 | 15.01 worse |
| BE leading-muon pT | 128.25 | 140.26 | 12.02 worse |
| BB Z pT | 79.55 | 83.38 | 3.84 worse |
| BE Z pT | 153.67 | 153.58 | essentially unchanged |

Variables that a momentum correction does not directly change, such as eta,
phi, and vertex multiplicity, remain essentially unchanged apart from small
event-migration effects at the mass-window boundary. Most notably, the vertex
deviances remain 6789/49 (BB), 5766/49 (BE), and 1708/49 (EE) after correction.
This confirms that pileup modelling is a separate, dominant problem and cannot
be repaired by muon momentum scale or resolution smearing.

The expanded suite therefore supports a deliberately narrow conclusion: the
fitted correction improves the Z-mass response, and it repairs a large EE
subleading-pT discrepancy, but it is not a general-purpose reweighting of event
or Z-boson kinematics. Several pT and Z-kinematic distributions are unchanged or
slightly worse and should be revisited after pileup and efficiency weighting.

## Downstream correction API

`muon_corrections.py` has no ROOT dependency. `Calibration.from_summary` reads
either `calibration.json` or the full scan summary. `correct_muon` and
`correct_dimuon` accept tuples, mappings, or objects exposing `pt`, `eta`, `phi`,
and `mass`. The API returns the applied factor, region, and random draw together
with the corrected momentum and reconstructed mass.

The RNG is deterministic per `(seed, event_key, muon_key)`, rather than a
sequential stream. A unique event key is required to prevent every event from
receiving the same smearing draw. The two muons use distinct muon keys.

## Interpretation and next steps

These plots use the same full sample that determined the parameters, so they are
fit diagnostics rather than independent validation. The next robust iteration
should use a deterministic holdout or statistically independent file split,
then inspect ordered BE and EB subsets separately. If those ordered validation
subsets disagree, the likely next dependency is muon `pT` or a selection/model
effect—not a leading/subleading detector correction.

The earlier pileup and Z-pT discrepancies also remain conceptually separate from
this mass calibration. More MC reduces template statistical uncertainty but
does not substitute for pileup, efficiency, or generator-kinematic reweighting.
