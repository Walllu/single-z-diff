# Muon Scale and Resolution Correction Validation

## Purpose and lineage

This is the first validation pass for the signed eta-pair correction derived in
the `tiered_v1` Poisson scan. It follows the original Z-to-dimuon event
selection, signed eta-pair scan, expanded-data fine scan, Poisson cross-check,
and user-selected tiered mass binning documented in
`scans/reports/SIGNED_ETA_PAIR_SCALE_RESOLUTION_SCAN.md`.

The goal here is to compare collision data with both the nominal uncorrected
simulation and corrected simulation, and to expose the dependence of the
stochastic correction on its random-number seed. This is a shape-only
validation: every simulated histogram is independently normalized to the data
integral.

## Implementation

`validate_muon_corrections.py` reads the same ROOT inputs and uses the same
selection and signed eta-pair definitions as `scan_zmumu_scale_resolution.py`.
It loads the best scale and added-resolution values from:

```text
scans/plots/scan_poisson_seed_nominal_s21_m0p005_p0p005_r21_p0p000_p0p050_bins_tiered_v1/summary.json
```

For every retained eta-pair region, it makes the six existing distributions:

- dimuon mass;
- leading- and subleading-muon transverse momentum;
- muon pseudorapidity;
- dimuon transverse momentum;
- reconstructed-vertex multiplicity.

The upper panel overlays data, uncorrected MC, and corrected MC for each seed.
The lower panel overlays the corresponding data/MC ratios. Both PNG and PDF
files are written beneath the abbreviated region directories (`NE_NE`,
`B_B`, `PT_PE`, and so on). The eight scan-excluded opposite-side transition/
endcap combinations have no fitted correction and are therefore not plotted.

The correction is applied only to simulation. Its multiplicative factor for
each muon is

```text
1 + fitted_scale + fitted_resolution * Gaussian(0, 1).
```

The two muons receive independent Gaussian values. A stable event identifier,
the chosen seed, and the muon index determine the draws, so rerunning the script
or changing file order does not change a given event's result.

## Full-run configuration

- Data: all 31 locally available 2016H files, 4,389,502 entries.
- Simulation: all 11 usable Z-to-dimuon files, 728,000 entries.
- Seeds: `314159`, `271828`, `161803`, and `141421`.
- Calibrated eta-pair regions: 17.
- Fit/validation mass range: 70--110 GeV.
- Output: `plots/tiered_v1_seeds_4`.
- Runtime: 620.7 seconds (10 minutes 21 seconds).
- Products: 204 plots, `validation_histograms.root`, and
  `validation_summary.json`.

No zero-length output files were found. A smaller two-seed, one-region test is
stored in `smoke-tests/tiered_v1_seeds_2`.

## Dimuon-mass validation

The table reports the uncorrected binned Poisson deviance, the mean corrected
deviance across the four seeds, the best corrected deviance and its seed, and
the number of corrected seeds that improve on the uncorrected template. These
numbers are validation diagnostics, not refitted parameters.

| Region | Uncorrected | Corrected mean | Corrected best | Best seed | Seeds improved |
|---|---:|---:|---:|---:|---:|
| NE_NE | 321.642 | 105.605 | 30.477 | 141421 | 4/4 |
| NE_NT | 46.413 | 39.857 | 36.573 | 314159 | 4/4 |
| NE_B | 84.019 | 81.192 | 70.099 | 161803 | 2/4 |
| NT_NE | 33.957 | 30.394 | 30.394 | 314159 | 4/4 |
| NT_NT | 103.117 | 77.261 | 63.190 | 141421 | 4/4 |
| NT_B | 142.650 | 147.218 | 138.867 | 141421 | 2/4 |
| B_NE | 74.842 | 70.492 | 61.195 | 161803 | 2/4 |
| B_NT | 205.422 | 138.963 | 138.963 | 314159 | 4/4 |
| B_B | 230.888 | 135.319 | 128.471 | 141421 | 4/4 |
| B_PT | 169.032 | 149.966 | 118.977 | 161803 | 2/4 |
| B_PE | 67.827 | 54.269 | 44.534 | 314159 | 4/4 |
| PT_B | 142.745 | 120.644 | 120.644 | 314159 | 4/4 |
| PT_PT | 102.472 | 75.676 | 66.300 | 161803 | 4/4 |
| PT_PE | 34.639 | 41.804 | 31.830 | 314159 | 1/4 |
| PE_B | 69.848 | 61.017 | 56.560 | 161803 | 4/4 |
| PE_PT | 21.381 | 28.040 | 22.423 | 141421 | 0/4 |
| PE_PE | 26.180 | 20.875 | 16.413 | 141421 | 3/4 |

Across the complete comparison, 52 of 68 seed-region combinations improve the
mass deviance. The mean corrected deviance is better than the uncorrected value
in 14 of 17 regions. Improvements are especially clear in `B_B`, `B_NT`, and
`NE_NE`; however, the very large `NE_NE` spread also demonstrates that a sparse
region with a sizeable fitted resolution is strongly seed-dependent.

The three regions whose seed-mean deviance worsens are `NT_B`, `PT_PE`, and
`PE_PT`. `PE_PT` improves for none of the four seeds, and `PT_PE` improves for
only one. These should not be treated as validated corrections without further
study. Regions whose fitted resolution is zero are seed-independent, as
expected, because their correction is purely a deterministic scale shift.

## Important reproducibility caveat

The grid scan used one sequential random-number stream while iterating through
its in-memory events. The reusable correction and this validation instead use
event-keyed random values, which are stable under reruns, file reordering, and
event filtering. Consequently, even the nominal seed does not reproduce the
exact random template used to choose the scan minimum. That explains why a
validation deviance can differ materially from the stored best-point scan
deviance.

Event-keyed random numbers are the safer analysis interface, but an exact
closure test should refit the grid with the same event-keyed RNG convention.
The present plots are therefore useful as a robustness diagnostic, not yet a
final closure claim. A later calibration pass should also consider independent
fit and validation samples to avoid assessing the correction only on the data
used to derive it.

## Running the validator

From `scans/scan-validation` with the `mg5` environment active:

```bash
python validate_muon_corrections.py
```

Regions, seeds, input limits, parameter JSON, and output location can be
directed through command-line options. For example:

```bash
python validate_muon_corrections.py \
  --regions barrel__barrel \
  --seeds 314159 271828 \
  --max-files 1 --max-events 5000 \
  --output-dir smoke-tests
```
