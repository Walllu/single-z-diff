# Signed Eta-Pair Muon Scale and Resolution Scan

**Iteration:** 2, following `FIRST_PASS_ZMUMU_RESOLUTION_STUDY.md`  
**Date:** 2026-07-14  
**Location:** `single-z-diff/muon-resolution/scans`

## Purpose and lineage

The first iteration established an inclusive, shape-normalized Z→μμ data/MC
comparison. This iteration adds ordered signed-eta categories and a reproducible
template scan of a common muon momentum-scale shift and additional simulation
resolution. It is an exploration framework, not yet a final calibrated
correction or uncertainty measurement.

## Signed eta-pair definition

The five configured regions are:

| Name | Eta interval |
|---|---|
| `neg_endcap` | `[-2.5, -2.1)` |
| `neg_transition` | `[-2.1, -1.4)` |
| `barrel` | `[-1.4, 1.4)` |
| `pos_transition` | `[1.4, 2.1)` |
| `pos_endcap` | `[2.1, 2.5]` |

The Cartesian product gives 25 ordered names such as
`neg_endcap__barrel`, where the first label is the nominal leading muon and the
second the nominal subleading muon. These are user-requested analysis bins; the
labels should not be interpreted as an official CMS subsystem-boundary
definition. In particular, standard CMS reconstructed-muon fiducial selections
often stop at `|eta| < 2.4`; this exploratory configuration intentionally uses
the requested track-motivated outer boundary of 2.5.

## Selection

All selection values are in the `SELECTION` dictionary near the script top.
The default requires:

- `IsoMu24` (retained as a configurable requirement even though the samples are
  expected to be triggered);
- at least two tight PF muons with `pT > 15 GeV` and `|eta| < 2.5`;
- delta-beta corrected relative isolation below 0.15;
- `|dxy| < 0.05 cm` and `|dz| < 0.10 cm` for prompt quality;
- an opposite-sign pair, choosing the pair closest to the nominal Z mass;
- nominal leading-muon `pT > 26 GeV`;
- nominal dimuon mass in 70–110 GeV.

The clean sample and its eta-pair category are fixed using nominal kinematics.
Scale and smearing are applied afterward. A shifted template is again restricted
to 70–110 GeV, so mass-window migration is represented while selection and eta
category migration are not.

## Scale and smearing model

For each simulated muon independently,

`pT' = pT * (1 + scale + resolution * G)`, with `G ~ N(0,1)`.

Both muons receive independent fixed Gaussian draws. Eta, phi, and muon mass
remain unchanged, after which the dimuon four-vector is reconstructed. The same
scale and added-resolution parameters are used for both legs at a grid point.
This is preferable to modifying only the leading leg: it avoids an asymmetric
template bias while keeping the first scan two-dimensional. A later fit can
assign separate per-eta parameters and solve the coupled regional system.

The coarse grid is configured in `SCAN_GRID`:

- scale: -5% to +5% in 1% steps;
- added resolution: 0% to 1% in 0.2% steps.

Four explicit seeds are defined in `RNG_SEEDS`; one is selected using
`--seed-name`. The same stored draws are reused across every grid point, avoiding
point-to-point fluctuations from newly generated random numbers.

## Histogram and objective configuration

`DEFAULT_HISTOGRAMS` defines all six comparisons: dimuon mass, leading and
subleading pT, both-muon eta, Z pT, and vertex multiplicity.
`REGION_HISTOGRAM_OVERRIDES` can replace any histogram definition for a named
eta-pair. The default mass fit uses 40 one-GeV bins over 70–110 GeV.

Two objectives are available:

- `chi2`: shape-normalized Pearson-style comparison using
  `sigma_data^2 + alpha^2 sigma_MC^2`, explicitly including finite MC
  statistical uncertainty;
- `poisson`: binned Poisson deviance for data conditional on the normalized MC
  expectation. This first-pass Poisson implementation does not profile finite
  template statistics (a Barlow–Beeston treatment would be a future upgrade).

One degree of freedom is subtracted for the fitted shape normalization. Complete
grid values, normalizations, accepted yields, and best point are saved to JSON.

## Modes and directed running

An unsmeared point for one category:

```bash
conda run -n mg5 python scan_zmumu_scale_resolution.py \
  --mode point --regions barrel__barrel --scale 0 --resolution 0 --metric chi2
```

A grid scan with the Poisson objective:

```bash
conda run -n mg5 python scan_zmumu_scale_resolution.py \
  --mode scan --regions barrel__barrel --metric poisson
```

Multiple region names can follow `--regions`; `--regions all` selects all 25.
Outputs include PNG/PDF plots, `histograms.root`, and `summary.json`. Scan mode
also creates a two-dimensional objective surface and plots at the best grid
point.

The later `--region-subdirectories` addition creates a configuration-labeled
top directory and one abbreviated directory per eta pair (`NE`, `NT`, `B`,
`PT`, and `PE`). This is intended for quickly reviewing bin populations and
deciding where region-specific histogram overrides are needed.

The feature was exercised on the complete usable inputs at scale 0 and added
resolution 0. The resulting `plots/scale_p0p000_resolution_p0p000` hierarchy
contains all 25 region directories and 300 plot files (six observables in PNG
and PDF for every region), plus the combined ROOT histogram file and JSON
summary. Six categories contain zero events in both samples after the nominal
selection: `NE_PT`, `NE_PE`, `NT_PE`, `PT_NE`, `PE_NE`, and `PE_NT`. Their
explicitly yield-labeled empty plots are retained because they are useful when
assessing whether the ordered region scheme is physically/statistically useful.

## Tests and produced outputs

The script passed Python compilation under the `mg5` conda environment with
PyROOT 6.40.02.

The requested default smoke test used one file and 10,000 entries per sample,
`barrel__barrel`, scale 0, resolution 0, and chi-square:

- 155 selected data events;
- 1,554 selected simulated events;
- chi-square 84.83 for 37 reported degrees of freedom;
- all six PNG/PDF comparisons, ROOT histograms, and JSON summary written to
  `smoke-tests/default_zero_point`.

A bounded Poisson grid test over all 66 configured points also completed. Its
smoke-sample minimum was scale 0 and added resolution 0.006; this is only a
software validation result and must not be interpreted as a calibration. It is
stored in `smoke-tests/poisson_grid`.

The corresponding full-sample zero-point plots are in `plots`. The full usable
inputs contained 865,402 data and 728,000 MC entries. In `barrel__barrel`, 13,256
data and 116,136 MC candidates passed, giving chi-square 56.90 for 39 reported
degrees of freedom. As in iteration 1, one MC file without `pfExtractor/pfTree`
was detected and skipped.

## Known next steps

The coarse 1% scale spacing is much wider than the precision normally sought
from a Z peak and is intended only to locate the neighborhood of the minimum.
The next iteration should refine the grid locally, study alternative seeds and
binning, add finite-template treatment to the Poisson likelihood, and formulate
a simultaneous fit of eta-dependent per-muon parameters across all 25 pair
categories. Closure tests on pseudo-data should precede use as a correction.

## Coarse 17-region chi-square scan test

An exclusion-list iteration added the top-level
`EXCLUDED_ETA_PAIR_ABBREVIATIONS` structure. Its default eight entries are
`NE_PT`, `NE_PE`, `NT_PT`, `NT_PE`, `PT_NE`, `PT_NT`, `PE_NE`, and `PE_NT`.
`--regions all` now scans the remaining 17 regions; the exclusion can be
temporarily bypassed with `--include-excluded-regions`.

The complete usable samples were scanned with the nominal seed, the configured
66-point coarse grid, the user's region-specific mass binning, shape
normalization, and the chi-square objective including MC statistical errors.

| Region | Data | MC | Scale | Added resolution | chi-square / ndof |
|---|---:|---:|---:|---:|---:|
| NE_NE | 149 | 1,292 | +0.010 | 0.004 | 32.557 / 19 |
| NE_NT | 375 | 3,223 | 0.000 | 0.002 | 20.895 / 19 |
| NE_B | 576 | 5,359 | 0.000 | 0.000 | 36.836 / 39 |
| NT_NE | 373 | 3,085 | 0.000 | 0.010 | 21.910 / 19 |
| NT_NT | 914 | 7,882 | 0.000 | 0.010 | 18.001 / 39 |
| NT_B | 2,207 | 18,838 | 0.000 | 0.004 | 45.471 / 39 |
| B_NE | 622 | 5,708 | 0.000 | 0.006 | 71.142 / 39 |
| B_NT | 2,251 | 19,556 | 0.000 | 0.002 | 64.430 / 39 |
| B_B | 13,256 | 116,136 | 0.000 | 0.000 | 56.899 / 39 |
| B_PT | 2,221 | 19,575 | 0.000 | 0.010 | 68.536 / 39 |
| B_PE | 607 | 5,648 | 0.000 | 0.010 | 44.743 / 39 |
| PT_B | 2,261 | 18,964 | 0.000 | 0.008 | 55.494 / 39 |
| PT_PT | 966 | 7,811 | 0.000 | 0.006 | 89.538 / 39 |
| PT_PE | 370 | 3,253 | 0.000 | 0.002 | 11.029 / 19 |
| PE_B | 596 | 5,287 | 0.000 | 0.000 | 120.218 / 39 |
| PE_PT | 357 | 3,090 | 0.000 | 0.002 | 19.665 / 19 |
| PE_PE | 141 | 1,336 | 0.000 | 0.004 | 22.546 / 19 |

The measured wall time was **217.1 seconds (3m37s)**: 127.7 seconds for the
single data/MC extraction and 89.4 seconds for all regional scans, best-point
histograms, and plotting. There were 17 × 66 = 1,122 grid evaluations and 238
PNG/PDF outputs. The scan inner loop fills only invariant mass; all six
diagnostics are filled once at the best point. Replacing repeated PyROOT vector
construction in that loop with the algebraically identical scalar invariant-
mass expression was validated over 1,000 random cases to a maximum absolute
difference of `2.3e-12 GeV` and reproduced the bounded reference scan exactly.

Outputs are stored under `plots/scan_chi2_seed_nominal`, including one directory
per region, a scan surface, best-point plots, ROOT histograms, and the full JSON
grid. These are exploratory minima. The 1% scale spacing is very coarse, and
`NT_NE`, `NT_NT`, `B_PT`, and `B_PE` minimize at the 1% resolution boundary;
their resolution values are therefore range-limited rather than measured
interior minima. `NE_NE` is the only category preferring a nonzero coarse scale
point (+1%). A finer, potentially wider scan and closure tests are required
before interpreting any entry as a correction.

## Expanded-data fine-scale/wide-resolution scan

After expanding the local 2016H subset from 6 to 31 files, the input schema and
file metadata were checked again for weights. Neither data nor simulation
contains an event/generator weight, pileup weight, cross section, luminosity
normalization object, or equivalent metadata. The branch named `lumi` is the
luminosity-section identifier, not a weight. The scan therefore continues to
make shape-only comparisons by normalizing each MC mass template to its data
integral.

The grid was changed to 21 scale points from -0.5% to +0.5% in 0.05% steps and
21 added-resolution points from 0% to 5% in 0.25% steps. This gives 441 points
per region and 7,497 evaluations across the 17 retained regions. The scan used
4,389,502 entries from all 31 data files and 728,000 entries from the 11 usable
MC files. Best points were:

| Region | Data | MC | Scale | Added resolution | chi-square / ndof |
|---|---:|---:|---:|---:|---:|
| NE_NE | 740 | 1,292 | +0.0045 | 0.0250 | 12.301 / 19 |
| NE_NT | 1,868 | 3,223 | -0.0005 | 0.0075 | 18.044 / 19 |
| NE_B | 3,018 | 5,359 | +0.0010 | 0.0200 | 38.265 / 39 |
| NT_NE | 1,844 | 3,085 | -0.0015 | 0.0000 | 18.061 / 19 |
| NT_NT | 4,805 | 7,882 | -0.0005 | 0.0100 | 35.944 / 39 |
| NT_B | 11,101 | 18,838 | -0.0005 | 0.0000 | 44.376 / 39 |
| B_NE | 3,191 | 5,708 | -0.0020 | 0.0175 | 31.057 / 39 |
| B_NT | 11,518 | 19,556 | -0.0025 | 0.0025 | 36.465 / 39 |
| B_B | 68,008 | 116,136 | -0.0010 | 0.0025 | 44.621 / 39 |
| B_PT | 11,365 | 19,575 | -0.0015 | 0.0025 | 48.859 / 39 |
| B_PE | 3,269 | 5,648 | -0.0015 | 0.0075 | 25.237 / 39 |
| PT_B | 10,944 | 18,964 | -0.0010 | 0.0075 | 37.121 / 39 |
| PT_PT | 4,682 | 7,811 | -0.0025 | 0.0075 | 47.196 / 39 |
| PT_PE | 1,900 | 3,253 | -0.0020 | 0.0050 | 16.286 / 19 |
| PE_B | 2,990 | 5,287 | -0.0020 | 0.0150 | 31.756 / 39 |
| PE_PT | 1,848 | 3,090 | -0.0010 | 0.0100 | 13.271 / 19 |
| PE_PE | 741 | 1,336 | +0.0005 | 0.0225 | 7.153 / 19 |

The measured total runtime was **519.5 seconds (8m40s)**: 262.4 seconds for
event extraction and 257.1 seconds for the 7,497 grid evaluations, best-point
histograms, ROOT/JSON serialization, and 238 plot files. Results are in
`plots/scan_chi2_seed_nominal_s21_m0p005_p0p005_r21_p0p000_p0p050`.

No best point reaches the 5% upper resolution boundary, so widening that axis
resolved the boundary problem seen in the coarse test. `NT_NE` and `NT_B`
prefer zero additional smearing, the physical lower boundary. No scale minimum
is exactly at ±0.5%; `NE_NE` at +0.45% lies close enough to the upper edge that
its local neighborhood should be checked in a follow-up. These remain separate
pair-region template fits with a common correction applied to both muons, not a
globally consistent solution for per-muon eta corrections.

## Poisson-likelihood cross-check and binning utility

The expanded-data 21×21 scan was repeated with the binned Poisson deviance,
holding inputs, selections, exclusions, histogram definitions, normalization,
RNG seed, and grid fixed. The Poisson implementation conditions on the
normalized MC template and does not include Barlow–Beeston finite-template
nuisances; the chi-square comparison does include MC bin errors.

| Region | Chi-square scale/res. | Poisson scale/res. | Poisson deviance / ndof |
|---|---:|---:|---:|
| NE_NE | +0.0045 / 0.0250 | +0.0045 / 0.0250 | 21.696 / 19 |
| NE_NT | -0.0005 / 0.0075 | -0.0005 / 0.0075 | 32.257 / 19 |
| NE_B | +0.0010 / 0.0200 | +0.0010 / 0.0200 | 60.550 / 39 |
| NT_NE | -0.0015 / 0.0000 | -0.0015 / 0.0000 | 30.394 / 19 |
| NT_NT | -0.0005 / 0.0100 | -0.0005 / 0.0100 | 58.554 / 39 |
| NT_B | -0.0005 / 0.0000 | -0.0005 / 0.0000 | 73.122 / 39 |
| B_NE | -0.0020 / 0.0175 | -0.0020 / 0.0175 | 45.431 / 39 |
| B_NT | -0.0025 / 0.0025 | -0.0025 / 0.0025 | 60.783 / 39 |
| B_B | -0.0010 / 0.0025 | -0.0010 / 0.0025 | 72.497 / 39 |
| B_PT | -0.0015 / 0.0025 | -0.0015 / 0.0025 | 78.570 / 39 |
| B_PE | -0.0015 / 0.0075 | -0.0010 / 0.0100 | 38.366 / 39 |
| PT_B | -0.0010 / 0.0075 | -0.0010 / 0.0075 | 61.904 / 39 |
| PT_PT | -0.0025 / 0.0075 | -0.0025 / 0.0075 | 74.590 / 39 |
| PT_PE | -0.0020 / 0.0050 | -0.0010 / 0.0075 | 28.271 / 19 |
| PE_B | -0.0020 / 0.0150 | -0.0020 / 0.0150 | 49.473 / 39 |
| PE_PT | -0.0010 / 0.0100 | -0.0010 / 0.0100 | 21.238 / 19 |
| PE_PE | +0.0005 / 0.0225 | +0.0005 / 0.0225 | 9.951 / 19 |

Fifteen of 17 minima are identical. `B_PE` moves by one scale and one
resolution step; `PT_PE` moves by two scale and one resolution step. The result
is therefore robust to the objectives at the current granularity, with local
sensitivity in two lower-statistics mixed regions. The Poisson run took 526.1
seconds (8m46s), including 261.0 seconds of extraction. Its 238 plots and full
grid are in
`plots/scan_poisson_seed_nominal_s21_m0p005_p0p005_r21_p0p000_p0p050`.

The standalone `optimize_mass_binning.py` utility was also added. It reports
Freedman–Diaconis, Scott, and Sturges bin counts, uses Freedman–Diaconis as its
primary density-estimation rule, and applies a soft average-population cap plus
configurable minimum/maximum bin counts. It uses selected data only and never
edits the scan configuration automatically.

A three-file test wrote current-versus-recommended mass plots and JSON to
`binning-tests/three-data-files`. The recommended counts were: NE_NE 5, NE_NT
8, NE_B 13, NT_NE 8, NT_NT 20, NT_B 54, B_NE 15, B_NT 52, B_B 100 (configured
cap), B_PT 53, B_PE 15, PT_B 52, PT_PT 22, PT_PE 9, PE_B 14, PE_PT 9, and
PE_PE 5. This demonstrates why “optimal” is an operating-point diagnostic:
sparse endcaps favor fewer bins, while density rules favor very fine binning in
the barrel. Final choices should also be checked for fit stability and
pseudo-experiment coverage.

## Tiered-v1 Poisson rerun and correction interface

The user-selected first-attempt final binning was labeled `tiered_v1`: 80 mass
bins by default, 40 in the medium-statistics pair regions, and 20 in the sparse
endcap combinations. The complete 31-file Poisson scan was rerun without
changing its selection, exclusions, seed, or 21×21 scale/resolution grid.

| Region | Mass bins | Best scale | Best resolution | Deviance / ndof |
|---|---:|---:|---:|---:|
| NE_NE | 20 | +0.0045 | 0.0250 | 21.696 / 19 |
| NE_NT | 20 | -0.0005 | 0.0075 | 32.257 / 19 |
| NE_B | 40 | +0.0010 | 0.0200 | 60.550 / 39 |
| NT_NE | 20 | -0.0015 | 0.0000 | 30.394 / 19 |
| NT_NT | 40 | -0.0005 | 0.0100 | 58.554 / 39 |
| NT_B | 80 | -0.0010 | 0.0025 | 126.300 / 79 |
| B_NE | 40 | -0.0020 | 0.0175 | 45.431 / 39 |
| B_NT | 80 | -0.0030 | 0.0000 | 138.963 / 79 |
| B_B | 80 | -0.0010 | 0.0025 | 134.866 / 79 |
| B_PT | 80 | -0.0010 | 0.0100 | 116.278 / 79 |
| B_PE | 40 | -0.0010 | 0.0100 | 38.366 / 39 |
| PT_B | 80 | -0.0010 | 0.0000 | 120.644 / 79 |
| PT_PT | 40 | -0.0025 | 0.0075 | 74.590 / 39 |
| PT_PE | 20 | -0.0010 | 0.0075 | 28.271 / 19 |
| PE_B | 40 | -0.0020 | 0.0150 | 49.473 / 39 |
| PE_PT | 20 | -0.0010 | 0.0100 | 21.238 / 19 |
| PE_PE | 20 | +0.0005 | 0.0225 | 9.951 / 19 |

Thirteen of 17 minima are unchanged from the preceding Poisson binning. The
four changes are all newly 80-bin categories: `NT_B`, `B_NT`, `B_PT`, and
`PT_B`. This quantifies the binning dependence rather than assuming it is zero.
The run took 511.0 seconds (8m31s), including 254.9 seconds of extraction.
Results are stored in
`plots/scan_poisson_seed_nominal_s21_m0p005_p0p005_r21_p0p000_p0p050_bins_tiered_v1`.

Each regional scan plot now displays delta deviance with the visual scale capped
at 30, black contours at 2.30, 6.18, and 11.83, and a solid red minimum marker;
the underlying JSON retains every uncapped objective value. The top-level
`best_fit_parameter_grid` is a two-panel 5×5 eta-pair map of best scale and
resolution in percent. It displays fitted zeros explicitly and masks all eight
excluded regions with gray `X` cells.

`muon_pair_correction.py` supplies the reusable `Muon` data class and
`correct_muon_pair` function. It defaults to this tiered-v1 summary, caches the
large JSON after first access, orders the two nominal input muons by transverse
momentum, finds their signed eta-pair region, and applies
`factor = 1 + scale + resolution*G` independently to both simulated muons. A
stable `event_key` combined with the selected `rng_seed` provides reproducible
event-specific Gaussian draws. `apply_correction=False` reports proposed
factors without changing momenta; excluded or out-of-acceptance pairs return
identity factors. The function must not be used to smear collision data.

## Multi-seed correction validation

The fitted correction was subsequently validated against uncorrected MC using
four deterministic event-keyed RNG seeds. The new validator, its full report,
204 plots across the 17 fitted eta-pair regions, ROOT histograms, and numerical
summary are under `scans/scan-validation`. Across the mass comparisons, 52 of
68 seed-region cases improve on uncorrected MC and the seed-mean deviance
improves in 14 of 17 regions. The validation also identifies several
seed-sensitive or non-improving regions and documents that the original scan's
sequential RNG convention prevents exact template closure with the reusable
event-keyed correction. This is a robustness pass rather than a final closure
claim.
