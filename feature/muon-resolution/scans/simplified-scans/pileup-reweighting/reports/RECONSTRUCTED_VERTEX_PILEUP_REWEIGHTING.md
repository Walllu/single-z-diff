# Regional Reconstructed-Vertex Pileup Reweighting

## Purpose and lineage

The simplified BB/BE/EE muon calibration exposed a very large disagreement in
the number of reconstructed vertices. Muon momentum scale and resolution
smearing cannot change the underlying pileup profile. This follow-up therefore
derives regional event weights from `nVertices` and evaluates their effects on
the already muon-corrected simulation.

The input ntuples do not store a true simulated interaction count or an official
luminosity-derived data pileup profile. This study is consequently an
approximate, selection-level reconstructed-vertex reweighting. It should not be
presented as the official CMS pileup prescription.

## Definition

For event category `R` and adaptive vertex group `g`, the event weight is

```text
w_R(g) = [N_data,R(g) / N_data,R(total)]
         ---------------------------------
         [N_MC,R(g)   / N_MC,R(total)]
```

Weights are derived separately for unordered BB, BE, and EE categories after
the established trigger, muon-quality, isolation, charge, pT, eta, and nominal
70--110 GeV mass selection. This region dependence follows the requested study
but can absorb category-specific selection differences in addition to pileup.

Raw per-vertex ratios fluctuate in the high-vertex tails. The production method
therefore builds contiguous groups requiring at least 100 data and 300 MC events.
The final deficient tail is merged backward. `nVertices >= 50` is mapped to the
overflow group. No arbitrary weight clipping is applied.

The weights are attached once per physical MC event. Three event-keyed muon
smearing replicas are averaged as before. Weighted template variances use the
squared per-event contribution, so Barlow--Beeston metrics reflect the reduced
effective MC statistics.

## Production inputs and weights

The run used all 31 data files and 37 usable MC files: 143,840 selected data and
647,096 selected MC events across the three regions. One file without the
required tree was skipped consistently with the calibration. Total runtime,
including 36 isolated canvas renders, was 486.2 seconds.

| Region | Groups | Weight range | Final tail group | Nominal MC | Effective weighted MC |
|---|---:|---:|---:|---:|---:|
| BB | 32 | 0.324--3.418 | 35--50 | 306,304 | 270,935 |
| BE | 32 | 0.311--3.522 | 36--50 | 259,701 | 229,343 |
| EE | 27 | 0.361--2.924 | 33--50 | 81,091 | 72,423 |

The sum of weights equals the nominal MC event count within floating-point
precision in every region. Reweighting retains roughly 88--89% of the original
effective MC statistics.

## Effects on validation distributions

All comparisons start from simulation with the fitted barrel/endcap muon scale
and resolution correction already applied. `before` means no vertex weight;
`after` adds the regional reconstructed-vertex weight.

| Region/observable | Before deviance | After deviance | Improvement |
|---|---:|---:|---:|
| BB nVertices | 6789.13 | 300.49 | 6488.64 |
| BB leading-muon pT | 90.94 | 86.75 | 4.19 |
| BB subleading-muon pT | 85.32 | 79.09 | 6.22 |
| BB Z pT | 83.38 | 76.08 | 7.31 |
| BB mass | 55.95 | 53.09 | 2.86 |
| BE nVertices | 5766.00 | 246.86 | 5519.15 |
| BE leading-muon pT | 140.26 | 122.71 | 17.55 |
| BE Z pT | 153.58 | 142.80 | 10.77 |
| BE mass | 91.37 | 89.98 | 1.39 |
| EE nVertices | 1707.66 | 68.12 | 1639.54 |
| EE leading-muon pT | 93.21 | 90.10 | 3.11 |
| EE subleading-muon pT | 100.97 | 95.85 | 5.12 |
| EE Z pT | 97.50 | 91.07 | 6.42 |
| EE mass | 38.90 | 40.17 | -1.28 |

The vertex objective does not reach zero because weights are constant within
adaptive groups rather than tuned independently in every bin, and corrected
mass-window migration differs slightly from the nominal sample used to derive
the groups. Its dramatic improvement is still partly closure by construction.

The more informative result is the correlated improvement in several pT
spectra. BE leading-muon pT and Z pT improve most clearly. BB and EE show smaller
but consistent gains in multiple pT observables. Some distributions worsen:
for example BE subleading pT changes by -1.16 deviance units and EE mass by
-1.28. The weights are therefore not a universal shape correction.

## Outputs

Each BB/BE/EE plot directory contains a regional weight curve plus 11 PNG/SVG/PDF
comparisons: dimuon mass; leading/subleading muon pT, eta, and phi; nVertices;
and Z-candidate pT, rapidity, and phi. Panels show data, muon-corrected MC before
pileup weighting, MC after weighting, both data/MC ratios, and weighted profile
pulls with +/-1, +/-2, and +/-3 sigma guides.

The exact adaptive groups and reusable constants are in `pileup_weights.json`.
`pileup_validation_summary.json` records every objective and effective event
count, and `pileup_histograms.root` permits plot-only regeneration.

## Limitations and recommended next step

The weights are derived and evaluated on the same selected data, so the vertex
closure is not independent. Reconstructed vertex count also includes vertex
reconstruction efficiency, hard-process selection, and detector effects. A
more defensible production correction would use the simulated true pileup count
and a luminosity-derived data pileup distribution, varied for the inelastic
cross-section uncertainty.

Before adopting these constants downstream, derive them on one file subset and
validate on another, compare regional weights against a single global weight,
and rerun the muon scale/resolution fit after pileup weighting. That last step is
important because the present muon correction was fitted before the pileup
profile was changed.

