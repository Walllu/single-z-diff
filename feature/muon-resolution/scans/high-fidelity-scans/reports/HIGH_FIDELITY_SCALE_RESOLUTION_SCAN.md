# High-Fidelity Signed Eta-Pair Muon Scale and Resolution Scan

## Context and lineage

This workflow is the next iteration of the 2016H CMS OpenData Z-to-dimuon
shape calibration. The original scripts established the event/object
selection, 25 ordered signed leading/subleading eta-pair categories, eight
excluded low-population opposite-side categories, fine scale/resolution grids,
chi-square and Poisson scans, tiered mass binning, and a reusable simulation-
only correction. A subsequent multi-seed validation exposed two limitations:
the scan's sequential random-number stream did not match the event-keyed
analysis correction, and fixed-template Poisson deviances could be dominated
by empty MC tail bins.

The code in `scans/high-fidelity-scans` is separate so all preceding outputs
remain reproducible. This report covers the first representative-region test,
not a final 17-region production calibration.

## Implemented workflow

### Event-keyed random numbers and seed averaging

Every Gaussian pair is derived from the global seed and stable
`run:lumi:event` key. A given event therefore receives the same Gaussian values
under reruns, file reordering, or earlier filtering. Leading and subleading
muons receive independent draws. The same draw is reused at every grid point,
which supplies common random numbers and keeps the likelihood landscape
correlated and smooth.

Each grid template averages a configurable number of seeds. The default is
three so exploratory scans remain quick. Importantly, seed replicas are used
to integrate the stochastic smearing kernel, not counted as new simulated
events. For each physical event the code first calculates its fractional bin
occupancy across replicas and then adds the square of that fraction to the MC
variance. This preserves the physical-event statistical unit and avoids an
artificial `1/sqrt(number of replicas)` reduction in MC uncertainty.

### Resonance-aware dynamic binning

The finest candidate widths are 2 GeV in the 70--82 and 100--110 GeV tails,
1 GeV in 82--84 and 98--100 GeV, and 0.5 GeV across 84--98 GeV. Adjacent bins
are merged outward-to-inward until the training subset contains at least ten
data events and twenty nominal MC events per bin by default. A deficient final
tail is merged into the preceding bin. Binning is derived only from training
events and then frozen for grid fitting and independent evaluation.

This directly removes the earlier `NE_NE` failure mode in which five data
events opposite an empty 108--110 GeV MC bin contributed 282 of 322 Poisson-
deviance units. All five representative fits have zero data-populated/MC-empty
bins at both the identity and fitted points.

### Barlow--Beeston MC statistical treatment

There is one simulated process and one profiled MC-statistical nuisance per
bin. For fractional replica templates, bin content `m` and physical-event
variance `v` are represented by the effective auxiliary count

```text
m_eff = m^2 / v.
```

For a fixed shape normalization `alpha`, the auxiliary-to-data exposure is
`tau = m_eff/(alpha*m)`. The profiled bin mean is

```text
mu_hat = (data + m_eff)/(1 + tau),
```

and the objective is the sum of the Poisson deviances for data about `mu_hat`
and effective MC about `tau*mu_hat`. With one MC source, this analytic
construction is equivalent to profiling the full per-bin Barlow--Beeston
nuisance. The overall comparison remains shape-only because `alpha` is fixed
from the data/MC integrals; one degree of freedom is removed for that
normalization.

### Independent deterministic split

Data and MC are each assigned independently to training or evaluation with a
SHA256 threshold of the sample label, split seed, and event key. The default is
80% training and 20% evaluation. The grid minimum is selected using training
events only. The identity point and frozen fitted point are then evaluated on
the untouched 20% subsets, both for the averaged template and every individual
seed.

This is more informative than validating on the fit sample, but it does not
create new MC information: both subsets are smaller than the complete sample.
For a final result, repeated folds or additional MC would improve stability.

### Correction-necessity utility

`assess_correction_need.py` compares the identity and fitted hypotheses using
three descriptive requirements:

- training improvement greater than four deviance units, an AIC-like `2k`
  cost for two fitted parameters;
- positive improvement on the independent evaluation subset;
- improvement for at least 60% of individual evaluation seeds.

A correction is `supported` only if all three pass. It is `prefer_identity` if
both independent evaluation and seed robustness fail; other cases are
`inconclusive`. These labels are decision aids, not p-values: resolution is
bounded at zero, parameters come from a grid, and the hypotheses are selected
after a scan.

### Validation plots

`plot_high_fidelity_validation.py` uses the frozen evaluation subset and makes
dimuon mass, leading/subleading muon transverse momentum, muon eta, dimuon
transverse momentum, and vertex multiplicity plots. The mass plot uses the
exact variable edges from the fit. Each plot contains:

1. data, uncorrected MC, corrected seed mean, individual seed lines, and the
   corrected seed envelope;
2. data/uncorrected-MC points and the corrected ratio envelope;
3. signed profiled Barlow--Beeston pulls for uncorrected MC and the corrected
   seed envelope.

Each MC shape is independently normalized to the data integral. Eta and vertex
multiplicity are not directly changed by a momentum-only correction, but are
retained to reveal selection-migration side effects.

The validator also has a `--render-only` mode that redraws from its saved ROOT
histogram archive. Seed-independent corrections are drawn as one central curve
without a zero-width envelope.

For visual studies where the 20% evaluation subset is too sparse, `--subset
all` recombines the deterministic partitions and writes a separately labeled
`validation-all` suite. These plots have better statistical precision but are
full-sample diagnostics, not independent closure tests, because they include
the training events used to select the correction.

Validation figures are written as PNG, SVG, and PDF. SVG/PDF are the preferred
inspection formats on macOS because PyROOT's bitmap backend can intermittently
drop glyphs even when the underlying vector canvas and ROOT histograms are
complete.

## Tests performed

All commands used `/opt/homebrew/Caskroom/miniforge/base/envs/mg5/bin/python`.
Syntax and CLI tests passed for all three scripts. Numerical checks established
that event-keyed draws reproduce exactly for the same key and change with
either key or seed, an exactly proportional data/MC template has zero profiled
deviance, and the vectorized invariant mass agrees with the scalar convention.

A one-file, two-seed, 5x5 smoke scan covered `NE_NE`, `NT_B`, and `B_B` and is
stored in `smoke-tests/core_smoke`. The necessity and all-observable plotting
utilities completed on that output. A separate truncated plot-style check is
stored in `smoke-tests/plot_style_test`.

## Representative all-input test

The main test used all 31 data files and all 11 usable MC files, an 80/20
split, seeds `314159`, `271828`, and `161803`, and the 11x11 test grid. It
covered sparse `NE_NE`, previously questionable `NT_B`, `PT_PE`, and `PE_PT`,
and high-statistics stable `B_B`. Event extraction and splitting took 267.3
seconds; the 605 multi-seed grid templates took only 12.4 seconds; total scan
time was 279.8 seconds.

| Region | Adaptive bins | Best scale | Best resolution | Training delta D | Evaluation delta D | Seeds improved | Assessment |
|---|---:|---:|---:|---:|---:|---:|---|
| NE_NE | 27 | +0.004 | 0.010 | +8.007 | +1.042 | 1/3 | inconclusive |
| NT_B | 43 | -0.001 | 0.000 | +4.104 | +4.397 | 3/3 | supported |
| PT_PE | 38 | 0.000 | 0.000 | 0.000 | 0.000 | 0/3 | prefer identity |
| PE_PT | 38 | -0.001 | 0.010 | +5.042 | -1.737 | 1/3 | prefer identity |
| B_B | 43 | -0.001 | 0.005 | +57.926 | +6.468 | 2/3 | supported |

Positive delta D means the fitted correction has the smaller profiled
deviance. Minimum training-bin populations were 13 data/20 MC for `NE_NE`,
26/41 for `NT_B`, 10/21 for `PT_PE`, 11/20 for `PE_PT`, and 143/233 for `B_B`.
No tested identity or best-point template contained a data-populated empty MC
bin.

The new treatment materially changes the interpretation of the earlier
problem regions:

- `NT_B` selects a pure -0.1% scale correction and zero added resolution. Its
  improvement survives the independent split and is seed-independent.
- `PT_PE` chooses the exact identity grid point, providing no evidence that a
  correction is needed in this test.
- `PE_PT` improves training but worsens independent evaluation, which is the
  expected signature of a fluctuation-driven minimum.
- `NE_NE` no longer has an empty-tail-bin pathology. Its averaged evaluation
  changes only modestly and individual seeds are not robust, so it remains
  inconclusive rather than showing the dramatic artificial improvement of the
  fixed-template scan.
- `B_B` is strongly preferred in training and improves the averaged evaluation
  result, although one of three seeds is marginally worse. More replicas are
  appropriate before a production decision.

The final independent evaluation event pass took 265.7 seconds and produced six
PNG/PDF pairs per region, a ROOT histogram archive, and a numerical JSON
summary under `plots/representative_test_r3_train80/validation-evaluation`.

## Output inventory

- Fit: `plots/representative_test_r3_train80/summary.json` and
  `high_fidelity_histograms.root`.
- Scan surfaces: abbreviated regional directories beneath that fit directory.
- Necessity analysis: `correction-need/correction_need.json`, CSV, and summary
  plots.
- Evaluation: `validation-evaluation`, with regional subdirectories, ROOT
  histograms, and `validation_summary.json`.
- Full-statistics diagnostic plots: `validation-all`, using the same fitted
  parameters and adaptive mass bins but all selected data and MC events.
- Smoke tests: `smoke-tests/core_smoke`, `smoke-tests/plot_style_test`, and the
  two-region batch-rendering regression in `smoke-tests/multiregion_plot_regression`.

## Limitations and next production step

The representative result uses an 11x11 grid, only three seeds, and only five
of the 17 retained eta-pair categories. It is deliberately a bounded test of
the new machinery. The scale step is 0.1% and the resolution step is 0.5%, so
grid values should not be quoted as final precision measurements.

The next production pass should use the 21x21 grid, increase replicas to at
least 10 after a timing check, and process all 17 retained regions. Repeated
cross-validation folds would use the finite MC more efficiently than relying
on a single 80/20 split. A final decision rule should be frozen before those
runs. The current pair-region model also continues to apply a common fitted
scale/resolution to both muons; a globally consistent per-muon eta calibration
would require a joint parameterization rather than 17 independent pair fits.
