# Reproduction of the Colleague Energy Gap Flow Notebook

## Scope

This production reproduces the analysis logic in `Analysis.ipynb` using a new
PyROOT event loop and newly written CMS-style plotting code. It is deliberately
separate from the BB/BE/EE muon-calibration workflow. No muon momentum scaling,
resolution smearing, candidate smearing, pileup reweighting, generator weight,
or luminosity normalization is applied.

A second set of productions now applies the fitted muon scale and resolution
to both MC samples. The original uncorrected outputs remain the direct notebook
reproduction; the corrected outputs test how portable that calibration is to
the notebook selection.

Two productions are provided:

- `notebook_file_counts_no_muon_corrections` uses the notebook's sample sizes:
  2 data files, 3 inclusive-DY files and 3 diffractive files. Files are sorted
  deterministically, whereas the notebook used the first entries returned by
  `glob`, so exact event identity may differ.
- `full_samples_no_muon_corrections` uses every locally usable file and is the
  preferred production for checking data/MC shapes.
- Corresponding `*_with_muon_corrections` directories use the same inputs and
  analysis definition but apply the fitted correction to simulated muons.

## Samples and their roles

| short name | local sample | role in this study |
|---|---|---|
| Data | `2016H` | Recorded CMS Open Data; shown as black points |
| Inclusive DY | `ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8` | Ordinary, non-diffractive Z/DY background simulation; blue solid line |
| Single-diffractive | `DYToMuMu_pomflux_Pt-30_TuneCP5_13TeV-pythia8` | Exploratory diffractive signal simulation; red dashed line |

The signal MC is overlaid for shape intuition. It is not expected to describe
the complete data distribution. The data contain overwhelmingly ordinary DY
plus other backgrounds, with at most a small diffractive component.

## Reproduced event and candidate selections

The event must pass `IsoMu24` and contain at least two muons satisfying
`muonPt > 22 GeV` and `muonIsLoose == 1`. As in the notebook, the first two
accepted muons in the stored collection are used. The audit found no event in
the notebook-sized run where this pair was not already pT ordered. The dimuon
mass is calculated in the massless-muon approximation and must satisfy
80 < mμμ < 100 GeV. The vector-sum dimuon pT must be below 15 GeV. The notebook
does not require opposite charge, so this reproduction does not add that cut.

In corrected productions, the MC muon momenta are corrected before the 22 GeV
threshold, dimuon-mass window and dimuon-pT cut. This allows events to migrate
across those boundaries. Data momenta are unchanged. Three deterministic seed
replicas are averaged, matching the seeds used in the calibration fit.

Candidates must have `candPt > 0.5 GeV` and `candFromPV == 3`. Every candidate
within ΔR < 0.01 of either selected muon is removed. This is intended to leave
primary-vertex charged activity other than the Z decay muons.

## What each distribution plot contains

Every distribution plot contains all three samples and a data/inclusive-DY
ratio panel. Data retain Poisson √N error bars. Each MC histogram is scaled
independently to the data integral inside the plotted range. Consequently these
are shape comparisons, not cross-section predictions; the vertical label says
that MC is normalized to data. The signal normalization carries no statement
about the diffractive rate.

### Cleaned candidate pT

This histogram contains one entry per accepted candidate, not one per event,
in 29 equal bins from 0 to 2 GeV. Candidates above 2 GeV are outside the plotted
range, matching the notebook. Because multiple candidates come from the same
event, the displayed √N uncertainties do not include within-event correlations.

### Cleaned candidate η

This contains one entry per accepted candidate in 29 equal bins from −4 to 4.
The candidate cleaning itself does not impose an η cut. The EFG and rapidity-gap
calculations below subsequently use only |η| < 2.5.

### Maximum forward rapidity gap

For the cleaned candidates with pT > 0.5 GeV and |η| < 2.5, the observable is

`max(2.5 − max(η), 2.5 + min(η))`.

It is the larger empty interval between the outermost accepted candidate and
the positive or negative tracking edge. Higher values are more gap-like. To
reproduce the notebook, an event without any accepted candidate is assigned
zero here, even though a full detector gap would be a physically reasonable
alternative convention.

### Energy Gap Flow score

The EFG algorithm is evaluated separately after orienting the event toward the
negative and positive detector sides. Rather than depending only on the
outermost candidate, it continuously suppresses a putative gap according to
the cumulative transverse momentum beyond each candidate. The parameters are
β = 0.1, Q0 = 2 GeV, pT > 0.5 GeV and −2.5 < η < 2.5. The plotted event score is
`max(G−, G+)`. It lies between zero and one; values near one are strongly
gap-like. An event with no accepted candidate has score one by construction.

### Threshold-efficiency plots

For both event-level scores, these plots show the fraction of data, inclusive
DY and diffractive simulation remaining above every possible threshold. Unlike
the ROC, they can and do include data. They directly expose whether a proposed
gap cut has a similar efficiency in data and inclusive-DY simulation.

## ROC curve and AUC

The ROC curve necessarily uses only the two MC samples because it needs known
class labels: diffractive simulation is labelled signal and inclusive DY is
labelled background. Data have no per-event truth label and cannot validly be
inserted into this curve.

To mirror the colleague's presentation, the horizontal coordinate is the
single-diffractive signal efficiency and the logarithmic vertical coordinate is
the inclusive-DY background efficiency. A useful curve therefore bends toward
the lower-right corner. The AUC quoted in the legend is the conventional
`∫ signal efficiency d(background efficiency)` obtained by scanning from high
to low score. AUC = 0.5 corresponds to random ordering and AUC = 1 to perfect
ordering.

In the notebook-sized reproduction:

| observable | AUC | DY efficiency at 50% signal | at 80% signal | at 90% signal |
|---|---:|---:|---:|---:|
| Maximum rapidity gap | 0.770 | 0.089 | 0.396 | 0.822 |
| Energy Gap Flow | 0.895 | 0.059 | 0.154 | 0.246 |

Thus EFG provides materially better global ranking. At 80% simulated signal
efficiency, for example, it retains about 15% of inclusive DY rather than about
40% for the hard gap. These are simulation-only expectations, not validated
data efficiencies.

## Notebook-sized reproduction result

The deterministic local subset processed 278,942 data entries, 150,000
inclusive-DY entries and 62,000 diffractive entries. After selection it retained
5,466 data, 29,631 inclusive-DY and 18,962 diffractive events. Its complete
configuration, cutflows, candidate counts, normalization factors, AUCs and
working points are stored in `colleague_reproduction_summary.json` beside the
plots.

The saved notebook itself reports 5,966 selected data, 22,133 inclusive-DY and
18,630 diffractive events. Those counts demonstrate that its unsorted EOS glob
selected different physical files from the deterministic local subset. Without
an explicit filename manifest in the notebook, exact event-for-event replay is
not defined. The local reproduction therefore claims equivalence of the code
path, selection, file counts and observable definitions—not identity of the
input events.

## Full-sample result

The full production uses 31 data files, 37 usable inclusive-DY files and all 6
diffractive files. One additional inclusive-DY file does not contain
`pfExtractor/pfTree` and is skipped by input validation. It processes 4,389,502
data, 1,914,000 inclusive-DY and 200,000 diffractive entries, retaining 84,447,
382,104 and 61,830 events respectively after the notebook selection.

The full-sample ROC results are:

| observable | AUC | DY efficiency at 50% signal | at 80% signal | at 90% signal |
|---|---:|---:|---:|---:|
| Maximum rapidity gap | 0.770 | 0.089 | 0.394 | 0.818 |
| Energy Gap Flow | 0.895 | 0.059 | 0.157 | 0.249 |

Their near identity to the deterministic notebook-sized result shows that the
relative EFG performance is not driven by the larger sample or by a fortunate
three-file choice. The full extraction and plot production takes approximately
198 seconds on the local machine.

The full data comparisons also make an important qualification visible. Data
generally contain more gap-like activity than the inclusive-DY simulation in
parts of the score range. The EFG threshold-efficiency curve for data lies
between inclusive DY and the diffractive signal over much of the range. That is
qualitatively compatible with either a diffractive contribution or imperfect
modelling of pileup/charged activity, but these shape plots alone cannot
separate those explanations.

## Plotting implementation

The notebook mixed an ATLAS Matplotlib style with a CMS label and normalized
the two MC histograms through a compact helper. The new implementation does not
reuse that plotting code. It builds the figures independently with ROOT and a
consistent CMS/Open Data visual design: black data markers, blue inclusive-DY,
red dashed signal, restrained typography, logarithmic event-score panels, and
an explicit data/DY ratio. ROOT was also the reliable choice in the requested
`mg5` environment because its installed Matplotlib requires a newer NumPy than
the environment currently provides. The final plotting workflow writes SVG
vector files only; PNG and PDF copies are intentionally not produced.

## Interpretation cautions

The data/inclusive-DY differences are essential. Candidate activity is
sensitive to pileup, vertex association, tracking and the underlying-event
model. A good signal/background AUC does not compensate for a badly modelled
data distribution. Before turning EFG into a final selection, the data/DY ratio
and threshold-efficiency curves should be used to define a tolerable operating
region and associated modelling uncertainty.

## Does the muon correction apply in this phase space?

Partly. The correction is a per-muon function of whether |η| is below or above
1.4, so it is not intrinsically tied to an event being labelled BB, BE or EE.
It can be evaluated for any simulated muon with |η| < 2.5. However, the fit was
derived with tight, isolated PF muons, pT > 15 GeV, a leading-muon threshold of
26 GeV and 70–110 GeV dimuon mass. The notebook instead uses loose muons above
22 GeV, 80–100 GeV mass, pT(Z) < 15 GeV and no explicit opposite-sign or muon-η
requirement. Applying the fitted values here is therefore a portability test
with overlapping but non-identical phase space, not a new calibration fit.

The notebook's lack of a muon-η cut also admits a small number of muons outside
the calibrated |η| < 2.5 domain. The code does not extrapolate: such muons are
left uncorrected and counted in each JSON summary.

## Corrected-production result

The applied calibration is a −0.125% barrel scale shift with 0.375% added
Gaussian resolution and a −0.150% endcap scale shift with 0.875% added
resolution. Both simulated samples use the three deterministic calibration-fit
seeds 314159, 271828 and 161803. Each event contributes the fraction of replicas
that passes the corrected selection. Data are unchanged.

On the full samples, the effective selected inclusive-DY yield changes from
382,104 to 382,121 and the diffractive yield from 61,830 to 61,840. The EFG AUC
changes from 0.8948737 to 0.8948751, while the conventional-gap AUC changes from
0.7696533 to 0.7696944. The corrected operating points are:

| observable | AUC | DY efficiency at 50% signal | at 80% signal | at 90% signal |
|---|---:|---:|---:|---:|
| Maximum rapidity gap | 0.770 | 0.089 | 0.394 | 0.818 |
| Energy Gap Flow | 0.895 | 0.059 | 0.157 | 0.249 |

The normalized-shape total-variation distance between corrected and nominal
inclusive-DY EFG is 0.00017 (0.017%); for diffractive MC it is 0.00051 (0.051%).
The candidate-pT, candidate-η and hard-gap differences are similarly below
0.05%. Thus the corrected plots are practically identical to the uncorrected
ones. This is not a failure of the calibration: EFG is calculated from cleaned
non-muon candidates, and the muon correction can affect it only by migrating a
small number of events through the selection boundaries.

Across the full preselection traversal, 6,156 inclusive-DY and 735 diffractive
loose-muon replica instances lie outside |η| < 2.5 and are deliberately left
uncorrected. These are replica-level diagnostic counts, not selected-event
counts.

The corrected SVGs and machine-readable products are in
`notebook_file_counts_with_muon_corrections` and
`full_samples_with_muon_corrections`. The full corrected production takes about
490 seconds because each loose simulated muon receives three independently
keyed Gaussian evaluations.
