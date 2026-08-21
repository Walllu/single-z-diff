# Energy Flow Gap Observables in the Simplified Z→μμ Analysis

## Context and lineage

This is an extension of the week-scale diffractive Z-to-dimuon analysis and of
the preceding simplified BB/BE/EE muon scale/resolution work. The immediate
goal is to reproduce a colleague's Energy Gap Flow observable in the local
`pfExtractor/pfTree` samples, compare data with inclusive DY simulation, and
provide first diagnostics for separating the locally available diffractive
`pomflux` simulation.

An additional, deliberately separate reproduction of the colleague notebook's
own loose-muon selection and CMS-style plots is documented in
`COLLEAGUE_NOTEBOOK_REPRODUCTION.md`. That production applies no muon scale or
resolution correction and should be used for direct comparison with the
notebook; this report continues to describe the tighter calibration-selection
study.

Two colleague inputs are retained beside the implementation:

1. `gistfile1.txt` supplies the numerical EFG definition.
2. `Analysis.ipynb` supplies the operational candidate cleaning and the idea
   of signal/background shape and ROC comparisons.

The notebook and this study do not use identical Z selections. The notebook
uses two loose muons above 22 GeV, 80–100 GeV dimuon mass, and Z pT below
15 GeV. This implementation deliberately retains the established calibration
selection: `IsoMu24`, tight PF and isolated muons, pT above 15 GeV, leading pT
above 26 GeV, opposite charge, and 70–110 GeV mass. This keeps the new plots
directly comparable with the existing BB/BE/EE calibration products.

## Implemented definition

For each selected event, candidates must satisfy:

- `candPt > 0.5 GeV`;
- `candFromPV == 3`;
- `-2.5 < candEta < 2.5`;
- no overlap within ΔR < 0.01 of either selected Z muon.

The last three numerical choices reproduce the notebook workflow. The EFG
parameters are β = 0.1 and Q0 = 2 GeV. Scores are calculated independently
from the negative and positive detector sides, and the event discriminator is
their maximum. The code additionally stores the effective EFG gap, the hard
forward-gap analogue, both one-sided scores, cleaned charged multiplicity and
cleaned scalar pT sum. An inclusive EFG variant, before Z-muon overlap
removal, makes the cleaning effect explicit.

An event with no accepted candidate follows the gist convention: EFG score 1
and full effective/hard gap. This differs from the notebook's provisional
conventional-gap fallback of zero and is recorded in every JSON summary.

## Data/MC and signal diagnostics

The main plots compare shape-normalized data to inclusive DY simulation before
and after the previously derived BB/BE/EE reconstructed-vertex weights. The
lower panels show data/MC and Beeston–Barlow profile pulls with ±1/2/3σ guides.

The locally available diffractive `DYToMuMu_pomflux` sample is processed with
the same selection and candidate definition. Separate normalized signal versus
inclusive-DY shapes and ROC curves follow the notebook's exploratory design.
High values are treated as signal-like for gap observables; low values are
treated as signal-like for candidate multiplicity and scalar pT sum.

## Reproducibility and outputs

The standalone algorithm is in `energy_flow_gap.py`; the directed event loop
is `analyze_energy_flow_gap.py`. The latter writes all templates to a ROOT file
and all configuration, cutflows, event counts, fit metrics, and timing to JSON.
The one-file smoke test is under
`smoke-tests/notebook_candidate_definition`. The full production is under
`plots/full_production_notebook_definition`, with one subdirectory per BB/BE/EE
region.

## Full-production results

The completed run processed 4,389,502 data entries from 31 files, 1,914,000
inclusive-DY entries from 37 usable files, and 200,000 diffractive entries from
6 files. One inclusive-DY ROOT file lacked `pfExtractor/pfTree` and was skipped
by the existing input validation. The extraction and histogramming stage took
546.6 seconds; the complete run including plot rendering took 576.1 seconds on
the local machine.

After the established Z selection the regional sample sizes were:

| region | data | inclusive DY | diffractive simulation |
|---|---:|---:|---:|
| BB | 68,008 | 306,304 | 47,846 |
| BE | 57,396 | 259,701 | 34,710 |
| EE | 18,436 | 81,091 | 13,217 |

There were 445 data, 2,241 inclusive-DY, and 5,263 diffractive selected events
without a remaining accepted candidate. These are a small but non-negligible
population and are represented using the full-gap convention described above.

The regional reconstructed-vertex weights improve the EFG-score
Beeston–Barlow objective in every category: 1959.5→1769.6 in BB,
1909.5→1738.3 in BE, and 774.1→707.4 in EE (49 nominal degrees of freedom).
The objectives remain much larger than the nominal degrees of freedom, so the
weighting helps but plainly does not make the inclusive-DY candidate activity
an adequate description of data. This is an important systematic issue for a
future EFG selection, rather than evidence to tune the observable itself.

As an exploratory simulation-only discriminator, the cleaned maximum EFG
score gives histogram ROC AUC values of approximately 0.884 (BB), 0.880 (BE),
and 0.876 (EE) for diffractive signal against inclusive DY. This strong shape
separation is encouraging, but must be considered together with the observed
data/MC mismodelling.

The bounded smoke test processed one file per sample (30,000 entries per main
sample and the available 10,000 signal entries) and generated 47 files. The
full production generated 218 files (4.8 MB), including ROOT/JSON products and
PNG/PDF/SVG plots.

## Interpretation cautions

- The reconstructed-vertex weights are a pragmatic regional correction, not a
  truth-pileup prescription. Both weighted and unweighted DY are therefore
  shown.
- EFG is sensitive to PF-candidate reconstruction, tracking efficiency,
  primary-vertex association, and pileup. Data/MC disagreement should not be
  interpreted immediately as diffractive physics.
- The diffractive sample is useful for shape intuition and cut design, but its
  normalization is not used here.
- ROC curves are simulation-to-simulation diagnostics and do not establish a
  calibrated signal efficiency in data.
