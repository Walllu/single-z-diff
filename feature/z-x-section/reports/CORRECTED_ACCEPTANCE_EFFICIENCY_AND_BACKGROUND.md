# Corrected Dressed-Z Acceptance, Efficiency, and Background Prototype

## Scope

This iteration makes the local inclusive-Z validation structurally complete
while deliberately leaving unavailable external inputs unset. It implements
both fiducial and acceptance-extrapolated interpretations, generator truth,
reconstruction matching, muon corrections and variations, an optional toy
pileup correction, and a provisional smooth sideband model.

No absolute cross section is quoted. The local file subset has no certified
luminosity, and the MC sample lacks the external cross section and production
normalization metadata needed for prompt-background subtraction.

## Physics definition

The configurable default is:

- stable generator muons descended from a PDG-23 particle;
- stable PDG-23-descended photons added to the nearest selected muon within
  `deltaR < 0.1`;
- opposite-sign dressed pair closest to the Z pole;
- `80 <= m(dressed mumu) < 100 GeV`;
- both dressed muons with `pT > 25 GeV` and `|eta| < 2.4` for the fiducial
  numerator;
- identical reconstructed kinematic bounds, plus quality, trigger matching,
  opposite sign and isolation requirements;
- reconstructed-to-dressed matching within `deltaR < 0.1` with compatible
  charge assignments.

The 80-100 GeV truth window was chosen instead of 90-100 GeV because the latter
is asymmetric about the 91.19 GeV pole and removes much of the physical low-side
FSR tail. A 50-120 GeV alternative is recorded in the configuration but should
only be adopted together with a validated background treatment.

PDG-23 ancestry provides the requested operational “Z only” MC definition.
It does not make Z, virtual-photon, and interference contributions separately
observable in data; a publication-level result should normally state a
mass-window `Z/gamma*` definition.

## Corrections and weights

The MC momentum correction is applied before muon ordering, pT selection,
relative-isolation calculation and reconstructed mass. Stable event-keyed
Gaussian draws make the smearing independent of processing order.

The nominal selected-event efficiency weight is

```text
genWeight * normalization * pileup
* Reco+ID(mu1) * Reco+ID(mu2)
* Iso(mu1) * Iso(mu2)
* Trigger(highest-pT matched muon).
```

The collaborator reconstruction map is provisionally interpreted as including
tight ID. Its final column is repeated through the missing 75-100 GeV bins.
All supplied reconstruction, isolation, and trigger statistical/systematic
maps are implemented as coherent up/down variations. This simple correlation
model and the matched-muon trigger prescription require confirmation.

The experimental BB/BE/EE reconstructed-vertex pileup weights are off by
default. They can be enabled as an impact study, but are not an official
true-pileup correction.

## Full local result

Every locally readable file was used: 31 data files and 31 MC files with
`genWeight`. Six MC files without `genWeight` and one file without the event
tree remain excluded.

| Configuration | Acceptance | Efficiency | A x efficiency |
|---|---:|---:|---:|
| No momentum or efficiency correction | 0.38453 | 0.83087 | 0.31949 |
| Momentum only | 0.38453 | 0.82980 | 0.31908 |
| Efficiency SFs only | 0.38453 | 0.78866 | 0.30326 |
| Nominal momentum + SFs | 0.38453 | 0.78747 | 0.30280 |
| Nominal + experimental vertex PU | 0.38453 | 0.78321 | 0.30117 |

Acceptance is unchanged because it is defined entirely at dressed truth level.
Momentum corrections affect reconstructed migrations at the per-mille level in
this setup. The preliminary efficiency SFs lower the total efficiency by about
5.2%. The toy pileup correction produces a further 0.54% relative reduction.

The coherent SF variations shift the selected efficiency by approximately:

| Source | Up | Down |
|---|---:|---:|
| Reco+ID systematic | +3.85% | -3.76% |
| Reco+ID statistical | +5.55% | -5.40% |
| Isolation systematic | +0.86% | -0.86% |
| Isolation statistical | +3.43% | -3.37% |
| Trigger systematic | +0.91% | -0.91% |
| Trigger statistical | +3.68% | -3.68% |

These large preliminary statistical variations agree with the provider's note
that the maps need a larger derivation sample.

## Background prototype

The existing nominal two-anti-isolated ABCD count in 80-100 GeV gives 29.25 raw
events but has no prompt subtraction and no validated closure. The local
constant closure-factor study also requires a factor much larger than one.

An additional positive exponential residual was fitted to the corrected MC and
data selections in 60-75 and 105-120 GeV after provisionally shape-normalizing
DY in the signal peak. It interpolates 425.5 events into 80-100 GeV. This value
is **not accepted as a background estimate**: fitting only the low sideband
predicts the high sideband with a deviance of about 7,693, while the high-only
fit predicts 43.3 signal-window events and the low-only fit predicts 2,333.9.
The failed holdout demonstrates that a single smooth residual is not stable
with the current prompt normalization.

The model implementation is retained because it becomes meaningful once DY is
externally normalized. A simple background model is reasonable for a fast
80-100 GeV measurement only after low/high sideband compatibility and ABCD
closure are demonstrated.

## Remaining external inputs

- certified luminosity and lumimask;
- DY cross section, filter efficiency, k-factor decision, and total sum of
  generator weights;
- replacement or metadata recovery for six weightless MC files;
- official pileup weights and variations;
- confirmation of the SF binning/correlations and event-trigger prescription;
- a decision about whether additional prompt backgrounds are negligible or
  require explicit samples and normalization uncertainties.

Once those values are placed in `measurement_config.json`, the code can report

```text
sigma_fid  = (Ndata - Nbkg) / (efficiency * luminosity)
sigma_full = (Ndata - Nbkg) / (acceptance * efficiency * luminosity).
```

The current code intentionally leaves both values null rather than filling them
with an arbitrary luminosity.
