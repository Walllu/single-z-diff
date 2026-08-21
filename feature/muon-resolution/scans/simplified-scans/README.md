# Simplified BB/BE/EE muon calibration

This directory contains a simultaneous, shape-only dimuon-mass calibration of
shared per-muon barrel and endcap scale/resolution parameters. Event categories
are unordered: `BE` means exactly one barrel and one endcap muon, irrespective
of which muon is leading in transverse momentum.

## Contents

- `simplified_scan.py`: reads the ROOT samples, applies the established muon and
  event selection, constructs BB/BE/EE templates, performs the joint grid fit,
  and writes fit landscapes plus mass ratio/pull plots.
- `muon_corrections.py`: ROOT-independent downstream correction API.
- `plot_simplified_validation.py`: regenerates plots from a saved summary and
  ROOT histogram file without rereading the ntuples or refitting.
- `plot_simplified_observables.py`: builds the full muon, event, and Z-kinematic
  validation suite from the ntuples, then supports render-only regeneration
  from its saved ROOT histograms.
- `pileup-reweighting`: regional reconstructed-vertex weight derivation,
  reusable event-weight API, and before/after kinematic validation plots.
- `../../calibrations/muon_momentum_2016H_bb_be_ee.json`: tracked compact
  calibration payload from the full 31-file data and 37-file usable-MC fit.
- `plots` and `smoke-tests`: regenerated output directories; not versioned.
- `reports/SIMPLIFIED_BB_BE_EE_SCALE_RESOLUTION_SCAN.md`: implementation and
  result report, including limitations and lineage.

## Run the fit

```bash
conda run -n mg5 python simplified_scan.py \
  --grid-profile production --replicas 3 \
  --output-dir plots --label full_production_r3
```

The defaults locate the workspace-level `HackathonDataRaw/2016H` and
`HackathonDataRaw/ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8`
directories. Explicit input paths can also be supplied on the command line.

## Regenerate plots

```bash
conda run -n mg5 python plot_simplified_validation.py \
  plots/full_production_r3/summary.json
```

Build the complete full-sample observable suite:

```bash
conda run -n mg5 python plot_simplified_observables.py \
  plots/full_production_r3/summary.json
```

The suite includes dimuon mass; leading/subleading muon pT, eta, and phi;
reconstructed vertex multiplicity; and Z-candidate pT, rapidity, and phi. Every
plot contains data, nominal MC, corrected MC, both ratios, signed profile pulls,
and labeled +/-1, +/-2, and +/-3 sigma pull guides. Numerical metrics are stored
in `observable_validation_summary.json`, and templates in
`observable_histograms.root`.

To regenerate one plot without rereading the ROOT inputs:

```bash
conda run -n mg5 python plot_simplified_observables.py \
  plots/full_production_r3/summary.json --render-only \
  --regions BE --observables z_pt
```

Individual regions or landscapes can be isolated if Cocoa ROOT develops a
canvas backing-store problem:

```bash
conda run -n mg5 python plot_simplified_validation.py SUMMARY --regions BE
conda run -n mg5 python plot_simplified_validation.py SUMMARY \
  --landscape-only fine_resolution
```

## Apply the calibration downstream

```python
from muon_corrections import Calibration, correct_dimuon

calibration = Calibration.from_summary(
    "../../calibrations/muon_momentum_2016H_bb_be_ee.json"
)

mu1_corr, mu2_corr, corrected_mass = correct_dimuon(
    (mu1_pt, mu1_eta, mu1_phi, mu1_mass),
    (mu2_pt, mu2_eta, mu2_phi, mu2_mass),
    calibration,
    seed=314159,
    event_key=f"{run}:{lumi}:{event}",
    apply=True,
)
```

Pass a unique, stable event key. The two muons receive independent hash-keyed
Gaussian draws, and rerunning with the same seed and keys reproduces exactly the
same correction. Use `apply=False` for an identity correction. Apply this to
simulation only; it is not a data correction.
