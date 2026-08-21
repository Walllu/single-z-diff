# First-Pass Z→μμ Muon Resolution and Scale Study

**Iteration:** 1 (baseline reconstruction and data/simulation comparison)  
**Date:** 2026-07-14  
**Analysis directory:** `single-z-diff/muon-resolution`

## Context and objective

This is the initial analysis pass toward a muon momentum scale and resolution
smearing correction. It establishes a reproducible Z→μμ selection, reconstructs
the dimuon resonance in the supplied 2016H data subset and Powheg+Pythia8
Z→μμ simulation subset, and produces normalized data/simulation comparisons.
It intentionally does **not** derive or apply a correction yet: the baseline
comparison must be validated before fitting scale shifts or extra MC smearing.

## Inputs and schema

- Data: `HackathonDataRaw/2016H` (6 ROOT files)
- Simulation: `HackathonDataRaw/ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8`
  (12 ROOT files found, 11 usable)
- Tree: `pfExtractor/pfTree`
- Runtime tested: the user's `mg5` conda environment, PyROOT/ROOT 6.40.02

The files use a custom `pfExtractor` schema. The script enables only required
branches, avoiding I/O for the much larger PF-candidate and jet collections.
The simulation file `14FF8F13-51E6-784E-BD14-B3C4417B3565.root` has no
`pfExtractor/pfTree`; it is detected and skipped with a warning.

## Baseline event and object selection

Events must pass `IsoMu24`. Reconstructed muons must satisfy:

- `pT > 15 GeV` and `|eta| < 2.4`;
- tight ID and PF-muon flags;
- `|dxy| < 0.05 cm` and `|dz| < 0.10 cm`;
- delta-beta corrected relative PF isolation below 0.15:
  `(charged + max(0, neutral + photon - 0.5*PU))/pT`.

All opposite-sign pairs are formed. If several exist, the pair whose invariant
mass is closest to 91.1876 GeV is retained. Its leading muon must have
`pT > 26 GeV`, and the final candidate must satisfy `60 < m(mumu) < 120 GeV`.
Four-vectors are constructed from the stored muon `pT`, `eta`, `phi`, and mass.

## Full-run results

| Cut | Data | Simulation |
|---|---:|---:|
| All processed | 865,402 | 728,000 |
| `IsoMu24` | 551,717 | 444,519 |
| At least two selected muons | 31,753 | 257,730 |
| Opposite-sign pair | 31,736 | 257,719 |
| Leading `pT > 26 GeV` | 31,046 | 253,967 |
| `60 < m(mumu) < 120 GeV` | **28,842** | **250,470** |

The simulation is normalized to the selected data integral independently in
each plot. This is a shape comparison, not an absolute-rate prediction: these
ntuples expose no event-weight/cross-section/luminosity bookkeeping sufficient
for a defensible luminosity normalization. The sharply different raw selection
efficiencies therefore should not be interpreted as a data/MC rate discrepancy.

The baseline mass spectrum shows a visible but modest localized data/simulation
shape difference around the Z peak, which is the intended handle for subsequent
scale and resolution extraction. No quantitative correction is claimed in this
iteration.

## Deliverables

- `analyze_zmumu.py`: command-line PyROOT event loop, selection, reconstruction,
  plotting, ROOT histogram output, and JSON cutflow output.
- `plots/`: full-sample PNG/PDF comparisons for dimuon mass, leading and
  subleading muon pT, muon eta, Z pT, and vertex multiplicity; also
  `histograms.root` and `summary.json`.
- `smoke-tests/first-pass/`: 5,000-event-per-sample smoke test (one file each),
  with 160 selected data and 1,664 selected simulation candidates.
- `smoke-tests/two-files/`: stronger all-event, two-input-file validation, with
  8,890 selected data and 17,851 selected simulation candidates.

Run the full analysis from `muon-resolution` with:

```bash
conda run -n mg5 python analyze_zmumu.py
```

## Validation performed

- Python syntax compilation passed.
- Both smoke tests completed without empty histograms.
- The complete usable input set completed and wrote all expected artifacts.
- Output histograms carry statistical uncertainties (`Sumw2`).
- Each comparison includes a data/MC ratio panel and MC is shape-normalized.
- The dimuon mass PNG was visually inspected for axes, legend, peak, and ratio.

## Limitations and next iteration

This first pass does not apply pileup, trigger, ID, isolation, acceptance, or
generator weights. The vertex comparison can help diagnose pileup mismodelling.
The next correction-focused iteration should split the sample into muon eta and
charge bins and fit the Z line shape (or use a template likelihood) to extract
momentum-scale shifts and the additional resolution required in simulation.
That work should include fit-model/systematic variations, deterministic random
seeds for smearing, closure tests, and correction parameter serialization.
