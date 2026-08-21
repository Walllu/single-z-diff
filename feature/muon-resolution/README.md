# Muon calibration and diffractive-Z features

The active calibration uses a simultaneous BB/BE/EE grid fit. Its compact payload is `calibrations/muon_momentum_2016H_bb_be_ee.json`; the reusable downstream API is `scans/simplified-scans/muon_corrections.py`.

`scans/scan_zmumu_scale_resolution.py` and `scans/high-fidelity-scans/high_fidelity_scan.py` are retained because the simplified scan imports their event selection, dynamic binning, and Barlow--Beeston implementation.

The reconstructed-vertex pileup payload is stored under
`calibrations/experimental/` because it is an exploratory correction rather than a final true-pileup calibration.

Input ROOT files are not versioned. The scripts locate `HackathonDataRaw` in a parent workspace directory or accept explicit `--data-dir` and `--mc-dir` arguments. Generated products default to local `plots/` directories, which should remain ignored by Git.
