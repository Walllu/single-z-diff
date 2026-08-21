# Muon-correction validation

This directory compares data with both uncorrected simulation and corrected
simulation for multiple reproducible RNG seeds. It uses the fixed nominal event
selection and eta-pair categories from the scan, and applies the fitted
simulation-only correction through `muon_pair_correction.correct_muon_pair`.

Full run:

```bash
conda run -n mg5 python validate_muon_corrections.py
```

Quick test:

```bash
conda run -n mg5 python validate_muon_corrections.py \
  --regions barrel__barrel --seeds 314159 271828 \
  --max-files 1 --max-events 5000 \
  --output-dir smoke-tests
```

Each MC histogram is normalized independently to the data integral. Ratio
panels show data/uncorrected-MC and data/corrected-MC for every requested seed.

The completed four-seed run is in `plots/tiered_v1_seeds_4`. It contains six
observables in PNG and PDF form for each of the 17 calibrated eta-pair regions,
plus `validation_histograms.root` and `validation_summary.json`. See
`reports/CORRECTION_VALIDATION.md` for the configuration, results, and caveats.
