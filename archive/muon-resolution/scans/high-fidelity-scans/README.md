# High-fidelity muon scale/resolution scans

This directory is an independent successor to the original signed eta-pair
scan. It adds event-keyed random numbers, multi-seed template averaging,
resonance-aware adaptive mass bins, analytic single-source Barlow--Beeston MC
statistical profiling, deterministic training/evaluation splits, correction-
necessity assessment, and validation plots with seed envelopes and pulls.

The default replica count is deliberately small (`3`). The available grids are:

- `smoke`: 5x5;
- `test`: 11x11;
- `production`: 21x21, matching the preceding fine scan ranges.

Example representative scan:

```bash
conda run -n mg5 python high_fidelity_scan.py \
  --grid-profile test --replicas 3 \
  --regions neg_endcap__neg_endcap neg_transition__barrel \
            pos_transition__pos_endcap pos_endcap__pos_transition \
            barrel__barrel \
  --label representative_test_r3_train80
```

Assess whether the fitted corrections are supported over the identity point:

```bash
conda run -n mg5 python assess_correction_need.py \
  plots/representative_test_r3_train80/summary.json
```

Make independent-evaluation plots for all six observables:

```bash
conda run -n mg5 python plot_high_fidelity_validation.py \
  plots/representative_test_r3_train80/summary.json
```

Make higher-statistics diagnostic plots using the complete selected samples:

```bash
conda run -n mg5 python plot_high_fidelity_validation.py \
  plots/representative_test_r3_train80/summary.json --subset all
```

The `validation-all` plots include the training events and therefore should not
be interpreted as an independent closure test.

Validation figures are written in PNG, SVG, and PDF. SVG/PDF are recommended
for inspection on macOS if the PyROOT bitmap preview drops glyphs.

Restyle plots without rereading ROOT event trees:

```bash
conda run -n mg5 python plot_high_fidelity_validation.py \
  plots/representative_test_r3_train80/summary.json --render-only
```

For stable macOS ROOT batch output, multi-region rendering is automatically
isolated into one short-lived subprocess per region.

Use `--seeds` to provide an explicit seed list, or `--replicas N` to take the
first `N` entries from the seed table. Input directories, file/event limits,
split configuration, bin population thresholds, requested regions, and output
locations can all be directed from the command line.

The complete implementation and representative-test findings are documented
in `reports/HIGH_FIDELITY_SCALE_RESOLUTION_SCAN.md`.
