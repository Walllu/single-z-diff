# Energy Flow Gap study

This directory implements the Energy Gap Flow (EFG) observable from the
colleague material under `../../../../../reference/energy-flow-gap/` in the
existing simplified BB/BE/EE Z-to-dimuon workflow.

## Files

- `energy_flow_gap.py`: reusable numerical implementation of the EFG score,
  effective gap, conventional forward gap, and an event-level wrapper.
- `analyze_energy_flow_gap.py`: ROOT event loop, regional pileup weighting,
  histogram production, ratio/pull plotting, diffractive-signal shape checks,
  and ROC plotting.
- `reproduce_colleague_efg.py`: independent reproduction of the notebook's
  loose-muon selection, candidate plots, event scores, ROC/AUC, and CMS-style
  data/MC comparisons, with no muon corrections.
- `../../../../../reference/energy-flow-gap/`: colleague notebook and source
  definition used as attributed workflow references.
- `smoke-tests/` and `plots/`: regenerated outputs; not versioned.
- `reports/`: implementation and result notes.

## Default definition

The Z selection is intentionally the same tight, opposite-sign selection used
by the simplified muon-calibration analysis. Candidate cleaning follows the
notebook: `candPt > 0.5 GeV`, `candFromPV == 3`,
`-2.5 < candEta < 2.5`, and removal of every candidate within
`deltaR < 0.01` of either selected Z muon. The primary event score is the
larger of the positive- and negative-side EFG scores.

## Running

From the repository root:

```bash
conda run -n mg5 python \
  single-z-diff/feature/muon-resolution/scans/simplified-scans/EFG/analyze_energy_flow_gap.py \
  --label full_production_notebook_definition
```

Useful controls include `--max-files`, `--max-events`, `--regions`,
`--observables`, `--skip-signal`, and `--render-only`. Each production writes
`efg_histograms.root`, `efg_summary.json`, and BB/BE/EE plot subdirectories.

To reproduce the notebook file counts without any muon correction:

```bash
conda run -n mg5 python \
  single-z-diff/feature/muon-resolution/scans/simplified-scans/EFG/reproduce_colleague_efg.py \
  --notebook-file-subset \
  --label notebook_file_counts_no_muon_corrections
```

Omit `--notebook-file-subset` to use every locally available file. Reproduction
plots and their JSON/ROOT products are stored under `plots/colleague-reproduction`.
Plots from this workflow are written as SVG only.

To apply the fitted BB/BE/EE muon correction to both MC samples, add
`--apply-muon-corrections`. The default calibration is
`../../../calibrations/muon_momentum_2016H_bb_be_ee.json`, and the default
replicas use the three fit seeds. Data are never corrected.
