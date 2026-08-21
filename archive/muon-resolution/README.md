# Muon resolution and scale study

First-pass reconstruction and shape comparison of the Z→μμ resonance in the
provided 2016H data and Powheg+Pythia8 simulation subsets.

Run the full local sample from this directory with:

```bash
conda run -n mg5 python analyze_zmumu.py
```

For a quick test:

```bash
conda run -n mg5 python analyze_zmumu.py \
  --max-files 1 --max-events 5000 --output-dir smoke-tests/quick
```

Use `python analyze_zmumu.py --help` for all options. The output contains PNG
and PDF comparisons, a ROOT histogram file, and a machine-readable JSON summary.
