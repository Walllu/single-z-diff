# Signed eta-pair scale and resolution scans

The 25 regions are ordered by the nominal leading and subleading muon eta bins.
Run one unsmeared point:

```bash
conda run -n mg5 python scan_zmumu_scale_resolution.py \
  --mode point --regions barrel__barrel --scale 0 --resolution 0
```

Run the configured grid with the chi-square objective:

```bash
conda run -n mg5 python scan_zmumu_scale_resolution.py \
  --mode scan --regions barrel__barrel --metric chi2
```

The current fine/wide grid contains 21 scale points from -0.5% to +0.5% in
0.05% steps and 21 added-resolution points from 0% to 5% in 0.25% steps.
Scan output directory names encode these bounds and point counts so earlier
coarse results are not overwritten.

Use `--regions all` for all 25 categories and `--metric poisson` for the binned
Poisson deviance. See `--help` and the editable configuration blocks near the
top of the script.

`--regions all` normally removes the categories listed in the top-level
`EXCLUDED_ETA_PAIR_ABBREVIATIONS` set. Use `--include-excluded-regions` to
restore them for diagnostic runs.

To create a browsable directory for every eta-pair beneath a directory labeled
with the point configuration:

```bash
conda run -n mg5 python scan_zmumu_scale_resolution.py \
  --mode point --regions all --scale 0 --resolution 0 \
  --region-subdirectories --output-dir plots
```

This produces, for example,
`plots/scale_p0p000_resolution_p0p000/PT_PT/`. The abbreviations are NE
(negative endcap), NT (negative transition), B (barrel), PT (positive
transition), and PE (positive endcap).

## Mass-binning utility

`optimize_mass_binning.py` recommends mass bin counts from selected data using
the Freedman–Diaconis rule, reports Scott and Sturges alternatives, and limits
the recommendation using a configurable minimum average population per bin.
For a three-file study:

```bash
conda run -n mg5 python optimize_mass_binning.py \
  --max-files 3 --minimum-events-per-bin 20 \
  --output-dir binning-tests/three-data-files
```

The recommendation is a starting point for fit-stability studies, not an
automatic replacement for `REGION_HISTOGRAM_OVERRIDES`.

## Applying the fitted correction

`muon_pair_correction.py` provides `Muon` and `correct_muon_pair`. It reads a
completed scan `summary.json`, selects the ordered nominal leading/subleading
eta-pair region, and applies independent Gaussian factors to the two simulated
muons. Supply a stable event key for event-to-event reproducibility:

```python
from muon_pair_correction import Muon, correct_muon_pair

result = correct_muon_pair(
    Muon(pt=45.0, eta=0.3, phi=1.2),
    Muon(pt=37.0, eta=-0.8, phi=-2.0),
    "plots/<scan>/summary.json",
    apply_correction=True,
    rng_seed=314159,
    event_key=f"{run}:{lumi}:{event}",
)
corrected_muons = result["muons"]
```

Excluded/out-of-acceptance pairs return identity factors. The correction is a
simulation smearing model and should not be applied to collision data.
