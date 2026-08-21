# Regional reconstructed-vertex pileup reweighting

This package derives approximate pileup weights from the selected `nVertices`
distributions in the unordered BB, BE, and EE event categories. It then applies
the weights to the existing best muon-scale/resolution-corrected simulation and
validates their effects on muon, event, and Z-candidate observables.

This is reconstructed-vertex reweighting, not official true-pileup reweighting.
The input ntuples do not expose a true number of simulated interactions or the
luminosity-derived data pileup profile.

## Contents

- `pileup_weights.py`: ROOT-independent downstream lookup API.
- `derive_and_validate_pileup.py`: weight derivation, weighted histogramming,
  numerical evaluation, and isolated ROOT rendering.
- `../../../calibrations/experimental/pileup_nvertices_2016H_bb_be_ee.json`:
  tracked exploratory weight table.
- `plots`: regenerated numerical effects, templates, and figures; not
  versioned.
- `reports/RECONSTRUCTED_VERTEX_PILEUP_REWEIGHTING.md`: method and results.
- `smoke-tests`: bounded one-file validation products.

## Derive and validate

Run from the `mg5` environment:

```bash
conda run -n mg5 python derive_and_validate_pileup.py \
  ../plots/full_production_r3/summary.json
```

The production defaults require at least 100 selected data events and 300
selected MC events in each contiguous vertex group. Sparse final bins are
merged into the preceding group. Values at or above 50 vertices use the final
overflow weight.

To regenerate a plot without rereading the ntuples:

```bash
conda run -n mg5 python derive_and_validate_pileup.py \
  ../plots/full_production_r3/summary.json --render-only \
  --regions BE --observables z_pt
```

## Apply downstream

```python
from pileup_weights import PileupCalibration

pileup = PileupCalibration.from_json(
    "../../../calibrations/experimental/pileup_nvertices_2016H_bb_be_ee.json"
)

event_weight = pileup.weight("BE", nvertices, apply=True)
```

Multiply this event weight by other MC weights. Do not apply it to data. The
`region` must be the same unordered BB/BE/EE definition used during derivation.
