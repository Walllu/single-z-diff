# Toy inclusive Z cross-section study

This directory contains a local-first reconstruction and background-control
workflow for a fast inclusive `Z -> mu+ mu-` cross-section exercise. It now
contains reconstructed selection and background controls plus a configurable
dressed-generator-muon acceptance/efficiency pass. The distributed workflow
now produces additive inputs for a separate finalization stage. Absolute cross
sections remain disabled until approved luminosity and sample-normalization
metadata are provided.

## Why PyROOT RDataFrame

The input is a custom ROOT TTree containing `std::vector` branches. PyROOT's
`RDataFrame` reads it directly, supports lazy one-pass cutflow actions, optional
multithreading, and ROOT skims, and can move to an HPC installation without
changing the event model. The XGBoost tools described in the team slide deck
solve the later diffractive-signal classification problem; they do not improve
this inclusive-Z cutflow or the two-axis ABCD estimate.

The environment also contains `uproot`, `awkward`, `hist`, and `mplhep`. Those
are useful for later plotting and lightweight result inspection, but the first
event pass deliberately has only a PyROOT runtime dependency.

For the BAF CephFS/HTCondor array workflow, including environment packaging,
deterministic input manifests, per-job outputs, merging, and monitoring, see
[`hpc/baf/README.md`](hpc/baf/README.md).

## Environment

Create a separate environment from the repository root:

```bash
conda env create -f feature/z-x-section/environment.yml
conda activate z-xsec
```

The existing `mg5` environment is sufficient for the RDataFrame event pass,
but its Matplotlib installation requires a newer NumPy than the environment
contains. Use the dedicated `z-xsec` environment for the diagnostic renderer.

## Default selection

The configuration is collected near the top of `run_z_selection.py`:

- `IsoMu24` trigger;
- at least two reconstructed muons;
- two muons with `pT > 25 GeV` and `|eta| < 2.4`;
- tight and PF muon flags;
- `|dxy| < 0.05 cm` and `|dz| < 0.10 cm`;
- the two highest-pT quality muons, selected without using charge or isolation;
- at least one of those two muons matched to the trigger object;
- `80 < m(mumu) < 100 GeV`.

After this common selection, the pair is classified as:

| Region | Charge | Isolation |
|---|---|---|
| `SR` | opposite sign | both relative-isolation values below 0.15 |
| `B` | same sign | both relative-isolation values below 0.15 |
| `C` | opposite sign | both relative-isolation values above 0.25 |
| `D` | same sign | both relative-isolation values above 0.25 |

Mixed-isolation and transition events are counted as `unassigned`. The raw
nonprompt estimate is `B*C/D`. It is not yet corrected for prompt leakage, and
charge/isolation closure has not yet been demonstrated.

## Run locally

Bounded data and MC test:

```bash
python feature/z-x-section/run_z_selection.py \
  --sample both --max-files 1 --max-events 100000 \
  --output-dir feature/z-x-section/smoke-tests --label core_smoke
```

Process every locally available file with four RDataFrame threads:

```bash
python feature/z-x-section/run_z_selection.py \
  --sample both --threads 4 --label local_full
```

Only `summary.json` and `cutflow.csv` are written by default. To make a compact
ROOT skim containing assigned and mixed/transition events, add:

```bash
--write-skim --skim-regions SR B C D unassigned
```

Useful controls include `--mass-min`, `--mass-max`, `--isolated-max`,
`--anti-isolated-min`, and `--anti-isolation-mode at_least_one`. Changing the
anti-isolation definition should be treated as a background-method variation,
not silently optimized on the signal-region estimate. `--no-trigger-match` is
provided for diagnostics; disabling `--no-trigger` also disables pair matching.

## Sideband closure and plots

First make a diagnostic skim over the full mass interval, retaining mixed
isolation events so both anti-isolation definitions come from identical input:

```bash
python feature/z-x-section/run_z_selection.py \
  --sample both --threads 4 --mass-min 60 --mass-max 120 \
  --write-skim --output-dir feature/z-x-section/outputs \
  --label local_sideband_input
```

Then calculate the prompt-subtracted closure and render SVG diagnostics:

```bash
python feature/z-x-section/plot_abcd_diagnostics.py \
  feature/z-x-section/outputs/local_sideband_input \
  --output-dir feature/z-x-section/plots \
  --label local_sideband_closure
```

For the final production, replace the provisional peak normalization with an
independently derived prompt-MC normalization using `--prompt-scale VALUE`
(and, if available, `--prompt-scale-uncertainty VALUE`).

The sideband implementation evaluates

```text
A_nonprompt = kappa * (B_data-B_prompt) * (C_data-C_prompt)
                        / (D_data-D_prompt)
```

in low (`60-75 GeV`), signal (`80-100 GeV`), and high (`105-120 GeV`)
regions. A constant `kappa` is estimated from the combined low and high
sidebands. Since the local data subset has no effective luminosity, prompt DY
is normalized iteratively to the OS-isolated peak after the current ABCD
nonprompt subtraction. This is a provisional closure diagnostic, not the final
HPC normalization prescription.

The renderer writes SVG only: nominal/systematic mass overlays, leading and
subleading isolation, raw and prompt-subtracted anti-isolated OS/SS ratios, and
nominal/systematic closure-versus-mass figures. Numerical inputs and fitted
values are stored in `closure_summary.json`.

The interpretation of the current local closure result is recorded in
[`reports/SIDEBAND_ABCD_CLOSURE_STUDY.md`](reports/SIDEBAND_ABCD_CLOSURE_STUDY.md).

## Corrected acceptance and efficiency

`measurement_config.json` holds the truth/reconstruction definitions, null
external-normalization fields, calibration payloads, and explicit provisional
statuses. The default truth object is the opposite-sign stable-muon pair with
a PDG-23 ancestor, dressed with PDG-23-descended stable photons within
`deltaR < 0.1`. Both truth and reconstruction use `80-100 GeV`, `pT > 25 GeV`,
and `|eta| < 2.4` by default.

Certified-data filtering is enabled with the official CMS 2016 legacy JSON
`Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt` stored at the
repository root. Run2016H data naturally restrict it to certified runs
281613-284044. Do not disable this filter for production.

Run the complete local measurement pass with:

```bash
python feature/z-x-section/measure_z_cross_section.py \
  --sample both --threads 4 \
  --output-dir feature/z-x-section/outputs \
  --label local_corrected_measurement
```

Muon momentum corrections are applied before object ordering and all kinematic
cuts. Reconstruction+ID and isolation SFs are multiplied for both muons. The
trigger SF uses the highest-pT selected trigger-matched muon. Supplied
statistical/systematic maps are exposed as coherent up/down variations.

The reconstructed-vertex pileup payload is disabled by default. Its provisional
impact can be evaluated with `--enable-experimental-pileup`. Momentum and SF
corrections can be disabled independently with `--no-muon-momentum` and
`--no-muon-efficiency`.
Alternative definitions can be tested without editing the checked-in config,
for example `--truth-mass-window 50 120 --reco-mass-window 50 120`. Such a
wide reconstructed signal region requires a background model validated over a
still wider control range; the present 60-75/105-120 sidebands would overlap it.

The script reports both interpretations once luminosity is configured:

```text
sigma_fid  = (Ndata - Nbkg) / (efficiency * luminosity)
sigma_full = (Ndata - Nbkg) / (acceptance * efficiency * luminosity)
```

For distributed processing, always add `--defer-normalization`. Each job then
stores raw signed-weight sums and squared-weight sums; global acceptance and
efficiency ratios are formed only by `hpc/baf/merge_measurement_parts.py`.
Per-job acceptance/efficiency values must never be averaged.

## Luminosity and final cross section

The merged measurement payload includes every observed `(run, lumisection)`
pair and writes a compact `processed_lumis.json`. That inventory is the
event-side input to the luminosity calculation; it is not itself an integrated
luminosity. Intersect it with the same golden JSON used during event selection:

```bash
python feature/z-x-section/prepare_luminosity.py \
  merged/measurement_inputs/processed_lumis.json \
  --golden-json PATH_TO_2016_GOLDEN_JSON \
  --output-dir merged/luminosity
```

The output prints the `brilcalc lumi` command. Run that command with the
collaboration-approved normtag, then repeat with `--brilcalc-csv luminosity.csv`
to produce `luminosity_summary.json`. The finalizer accepts either its value via
`luminosity.summary_json` or a reviewed direct value in
`luminosity.integrated_pb_inverse`.

Copy `finalization_config.template.json`, replace every applicable placeholder,
and run:

```bash
python feature/z-x-section/finalize_z_cross_section.py \
  merged/measurement_inputs/measurement_inputs.json \
  --config finalization_config.json \
  --output-dir final/cross_section
```

The finalizer applies luminosity, sample cross sections, generator-weight
denominators, filter efficiencies, k-factors, muon-SF variations, acceptance
and efficiency, ABCD closure uncertainty, and the remaining configured
systematics. Prompt backgrounds have two explicit policies:

- `explicit`: process each available prompt sample and provide its cross
  section/filter efficiency/k-factor/normalization uncertainty;
- `missing_uncertainty`: subtract no prompt central value and propagate a
  documented absolute or fractional missing-component envelope.

Until explicit non-Z prompt samples are available, derive the unresolved
background envelope from the full-statistics closure result:

```bash
python feature/z-x-section/derive_missing_background_envelope.py \
  "$BUDDY/z-xsec/results/RUN_LABEL/plots/full_sideband_closure/closure_summary.json" \
  --output "$BUDDY/z-xsec/results/RUN_LABEL/merged/missing_background_envelope.json"
```

Set `prompt_backgrounds.missing_component.envelope_json` in the finalization
configuration to that output. The envelope is the maximum of the two
anti-isolation definitions, low-only/high-only closure transfers, and
width-scaled one-sided 95% residual bounds. It is a total unresolved-background
model uncertainty, so the same method shifts must not be added again as
independent systematics.

The checked-in signal metadata uses the official CMS Open Data generator value
`2116 pb` for the `ZToMuMu M=50-120 GeV` sample, with unit matching and filter
efficiencies. It normalizes expected-yield and control-region diagnostics only;
the measured cross section remains data-driven. The template has no arbitrary
fractional missing-background default: a full-statistics envelope JSON or
explicit prompt samples are required before interpreting the result.

The smooth sideband diagnostic is run with:

```bash
python feature/z-x-section/measure_z_cross_section.py \
  --sample both --threads 4 --reco-mass-window 60 120 \
  --write-selected-skim --output-dir feature/z-x-section/outputs \
  --label local_corrected_sideband_input

python feature/z-x-section/fit_sideband_background.py \
  feature/z-x-section/outputs/local_corrected_sideband_input \
  --output-dir feature/z-x-section/plots \
  --label local_corrected_functional_background
```

It fits an exponential residual in `60-75` and `105-120 GeV`, interpolates into
`80-100 GeV`, and performs low-to-high/high-to-low holdout tests. It must be
given an external `--prompt-scale` for non-circular production use.
It can read either the original `*_abcd_skim.root` inputs or corrected
`data_selected.root`/`mc_selected.root` files produced by a 60-120 GeV
measurement pass.

## Output interpretation

MC cutflows use signed `genWeight` sums and save `sum_weights_squared` for
statistical propagation. Data use unit weights. The reported ABCD number is a
raw diagnostic rather than a background-subtracted cross-section input.

An arbitrary subset of Run2016H files does not have a known luminosity simply
from its event or file fraction. A final cross-section result should run over
the complete certified data sample on the HPC, or over a demonstrably complete
set of luminosity sections with a separately calculated effective luminosity.
The golden JSON must also be active in `measurement_config.json` during the
event pass; intersecting lumisections after processing cannot remove events
that were already counted.
