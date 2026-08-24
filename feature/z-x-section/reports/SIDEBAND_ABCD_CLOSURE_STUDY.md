# Sideband ABCD Closure Study

## Purpose and status

This is the second iteration of the local toy inclusive
`Z/gamma* -> mu+ mu-` cross-section workflow. It adds mass sidebands,
prompt-DY subtraction, a fitted ABCD closure factor, and the requested SVG
diagnostics. The result is a method-development check, not yet a background
estimate suitable for a cross-section result.

The nominal anti-isolated definition requires **both** muons to have relative
isolation above 0.25. Requiring **at least one** anti-isolated muon is retained
as a method variation. In both cases, the isolated selection requires both
muons to have relative isolation below 0.15.

## Regions and estimator

The common reconstruction selection, including pair-level trigger matching, is
documented in the directory README.
After that selection, the four ABCD regions are

| Region | Charge | Isolation |
|---|---|---|
| A | opposite sign | both isolated |
| B | same sign | both isolated |
| C | opposite sign | anti-isolated |
| D | same sign | anti-isolated |

The mass windows are `60-75 GeV` (low sideband), `80-100 GeV` (signal), and
`105-120 GeV` (high sideband). The intervening gaps reduce direct migration
between the peak and sidebands. The prompt-subtracted estimator implemented is

```text
A_nonprompt = kappa * (B_data - B_prompt) * (C_data - C_prompt)
                      / (D_data - D_prompt).
```

Thus, the division by `(D-D_prompt)` applies to the complete product. A
constant `kappa` is fitted using the union of the low and high sidebands.

## Local input and prompt normalization

The diagnostic skim was made from every locally readable file over
`60 < m(mumu) < 120 GeV`. It contains 131,274 selected data events and 459,428
selected prompt-DY MC events across assigned and transition/mixed-isolation
categories. The nominal full-window data counts are A=131,146, B=18, C=89,
and D=21. The small B and D samples dominate the statistical fragility.

The inputs were regenerated after making pair-level trigger matching explicit.
Of 150,901 data events with two quality muons, 150,900 had at least one of the
selected pair marked `muonIsTrigMatched`; all 507,306 corresponding MC events
passed. The single rejected data event lay outside the subsequent mass/ABCD
selection, so all reported regional yields and closure values were unchanged.

The local file subset has no independently established effective luminosity.
For these plots only, prompt DY is therefore normalized iteratively to the
OS-isolated peak after subtracting the current ABCD nonprompt prediction. This
makes the peak comparison partly circular: the peak is not an independent
closure validation in this operating mode. The script also accepts
`--prompt-scale` and `--prompt-scale-uncertainty`, allowing the HPC production
to use a luminosity/cross-section-derived prompt normalization instead.

## Results

The local fits give

| Anti-isolation definition | Prompt scale | fitted kappa | approximate statistical uncertainty |
|---|---:|---:|---:|
| both muons (nominal) | 1.3462e-4 | 20.53 | 10.80 |
| at least one (variation) | 1.3472e-4 | 28.27 | 11.71 |

These closure factors are much larger than one and have roughly 40-50%
statistical uncertainties. They should not be interpreted as measured physics
corrections. They show that the current local ABCD construction does not close
without a very large empirical factor. Likely contributors are the extremely
sparse same-sign controls, peak-normalized prompt-DY shape differences in the
sidebands, missing non-DY prompt backgrounds, and a failure of charge and
isolation to factorize for the selected background mixture.

The `at least one` variation has much greater prompt contamination. In the
signal window it contains 3,551 OS anti-isolated data events, while the scaled
DY estimate is about 3,254, leaving only about 297 after subtraction. Small
changes to prompt modeling can consequently cause large changes in the ABCD
estimate. This supports using the two-anti-isolated definition as nominal and
the looser definition only as a stress test.

## Produced diagnostics

`plot_abcd_diagnostics.py` writes only SVG figures:

- A/B/C/D invariant-mass overlays for each anti-isolation definition, with
  dashed lines at 75, 80, 100, and 105 GeV;
- leading- and subleading-muon relative-isolation spectra;
- raw and prompt-subtracted OS/SS ratios versus mass in anti-isolated events;
- observed prompt-subtracted A and fitted `kappa BC/D` versus mass, with an
  observed/predicted panel, for nominal and variation definitions.

Exact yields, propagated variances, normalization values, and plot names are
stored beside the figures in `closure_summary.json`.

## Recommended next step

Run the same skim over the complete certified data set and all available MC,
using an independently computed luminosity normalization for prompt DY and
adding top, diboson, W+jets, and `Z -> tautau` prompt backgrounds where
available. Refit `kappa` separately in the low and high sidebands and test
their compatibility. Only if the sidebands show stable closure should one use
a constant factor in the signal window; otherwise the next iteration should
fit a simple mass-dependent background/transfer model and validate it in
additional control regions.
