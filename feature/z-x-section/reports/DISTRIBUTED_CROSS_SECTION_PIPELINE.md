# Distributed inclusive-Z cross-section pipeline

## Purpose and lineage

This extends the August 2026 local inclusive `Z -> mumu` selection and ABCD
study after the successful BAF pilot (`closure_pilot_v2`). The pilot proved that
one data and one MC file could be selected on a Rocky9 worker and persisted to
CephFS. This iteration makes the worker output sufficient for a later fiducial
and full cross-section calculation without normalizing independent chunks.

## Processing output

Every data job now records its selected `A/B/C/D` counts under both the nominal
two-anti-isolated and alternate at-least-one-anti-isolated definitions, plus
the unique observed run/luminosity sections. Every signal-MC job records raw
signed generator-weight sums and squared sums for the dressed truth denominator,
fiducial numerator, selected truth-matched numerator, muon-SF variations, and
ABCD regions. Optional prompt-background jobs use the same reconstructed
selection and are labeled by process.

`merge_measurement_parts.py` rejects incomplete jobs, duplicate input paths,
pre-normalized MC chunks, and inconsistent configurations. It sums additive
quantities first and only then forms acceptance and efficiency. It also writes
the CMS-style `processed_lumis.json` inventory.

## Final-stage model

`finalize_z_cross_section.py` evaluates

$$
\sigma_\mathrm{fid} =
\frac{N_A-N_\mathrm{nonprompt}-N_\mathrm{prompt}}
{\epsilon\,\mathcal{L}},
\qquad
\sigma_\mathrm{full} =
\frac{N_A-N_\mathrm{nonprompt}-N_\mathrm{prompt}}
{A\,\epsilon\,\mathcal{L}}.
$$

An explicitly processed MC sample is normalized with

$$
w_\mathrm{sample} =
\frac{\mathcal{L}\,\sigma_\mathrm{sample}\,
f_\mathrm{filter}\,f_\mathrm{matching}\,k}{\sum w_\mathrm{gen}}.
$$

Prompt backgrounds can either be included explicitly or represented by a
missing-component uncertainty. The latter is operationally useful while the
team determines which samples exist, but its placeholder size is not a final
physics prescription.

The production configuration now uses the official 2016 legacy certification
JSON and records the CMS Open Data generator metadata for the inclusive signal:
`2116 pb`, matching efficiency 1, and filter efficiency 1. The generator value
is used only for MC normalization diagnostics and optional control-region
leakage subtraction.

`derive_missing_background_envelope.py` implements the interim background
policy. It takes the largest signal-yield excursion from the nominal/alternate
anti-isolation definitions, low/high sideband closure transfers, and one-sided
95% sideband residual bounds. The finalizer can read this result directly with
`missing_component.envelope_json`.

## Local validation

The code was exercised with a complete 668 MB local Run2016H file and a
complete 343 MB local inclusive-Z MC file in the `mg5` environment. The global
merge returned acceptance `0.38155` and efficiency `0.78991` for this one-file
smoke sample. The nominal two-anti-isolated definition had no same-sign
anti-isolated events in one data file; the alternate definition had five and
therefore exercised the complete toy finalization path. A separately labeled
MC file was processed as a synthetic prompt-background role to validate the
explicit-normalization code path, including filter efficiency and k-factor.
These numerical results are software checks only.

The luminosity utility was tested by intersecting 27 observed luminosity
sections and parsing a synthetic BRIL CSV summary. All generated test artifacts
are under `feature/z-x-section/smoke-tests/`, which is excluded from version
control.

## Inputs still required before production interpretation

- verification that the checked-in 2016 legacy golden JSON is the certification
  agreed by the group;
- the approved BRIL normtag and integrated luminosity uncertainty;
- prompt-sample cross sections, generator-weight sums (or confirmation that the
  processed sum is complete), filter efficiencies, and k-factors;
- the final prompt-background sample list or a justified missing-component
  envelope;
- pileup and acceptance-modeling prescriptions and reviewed systematic sizes;
- a statistically supported ABCD closure factor and uncertainty.
