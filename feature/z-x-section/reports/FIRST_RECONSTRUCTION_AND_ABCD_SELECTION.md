# First Inclusive Z Reconstruction and ABCD Selection

## Purpose and project context

This is the first local implementation of the fast inclusive
`Z/gamma* -> mu+ mu-` cross-section exercise proposed by Matthias. Its purpose
is to freeze and inspect the reconstructed event selection and the proposed
charge/isolation ABCD regions before adding generator-level acceptance,
truth-matched efficiency, luminosity, and a final cross-section calculation.

The eventual inclusive measurement is distinct from the diffractive-Z signal
extraction. Inclusive Drell--Yan is the signal here. The ABCD construction is a
candidate estimate of the nonprompt/QCD contribution; it does not estimate the
inclusive Drell--Yan background to the separate diffractive analysis.

## Implementation choice

`run_z_selection.py` uses Python/PyROOT with `RDataFrame` and small declared C++
helpers for vector-branch object selection. This was chosen because the input
is already a custom ROOT TTree and because the same lazy event graph can run
locally or with implicit multithreading on an HPC system. The XGBoost workflow
in the team slides addresses diffractive topology discrimination and is not
needed for this cutflow.

`environment.yml` defines a separate `z-xsec` Conda environment containing
ROOT plus `uproot`, `awkward`, `hist`, and `mplhep` for later inspection and
plotting. Only PyROOT is required by the first event counter. The implementation
was tested in the existing `mg5` environment with Python 3.11.15 and ROOT
6.40.02.

## Common reconstruction selection

The common selection is applied before looking at charge or isolation:

1. `IsoMu24` passes.
2. At least two reconstructed muons exist.
3. At least two muons have `pT > 25 GeV` and `|eta| < 2.4`.
4. At least two additionally pass tight ID, the PF-muon flag,
   `|dxy| < 0.05 cm`, and `|dz| < 0.10 cm`.
5. The two highest-pT quality muons are chosen without using their charge or
   isolation.
6. At least one of the selected muons matches the `IsoMu24` trigger object.
7. Their reconstructed mass satisfies `80 < m(mumu) < 100 GeV`.

The pT, eta, charge, and mass choices follow the event selection summarized in
the August team slides. Tight/PF and impact-parameter requirements provide a
clean prompt reconstruction definition and match the earlier local muon study.
They remain configurable and should be aligned with the precise definitions of
Matthias's efficiency scale factors before the final measurement.

## Isolation and ABCD definition

Relative isolation is the delta-beta corrected PF definition

`(charged + max(0, neutral + photon - 0.5*pileup))/pT`.

The default isolated requirement is below 0.15 for both muons. The default
anti-isolated requirement is above 0.25 for both muons. The assigned regions
are:

| Region | Charge | Isolation |
|---|---|---|
| SR | opposite sign | both isolated |
| B | same sign | both isolated |
| C | opposite sign | both anti-isolated |
| D | same sign | both anti-isolated |

Mixed and transition events are counted but not used in the raw transfer
factor. The nominal raw estimate is `N_B*N_C/N_D`. Region name `SR` is used
instead of `A` to avoid confusion with detector acceptance.

The script also exposes `--anti-isolation-mode at_least_one` as a predefined
method variation. It must not be chosen solely because it gives a preferred
signal-region estimate; closure and composition studies should decide which
definition is defensible.

## Outputs

Every run writes:

- `summary.json`: exact configuration, input files, rejected files, unweighted
  and weighted cutflows, ABCD counts, sum of weights squared, and the raw
  estimate;
- `cutflow.csv`: a compact table with step and cumulative efficiencies;
- optionally, one compact ROOT skim containing run/lumi/event identifiers,
  pair kinematics, charge, isolation, mass, event weight, and ABCD region code.

MC yields use signed `genWeight`. The cutflow records both raw-event and
signed-weight efficiencies. A file with an incompatible schema is rejected
before constructing the chain rather than silently assigned unit weight.

## Smoke tests

The bounded test under `smoke-tests/core_smoke_r2` processed up to 100,000
entries from one file per sample. It selected 2,892 data and 14,671 MC events in
SR. The small test had no populated same-sign region, so its ABCD estimate is
correctly recorded as undefined.

`smoke-tests/skim_smoke` separately checked the lazy ROOT snapshot. It wrote
292 selected entries and all requested event, pair, isolation, and region-code
branches.

## Full local-input diagnostic

The four-thread local pass completed in about six seconds of reported event-loop
wall time.

### Data

- 31 files and 4,389,502 input events;
- 2,813,354 pass `IsoMu24`;
- 250,872 contain two kinematic muons;
- 150,901 contain two quality muons;
- 130,769 enter the mass window;
- 120,314 enter SR.

With the nominal two-anti-isolated definition:

| Region | Events |
|---|---:|
| SR | 120,314 |
| B | 3 |
| C | 39 |
| D | 4 |
| mixed/transition | 10,409 |

The raw estimate is `29.25 +/- 22.83 (stat.)` events, approximately 0.024% of
SR. This is not yet a background prediction: it has no prompt subtraction or
closure uncertainty, and its apparent precision is controlled by three B and
four D events.

As a diagnostic only, requiring at least one anti-isolated muon gives
`B=3`, `C=3551`, and `D=95`, corresponding to a raw estimate of
`112.14 +/- 65.78 (stat.)`. It improves the anti-isolated statistics but does
not fix the B-region limitation.

### Inclusive-DY simulation

Thirty-one schema-compatible MC files containing 1,496,000 entries were used.
The unweighted SR yield is 423,759, corresponding to 28.33% of all input
events; the signed-weight cumulative efficiency is 28.16%.

One local file has no valid `pfExtractor/pfTree`, and six readable files lack
`genWeight`; all seven are explicitly listed as rejected in the summary. The
missing-weight files must be regenerated or replaced before they can
contribute to weighted acceptance and efficiency. They should not be mixed
with the weighted files using an assumed unit weight.

The inclusive-DY simulation has only two B events and no D events, so its ABCD
estimate is undefined. That is expected qualitatively: this sample models the
prompt signal, not the nonprompt/QCD component for which the ABCD construction
is intended.

## Limitations and next steps

The current result is a selection and counting diagnostic. Before using the
ABCD number in a cross section:

1. study charge/isolation factorization in mass sidebands;
2. quantify prompt-Z leakage in the anti-isolated regions and prompt processes
   in the same-sign region;
3. consider a loose-not-tight or fake-factor method if the isolated same-sign
   region remains statistically limiting;
4. add invariant-mass and relative-isolation plots for SR/B/C/D;
5. define Born or dressed generator muons and calculate signed-weight
   acceptance and truth-matched efficiency;
6. confirm the reconstruction, isolation, and event-trigger interpretation of
   Matthias's scale factors;
7. run the frozen counter over the complete certified Run2016H data on the HPC
   and calculate its trigger-compatible integrated luminosity.

The arbitrary 31-file data subset is sufficient for implementation and
qualitative background diagnostics. It must not be normalized using the full
Run2016H luminosity.
