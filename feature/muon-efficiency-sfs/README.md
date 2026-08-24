# Collaborator-provided muon efficiency scale factors

These reconstruction, trigger, and isolation efficiency scale-factor maps were
provided by a professor collaborating on the analysis. They are copied here
verbatim and remain under the provider's scientific purview. They are distinct
from the BB/BE/EE momentum scale and resolution correction.

Before production use, confirm the intended validity period, sample, pT/eta
binning, and treatment outside the calibrated range with the provider.

The trigger and isolation maps have shape 10 x 15, matching the declared ten
eta and fifteen pT bins. The reconstruction maps have shape 10 x 10. Their
final column is identically 1.0 in the nominal map and 0.0 in both uncertainty
maps for every eta row. As an explicit provisional overflow policy, the helper
repeats that final reconstruction column through the five missing high-pT bins
(`75-100 GeV`). The collaborator-provided text files remain unchanged. The
loader validates all dimensions and resolves input files relative to the
module, so lookups do not depend on the caller's working directory.
