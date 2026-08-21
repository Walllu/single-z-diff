# Collaborator-provided muon efficiency scale factors

These reconstruction, trigger, and isolation efficiency scale-factor maps were
provided by a professor collaborating on the analysis. They are copied here
verbatim and remain under the provider's scientific purview. They are distinct
from the BB/BE/EE momentum scale and resolution correction.

Before production use, confirm the intended validity period, sample, pT/eta
binning, and treatment outside the calibrated range with the provider.

Known integration check: the trigger and isolation maps have shape 10 x 15,
while the reconstruction maps have shape 10 x 10. The supplied helper defines
15 pT bins for every map, so high-pT reconstruction lookups currently raise an
`IndexError`. The helper also resolves input text files relative to the process
working directory. These points are documented here rather than silently
changing collaborator-owned inputs.
