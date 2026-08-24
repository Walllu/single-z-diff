#!/bin/bash
# Run once inside a Rocky9 interactive CPU job after stage_code.sh.
set -euo pipefail

if [[ -z "${BUDDY:-}" ]]; then
    echo "ERROR: BUDDY is not set" >&2
    exit 2
fi
source /etc/profile >/dev/null 2>&1 || true
module load miniforge/24.7.1-0-py312

software="$BUDDY/z-xsec/software"
env_dir="/jwd/z-xsec-build-env"
tarball="$software/z-xsec-env.tar.gz"
if [[ -e "$env_dir" ]]; then
    echo "ERROR: build directory already exists: $env_dir" >&2
    exit 2
fi
conda env create --prefix "$env_dir" --file "$software/environment.yml" --yes
conda run --prefix "$env_dir" python -c 'import ROOT, numpy, uproot, mplhep; print(ROOT.gROOT.GetVersion())'
conda run --prefix "$env_dir" conda-pack --prefix "$env_dir" --output "$tarball" --force
echo "Created $tarball"
ls -lh "$tarball"
