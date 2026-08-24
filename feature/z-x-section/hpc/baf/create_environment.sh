#!/bin/bash
# Run once inside a Rocky9 interactive CPU job after stage_code.sh.
set -euo pipefail

if [[ -z "${BUDDY:-}" ]]; then
    echo "ERROR: BUDDY is not set" >&2
    exit 2
fi

# Do not source /etc/profile here: BAF's profile may terminate a non-login
# batch script before the environment build begins. Interactive Rocky9 jobs
# normally export the module function; otherwise initialize Lmod directly.
if ! type module >/dev/null 2>&1; then
    for module_init in \
        /etc/profile.d/modules.sh \
        /etc/profile.d/lmod.sh \
        /usr/share/lmod/lmod/init/bash
    do
        if [[ -r "$module_init" ]]; then
            source "$module_init"
            break
        fi
    done
fi
if ! type module >/dev/null 2>&1; then
    echo "ERROR: the BAF module command is unavailable" >&2
    exit 2
fi

echo "Loading the BAF miniforge module"
module load miniforge/24.7.1-0-py312

software="$BUDDY/z-xsec/software"
env_dir="/jwd/z-xsec-build-env"
tarball="$software/z-xsec-env.tar.gz"
if [[ -e "$env_dir" ]]; then
    echo "ERROR: build directory already exists: $env_dir" >&2
    exit 2
fi
echo "Creating the environment in $env_dir"
conda env create --prefix "$env_dir" --file "$software/environment.yml" --yes
echo "Running the import smoke test"
conda run --prefix "$env_dir" python -c 'import ROOT, numpy, uproot, mplhep; print(ROOT.gROOT.GetVersion())'
echo "Packing the relocatable environment"
conda run --prefix "$env_dir" conda-pack --prefix "$env_dir" --output "$tarball" --force
echo "Created $tarball"
ls -lh "$tarball"
