#!/bin/bash
# Run inside a BAF Rocky9 interactive CPU job after every array job succeeds.
set -euo pipefail

run_label="${1:?usage: merge_results.sh RUN_LABEL [PROMPT_SCALE]}"
prompt_scale="${2:-}"
if [[ -z "${BUDDY:-}" ]]; then
    echo "ERROR: BUDDY is not set" >&2
    exit 2
fi

software="$BUDDY/z-xsec/software"
result_root="$BUDDY/z-xsec/results/$run_label"
workdir="/jwd/z-xsec-merge-$run_label-$$"
mkdir -p "$workdir/source" "$workdir/env"
tar -xzf "$software/single-z-diff-z-xsec.tar.gz" -C "$workdir/source"
tar -xzf "$software/z-xsec-env.tar.gz" -C "$workdir/env"
source "$workdir/env/bin/activate"
conda-unpack

python -u "$workdir/source/feature/z-x-section/hpc/baf/merge_abcd_parts.py" \
    "$result_root/parts" \
    --output-dir "$result_root/merged" \
    --label full_sideband_input

plot_args=()
if [[ -n "$prompt_scale" ]]; then
    plot_args+=(--prompt-scale "$prompt_scale")
fi
python -u "$workdir/source/feature/z-x-section/plot_abcd_diagnostics.py" \
    "$result_root/merged/full_sideband_input" \
    --output-dir "$result_root/plots" \
    --label full_sideband_closure \
    "${plot_args[@]}"

echo "Merged skims: $result_root/merged/full_sideband_input"
echo "Closure outputs: $result_root/plots/full_sideband_closure"
