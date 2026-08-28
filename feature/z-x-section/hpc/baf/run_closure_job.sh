#!/bin/bash
# CPU-only BAF worker for closure skims and additive measurement inputs.
set -euo pipefail

sample="${1:?sample is required}"
manifest="${2:?manifest is required}"
cluster_id="${3:?cluster id is required}"
process_id="${4:?process id is required}"
run_label="${5:?run label is required}"

if [[ "$sample" != "data" && "$sample" != "mc" && "$sample" != background__* ]]; then
    echo "ERROR: sample must be data, mc, or background__NAME; got $sample" >&2
    exit 2
fi
if [[ -z "${BUDDY:-}" ]]; then
    echo "ERROR: BUDDY is not available inside the job" >&2
    exit 2
fi

job_id="${cluster_id}_${process_id}"
software="$BUDDY/z-xsec/software"
result_root="$BUDDY/z-xsec/results/$run_label"
part_root="$result_root/parts/$sample"
output_label="part_$job_id"
log_root="$result_root/logs"
mkdir -p "$part_root" "$log_root"
rm -f "$part_root/$output_label/SUCCESS"
out="$log_root/out.${sample}.${job_id}.log"
err="$log_root/err.${sample}.${job_id}.log"
exec 1>>"$out" 2>>"$err"

echo "[$(date --iso-8601=seconds)] starting $sample job $job_id on $(hostname)"
echo "manifest=$manifest"
echo "BUDDY=$BUDDY"

workdir="/jwd/z-xsec-$job_id"
source_dir="$workdir/source"
env_dir="$workdir/env"
mkdir -p "$source_dir" "$env_dir"

tar -xzf "$software/single-z-diff-z-xsec.tar.gz" -C "$source_dir"
tar -xzf "$software/z-xsec-env.tar.gz" -C "$env_dir"

# Conda's generated activation script probes variables such as CONDA_PREFIX
# before defining them. Temporarily disable nounset while activating and
# relocating the conda-pack environment, then restore strict worker settings.
set +u
source "$env_dir/bin/activate"
conda-unpack
set -u

input_option="--data-files-from"
measurement_sample="data"
mc_role="signal"
process_name="data"
if [[ "$sample" == "mc" ]]; then
    input_option="--mc-files-from"
    measurement_sample="mc"
    process_name="z_to_mumu"
elif [[ "$sample" == background__* ]]; then
    input_option="--mc-files-from"
    measurement_sample="mc"
    mc_role="background"
    process_name="${sample#background__}"
fi

if [[ "$sample" == "data" || "$sample" == "mc" ]]; then
    python -u "$source_dir/feature/z-x-section/run_z_selection.py" \
        --sample "$sample" \
        "$input_option" "$PWD/$manifest" \
        --threads 4 \
        --mass-min 60 --mass-max 120 \
        --write-skim --skim-regions SR B C D unassigned \
        --output-dir "$part_root" \
        --label "$output_label"
fi

python -u "$source_dir/feature/z-x-section/measure_z_cross_section.py" \
    --config "$source_dir/feature/z-x-section/measurement_config.json" \
    --sample "$measurement_sample" \
    "$input_option" "$PWD/$manifest" \
    --mc-role "$mc_role" \
    --process-name "$process_name" \
    --defer-normalization \
    --threads 4 \
    --output-dir "$part_root" \
    --label "$output_label"

touch "$part_root/$output_label/SUCCESS"
echo "[$(date --iso-8601=seconds)] completed $sample job $job_id"
