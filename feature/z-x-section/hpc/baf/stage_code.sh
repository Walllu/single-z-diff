#!/bin/bash
# Package only runtime code/payloads into persistent BAF CephFS storage.
set -euo pipefail

if [[ -z "${BUDDY:-}" ]]; then
    echo "ERROR: BUDDY is not set" >&2
    exit 2
fi
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../../../.." && pwd)"
software="$BUDDY/z-xsec/software"
mkdir -p "$software"

archive="$software/single-z-diff-z-xsec.tar.gz"
archive_tmp="$software/.single-z-diff-z-xsec.tar.gz.$$"
trap 'rm -f "$archive_tmp"' EXIT

# BAF home storage can update directory metadata while GNU tar traverses it,
# producing "file changed as we read it" and exit status 1 even when the
# regular files are stable. Treat that condition as non-fatal, then explicitly
# validate the archive and the runtime payloads before publishing it.
if tar --help 2>&1 | grep -q -- '--ignore-failed-read'; then
    tar --ignore-failed-read --warning=no-file-changed -czf "$archive_tmp" \
        -C "$repo_root" \
        --exclude='__pycache__' \
        --exclude='*.pyc' \
        --exclude='*.root' \
        --exclude='feature/z-x-section/outputs' \
        --exclude='feature/z-x-section/plots' \
        --exclude='feature/z-x-section/smoke-tests' \
        --exclude='feature/z-x-section/hpc/CephFSData' \
        --exclude='feature/z-x-section/hpc/dataset*' \
        --exclude='feature/z-x-section/hpc/dump.txt' \
        --exclude='feature/z-x-section/hpc/missing_from_cephfs.txt' \
        --exclude='feature/z-x-section/hpc/replace_on_cephfs.txt' \
        feature/z-x-section \
        feature/muon-efficiency-sfs \
        feature/muon-resolution/calibrations \
        Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
else
    tar -czf "$archive_tmp" \
        -C "$repo_root" \
        --exclude='__pycache__' \
        --exclude='*.pyc' \
        --exclude='*.root' \
        --exclude='feature/z-x-section/outputs' \
        --exclude='feature/z-x-section/plots' \
        --exclude='feature/z-x-section/smoke-tests' \
        --exclude='feature/z-x-section/hpc/CephFSData' \
        --exclude='feature/z-x-section/hpc/dataset*' \
        --exclude='feature/z-x-section/hpc/dump.txt' \
        --exclude='feature/z-x-section/hpc/missing_from_cephfs.txt' \
        --exclude='feature/z-x-section/hpc/replace_on_cephfs.txt' \
        feature/z-x-section \
        feature/muon-efficiency-sfs \
        feature/muon-resolution/calibrations \
        Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
fi

for member in \
    feature/z-x-section/run_z_selection.py \
    feature/z-x-section/measure_z_cross_section.py \
    feature/z-x-section/finalize_z_cross_section.py \
    feature/z-x-section/derive_missing_background_envelope.py \
    feature/z-x-section/hpc/baf/merge_measurement_parts.py \
    Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt \
    feature/z-x-section/environment.yml \
    feature/muon-efficiency-sfs/MuonPerformance_EfficiencyCorrections.py \
    feature/muon-resolution/calibrations/muon_momentum_2016H_bb_be_ee.json
do
    tar -xOf "$archive_tmp" "$member" >/dev/null
done
mv "$archive_tmp" "$archive"
trap - EXIT

cp "$repo_root/feature/z-x-section/environment.yml" "$software/environment.yml"
cp "$script_dir/create_environment.sh" "$software/create_environment.sh"
cp "$script_dir/merge_results.sh" "$software/merge_results.sh"
chmod +x "$software/create_environment.sh" "$software/merge_results.sh"
echo "Staged runtime package in $software"
ls -lh "$archive" "$software/environment.yml"
