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

tar -czf "$software/single-z-diff-z-xsec.tar.gz" \
    -C "$repo_root" \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='feature/z-x-section/outputs' \
    --exclude='feature/z-x-section/plots' \
    --exclude='feature/z-x-section/smoke-tests' \
    feature/z-x-section \
    feature/muon-efficiency-sfs \
    feature/muon-resolution/calibrations

cp "$repo_root/feature/z-x-section/environment.yml" "$software/environment.yml"
cp "$script_dir/create_environment.sh" "$software/create_environment.sh"
cp "$script_dir/merge_results.sh" "$software/merge_results.sh"
chmod +x "$software/create_environment.sh" "$software/merge_results.sh"
echo "Staged runtime package in $software"
ls -lh "$software/single-z-diff-z-xsec.tar.gz" "$software/environment.yml"
