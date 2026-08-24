# Running the inclusive-Z closure workflow on BAF

This directory packages the current CPU-only `Z -> mumu` sideband/ABCD pass as
a reproducible HTCondor array. It follows the local BAF rules:

- submit jobs from `$HOME`, never from `/cephfs`;
- keep large input and persistent output under `$BUDDY`;
- copy packaged code and the environment into `/jwd` for each job;
- request no GPU and use the Rocky9 container;
- request high CephFS I/O for the ROOT event pass;
- give every job a deterministic manifest and a collision-free output path.

The worker expects the existing custom ROOT format containing
`pfExtractor/pfTree`. CMS MINIAOD or NanoAOD files cannot be passed directly to
this workflow. They must first be converted with the collaboration extractor,
or the event reader must be adapted to their schema.

## Persistent layout

The scripts create the following layout:

```text
$BUDDY/z-xsec/
  software/                         packaged code and relocatable environment
  results/RUN_LABEL/
    parts/data/part_CLUSTER_PROC/   independent data outputs
    parts/mc/part_CLUSTER_PROC/     independent MC outputs
    logs/                            worker stdout/stderr
    merged/full_sideband_input/     merged ROOT skims and cutflow
    plots/full_sideband_closure/    closure JSON and SVG plots
```

## 1. Put inputs on CephFS

Choose stable directories containing only the intended custom ROOT inputs, for
example:

```bash
export ZXS_DATA="$BUDDY/z-xsec/inputs/2016H"
export ZXS_MC="$BUDDY/z-xsec/inputs/ZToMuMu_M-50To120"
```

Personal `$BUDDY` space is documented as 500 GB. A complete extracted Run2016H
sample may exceed that, so use an approved collaboration/shared CephFS area for
the large inputs if one exists. The manifest driver accepts any readable
absolute CephFS path; the inputs do not have to live below personal `$BUDDY`.

### Compare a partial CephFS copy with a complete source

Create an inventory wherever the known-complete source is visible:

```bash
python feature/z-x-section/hpc/baf/compare_file_inventories.py scan \
  /path/to/complete/2016H \
  --output complete_2016H.json
```

On BAF, inventory the candidate CephFS copy:

```bash
python feature/z-x-section/hpc/baf/compare_file_inventories.py scan \
  /cephfs/user/tsaala/Hackathon/data/Real/2016H \
  --output "$HOME/cephfs_2016H.json"
```

Place both small JSON files on the same machine and compare them, taking the
known-complete inventory as the reference:

```bash
python feature/z-x-section/hpc/baf/compare_file_inventories.py compare \
  complete_2016H.json "$HOME/cephfs_2016H.json" \
  --output-dir "$HOME/2016H_inventory_comparison"
```

The comparison matches unique ROOT basenames and requires exact byte sizes. It
writes separate missing, size-mismatch, and extra lists. The combined
`copy_or_replace_from_reference.txt` is suitable as an `rsync --files-from`
input when the complete source is mounted as a filesystem. Inspect the JSON
summary before copying; the utility refuses ambiguous duplicate basenames.

Large downloads should be performed from the designated BAF login node in a
`screen` session, not from an HTCondor worker. Apply the certified-run JSON in
the final measurement workflow; a directory or event fraction does not define
an integrated luminosity.

## 2. Stage code and build the environment once

From a clone under `$HOME` on the BAF login node:

```bash
cd ~/single-z-diff
bash feature/z-x-section/hpc/baf/stage_code.sh
```

Start a short CPU-only interactive Rocky9 job from `$HOME`:

```bash
condor_submit -interactive \
  -append '+ContainerOS = "Rocky9"' \
  -append '+CephFS_IO = "medium"' \
  -append '+MaxRuntimeHours = 2' \
  -append 'request_gpus = 0' \
  -append 'request_cpus = 4' \
  -append 'request_memory = 8000 MB'
```

Inside it, build and package the conda environment:

```bash
bash "$BUDDY/z-xsec/software/create_environment.sh"
exit
```

Re-run `stage_code.sh` whenever analysis code or payloads change. Rebuild the
environment only when `environment.yml` changes.

## 3. Prepare and submit the array

Run this from the repository under `$HOME`:

```bash
python feature/z-x-section/hpc/baf/prepare_submission.py \
  --data-dir "$ZXS_DATA" \
  --mc-dir "$ZXS_MC" \
  --submit-dir "$HOME/z-xsec-submit/full_sideband_v1" \
  --run-label full_sideband_v1 \
  --data-files-per-job 10 \
  --mc-files-per-job 10

cd "$HOME/z-xsec-submit/full_sideband_v1"
condor_submit submit.jdl
```

The generated JDL requests four CPUs, 8 GB RAM, no GPU, Rocky9, two hours, and
high CephFS I/O. Do not increase concurrency blindly: start with these bundles,
measure throughput, and avoid saturating shared CephFS.

Useful monitoring commands are:

```bash
condor_q -all -global
condor_q -better-analyze CLUSTER_ID
condor_history "$USER"
```

Worker logs are persistent under
`$BUDDY/z-xsec/results/full_sideband_v1/logs`. A successful part contains a
`SUCCESS` marker. Held jobs can safely be released because the same Condor
cluster/process identifiers overwrite their own output directory.

## 4. Check completeness, merge, and plot

Do not merge while array jobs are still running. Count expected jobs from the
preparation output and compare against `SUCCESS` markers:

```bash
find "$BUDDY/z-xsec/results/full_sideband_v1/parts" -name SUCCESS | wc -l
```

The stronger check compares every expected manifest path with the successful
summaries and detects missing, unexpected, or duplicated files:

```bash
python ~/single-z-diff/feature/z-x-section/hpc/baf/check_status.py \
  "$HOME/z-xsec-submit/full_sideband_v1"
```

It exits successfully only when the complete expected file set is covered
exactly once.

Start another short Rocky9 CPU interactive job and run:

```bash
bash "$BUDDY/z-xsec/software/merge_results.sh" full_sideband_v1
```

The merger rejects duplicate input-file paths, which prevents accidentally
double counting a separately resubmitted job. It sums event counts, signed
weights, squared weights, cutflows, and ABCD regions, then uses `hadd` for the
ROOT skims. Finally it regenerates the closure diagnostics.

The default closure plot still uses the provisional peak-derived DY scale. If
an independently derived prompt normalization is available, pass it as the
second argument:

```bash
bash "$BUDDY/z-xsec/software/merge_results.sh" full_sideband_v1 PROMPT_SCALE
```

## Scope of this first driver

This array deliberately prioritizes the full-statistics sideband closure pass.
The acceptance/efficiency script also accepts `--data-files-from` and
`--mc-files-from`, but its distributed MC summaries should not be combined
until a global generator-weight denominator and luminosity normalization are
frozen. Per-chunk normalization would be incorrect. Once closure is available,
the corrected measurement can either run as one larger RDataFrame job or gain a
separate reducer that combines its unnormalized sufficient statistics.
