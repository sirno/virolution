#!/usr/bin/env bash

set -euo pipefail

# flags
DEBUG=0
MODE="${1:-sequential}"

# Get a short commit hash and a timestamp
commit_sha=$(git rev-parse --short HEAD)
timestamp=$(date +%Y%m%d-%H%M%S)

# Create a directory for this commit/timestamp
dest="benchmarks/${commit_sha}/$MODE/${timestamp}"
mkdir -p "$dest"

args="--generations 1000 --n-compartments 2 --disable-progress-bar"

hargs="--warmup 1"

if [ $DEBUG -eq 1 ]; then
  hargs="$hargs --show-output"
fi

cfgs=( "epistasis" "multihost" "singlehost" )

for cfg in "${cfgs[@]}" ; do
  mkdir -p .${cfg}
  hyperfine $hargs --export-json "benchmarks/$cfg.json" "virolution $args \
    --settings tests/${cfg}/${cfg}.yaml \
    --sequence tests/${cfg}/ref.fasta \
    --outdir .${cfg}"
  rm -rf .${cfg}
done

