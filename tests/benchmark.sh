#!/usr/bin/env bash

DEBUG=0
ARGS="--generations 1000 --n-compartments 2 --disable-progress-bar"
HARGS="--warmup 1"
if [ $DEBUG -eq 1 ]; then
  HARGS="$HARGS --show-output"
fi

cfgs=( "epistasis" "multihost" "singlehost" )
mkdir -p benchmarks

for cfg in "${cfgs[@]}" ; do
  mkdir -p .${cfg}
  hyperfine $HARGS --export-json "benchmarks/$cfg.json" "virolution $ARGS \
    --settings tests/${cfg}/${cfg}.yaml \
    --sequence tests/${cfg}/ref.fasta \
    --outdir .${cfg}"
  rm -rf .${cfg}
done

