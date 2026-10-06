#!/bin/bash
# Alternating upstream/fast rounds. Run from a work directory that holds
# data/ (from make_data.R) and out/. Libraries: LIB_BASE (upstream susieR)
# and LIB_FAST (this fork), each an R library with susieR installed.
# Usage: bench.sh "<cases>" [rounds]   e.g. bench.sh "rss_p1000 ind_p1000" 5
cases="$1"; rounds="${2:-3}"
LIB_BASE="${LIB_BASE:-../lib_base}"; LIB_FAST="${LIB_FAST:-../lib_fast}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p out
echo "round,case,build,seconds,niter" > out/bench.csv
for r in $(seq 1 $rounds); do
  for c in $cases; do
    for build in base fast; do
      lib=$LIB_BASE; [ $build = fast ] && lib=$LIB_FAST
      line=$(Rscript "$HERE/bench_one.R" "$lib" $c 2>/dev/null | grep RESULT)
      set -- $line
      echo "$r,$c,$build,$4,$5" | tee -a out/bench.csv
    done
  done
done
