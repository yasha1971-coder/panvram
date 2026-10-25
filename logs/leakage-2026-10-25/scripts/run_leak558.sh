#!/usr/bin/env bash
# Wait for the window-law run to finish (its latencies must not be disturbed), then measure leakage on all 558
# assemblies: split by individual against split by locus, at W = 8 KiB and 32 KiB, and the accessor self-test on
# ranges that cover a non-ACGT run.
set -uo pipefail
D=/tmp/claude-1002/l2
V=/tmp/claude-1002/v3b
W1=$HOME/pubrepo/.wk/rv2b
REF=$HOME/golden/genome/t2t.fa
V1=$HOME/pubrepo/.wk/hprc_cohort/v1
until grep -q LAW_DONE "$V/chain.log" 2>/dev/null; do sleep 120; done
echo "LAW_DONE seen, starting the leakage runs"

echo "=== accessor self-test aimed at non-ACGT runs"
"$V/rv2dec_packed" accself "$REF" "$V1" "$(paste -sd, "$W1/anchors558_k8.txt")" 1200 20261003 2>/dev/null | grep ACCSELF

for W in 8192 32768; do
  for MODE in people locus; do
    echo "=== W $W mode $MODE"
    "$D/locate_leak" leak "$REF" "$V1" "$W1/names558.txt" "$W" 16 "$MODE" 20261003 10 2>/dev/null | grep '^LEAK'
  done
done
echo LEAK558_DONE
