#!/usr/bin/env bash
# run_law.sh <rv2dec> <t2t.fa> <v1 dir> <anchors.txt> <rv2 dir> <names.txt> <requests.tsv> <truth.tsv> <out dir> [label]
# Window law (model B, hw-apex protocol v1.1 rules) of the refrel v2 prototype: one thread (taskset -c 3), no nice,
# silence gate before the run (load average (1 min) < 0.5 AND no process above 5 % CPU); the gate state, the input
# hashes, the law log and the raw per-request timings are written to <out dir>.
# requests: research/refrel_v2/law_requests.py (seed 20261003, 10 000 per W, W = 1 / 1024 / 8192 / 65536);
# truth: `rv2dec lawtruth`, from the v1 archives fully decoded with every v1 check (independent of the v2 code).
set -euo pipefail
B=${1:?}; REF=${2:?}; V1=${3:?}; AN=${4:?}; AR=${5:?}; NM=${6:?}; RQ=${7:?}; TR=${8:?}; L=${9:?}; LB=${10:-n558}
mkdir -p "$L"; A=$(paste -sd, "$AN")
sha256sum "$B" "$RQ" "$TR" "$NM" "$AN" > "$L/${LB}_inputs.sha256"
until awk '{exit !($1<0.5)}' /proc/loadavg && [ -z "$(ps -eo pcpu= --sort=-pcpu | awk '$1>5' | head -1)" ]; do sleep 30; done
{ echo "# silence gate ${LB} $(date -u +%FT%TZ)"; echo "loadavg: $(cat /proc/loadavg)"; echo "mem: $(free -g | awk '/Mem:/{print $7" GB available of "$2}')"
  echo "swap: $(free -g | awk '/Swap:/{print $3" GB used of "$2}')"; echo "top processes by CPU:"; ps -eo pid,user,ni,pcpu,pmem,comm --sort=-pcpu | head -8; } > "$L/silence_${LB}.txt"
{ echo "# ${LB} start $(date -u +%FT%TZ) loadavg $(cat /proc/loadavg)"
  LAW_RAW=$L/raw_${LB}.tsv /usr/bin/time -v taskset -c 3 "$B" law rv2 "$REF" "$V1" "$A" "$AR" "$NM" "$RQ" "$TR" || echo "# law exit $?"
  echo "# end $(date -u +%FT%TZ) loadavg $(cat /proc/loadavg)"; } > "$L/${LB}_law.log" 2>&1
gzip -9 -f "$L/raw_${LB}.tsv"
grep -E "^LAW|^DQ_RESULT|^WIN|Maximum resident" "$L/${LB}_law.log"
