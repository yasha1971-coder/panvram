#!/usr/bin/env bash
# run_window_cost.sh <label> <reference.fa> <archive dir> <names.txt> <requests.tsv> <truth.tsv> <out dir>
# One measurement, one thread (taskset -c 3), no nice, after a silence gate: load average (1 min) < 0.5 and no process
# using > 5 % CPU. Written next to the result: <label>.log (window_cost output + /usr/bin/time -v), raw_<label>.tsv.gz
# (every timed request), silence_<label>.txt (host state at the gate), <label>_inputs.sha256 (inputs and the binary).
# Paths are recorded as given; call it with paths relative to the current directory to keep the log free of home paths.
set -euo pipefail
LABEL=${1:?}; REF=${2:?}; DIR=${3:?}; NAMES=${4:?}; REQ=${5:?}; TRUTH=${6:?}; OUT=${7:?}
H=$(cd "$(dirname "$0")" && pwd); BIN=${WINDOW_COST_BIN:-$H/window_cost}; CPU=${WINDOW_COST_CPU:-3}
mkdir -p "$OUT"
sha256sum "$REQ" "$TRUTH" "$NAMES" "$BIN" > "$OUT/${LABEL}_inputs.sha256"
until awk '{exit !($1<0.5)}' /proc/loadavg && [ -z "$(ps -eo pcpu= --sort=-pcpu | awk '$1>5' | head -1)" ]; do sleep 30; done
{ echo "# silence gate $LABEL $(date -u +%FT%TZ)"; echo "loadavg: $(cat /proc/loadavg)"; echo "df: $(df -h / | tail -1 | awk '{print $4" free of "$2}')"
  echo "mem: $(free -g | awk '/Mem:/{print $7" GB available of "$2}')"; echo "top processes by CPU:"; ps -eo pid,user,ni,pcpu,pmem,comm --sort=-pcpu | head -8; } > "$OUT/silence_$LABEL.txt"
{ echo "# $LABEL start $(date -u +%FT%TZ) loadavg $(cat /proc/loadavg)"
  WINDOW_COST_RAW="$OUT/raw_$LABEL.tsv" /usr/bin/time -v taskset -c "$CPU" "$BIN" run "$REF" "$DIR" "$NAMES" "$REQ" "$TRUTH" || echo "# run exit $?"
  echo "# end $(date -u +%FT%TZ) loadavg $(cat /proc/loadavg)"; } > "$OUT/$LABEL.log" 2>&1
gzip -9 -f "$OUT/raw_$LABEL.tsv"
grep -E "^WINDOW_COST|^DQ_RESULT|^WIN|^RESULT|Maximum resident|Elapsed" "$OUT/$LABEL.log"
