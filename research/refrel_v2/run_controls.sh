#!/usr/bin/env bash
# run_controls.sh <aceapex checkout at e902934> <t2t.fa> <v1 dir> <anchors.txt> <rv2 dir> <names.txt> <MANIFEST.tsv>
#                 <requests.tsv> <truth.tsv> <out dir>
# The judge of the v2 decoder on a subset, run with the clean decoder and with each of the five deliberate faults of
# build.sh, one after another: (1) `rv2dec judge`: every assembly decoded in full, FASTA bytes + FASTA XXH3 == MANIFEST.tsv;
# (2) `rv2dec law`: the window law with every answer's SHA-256 == truth and every full decode's XXH3 == truth.
# Expected: clean PASS in both; every fault FAIL in at least one, written to <out dir>/<variant>_{judge,law}.log.
set -uo pipefail
AX=${1:?}; REF=${2:?}; V1=${3:?}; AN=${4:?}; AR=${5:?}; NM=${6:?}; MF=${7:?}; RQ=${8:?}; TR=${9:?}; O=${10:?}
H=$(cd "$(dirname "$0")" && pwd); mkdir -p "$O/bin"; A=$(paste -sd, "$AN")
for F in clean break_lift break_abs break_bt3 break_ref break_self; do
  if [ $F = clean ]; then bash "$H/build.sh" "$AX" "$O/bin" > "$O/${F}_build.txt" 2>/dev/null; B=$O/bin/rv2dec
  else bash "$H/build.sh" "$AX" "$O/bin" $F > "$O/${F}_build.txt" 2>/dev/null; B=$O/bin/rv2dec_$F; fi
  /usr/bin/time -v "$B" judge "$REF" "$V1" "$A" "$AR" "$MF" "$NM" 16 > "$O/${F}_judge.log" 2>&1; J=$?
  /usr/bin/time -v taskset -c 3 "$B" law rv2 "$REF" "$V1" "$A" "$AR" "$NM" "$RQ" "$TR" > "$O/${F}_law.log" 2>&1; L=$?
  jt=$(grep '^JUDGE_TOTAL' "$O/${F}_judge.log" | cut -f2,3); lr=$(grep '^LAW_RESULT' "$O/${F}_law.log" | sed 's/.*verdict //')
  lb=$(grep '^WIN' "$O/${F}_law.log" | awk -F'\t' '{split($4,a," "); s+=a[2]} END {print s+0}')
  printf '%s\tjudge exit %d\t%s\tlaw exit %d\tlaw verdict %s\tbad answers %s\n' "$F" $J "${jt:-no total}" $L "${lr:-none}" "$lb" | tee -a "$O/SUMMARY.tsv"
done
