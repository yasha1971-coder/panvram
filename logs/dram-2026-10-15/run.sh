#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
task_out=${OUT_DIR:-"$PWD/DRAM_LAT_results_$(date -u +%Y%m%dT%H%M%SZ)_$$"}
mkdir "$task_out" # refuse an existing directory; preserve previous evidence
task_out=$(cd -- "$task_out" && pwd)
task_ok=0
trap 'if [[ $task_ok != 1 ]]; then printf "FAIL\n" > "$task_out/STATUS"; fi' EXIT
for task_cmd in gcc awk date sleep getconf lscpu uname sha256sum; do
    command -v "$task_cmd" >/dev/null || { echo "missing: $task_cmd" >&2; exit 1; }
done
uname -a > "$task_out/uname.txt"
lscpu > "$task_out/lscpu.txt"
gcc --version > "$task_out/gcc.txt"
sha256sum "$task_dir/dram_lat.c" "$task_dir/run.sh" > "$task_out/source_sha256.txt"
task_load=$(awk '{print $1}' /proc/loadavg)
printf 'load1=%s\n' "$task_load" > "$task_out/gate.txt"
awk -v l1="$task_load" 'BEGIN {exit !(l1+0 < 0.5)}' || { echo 'STOP: load1 >= 0.5' >&2; exit 1; }
# /proc stat fields after the final ") ": state=1, ppid=2, utime=12, stime=13, starttime=20.
# Use two tick samples, not ps %CPU (which is a lifetime average).
snapshot() {
    awk 'FNR==1 {line=$0;sub(/^.*\) /,"",line);n=split(line,a," ");
         if(n<20){exit 1}p=FILENAME;sub("/proc/","",p);sub("/stat","",p);print p,a[2],a[12]+a[13],a[20]}' \
         /proc/[0-9]*/stat # fail closed, including unreadable or disappearing processes
}
snapshot > "$task_out/process_before.tsv"
[[ -s "$task_out/process_before.tsv" ]] || { echo 'STOP: empty process sample' >&2; exit 1; }
task_start=$(awk '{print $1}' /proc/uptime)
sleep 1.1
snapshot > "$task_out/process_after.tsv"
[[ -s "$task_out/process_after.tsv" ]] || { echo 'STOP: empty process sample' >&2; exit 1; }
task_end=$(awk '{print $1}' /proc/uptime)
task_hz=$(getconf CLK_TCK)
awk -v start="$task_start" -v end="$task_end" -v hz="$task_hz" -v own="$$" '
    NR==FNR {ticks[$1]=$3;birth[$1]=$4;next}
    {pp[$1]=$2;now[$1]=$3;born[$1]=$4}
    END {
      seconds=end-start;if(seconds<1 || hz<=0){print "invalid CPU sample";exit 1}
      printf "sampling_seconds=%.6f\n",seconds
      bad=0
      for(pid in now){
        p=pid;ours=0;for(k=0;k<1024 && p>0;k++){if(p==own){ours=1;break}p=pp[p]}
        delta=(pid in ticks && birth[pid]==born[pid])?now[pid]-ticks[pid]:now[pid]
        if(delta<0){bad=1;continue}
        cpu=100*delta/hz/seconds
        if(!ours && cpu>=10){printf "foreign_pid=%s cpu_percent=%.6f\n",pid,cpu;bad=1}
      }
      exit bad
    }' "$task_out/process_before.tsv" "$task_out/process_after.tsv" >> "$task_out/gate.txt" \
    || { echo 'STOP: foreign process CPU >= 10% or invalid sample' >&2; exit 1; }
task_load=$(awk '{print $1}' /proc/loadavg)
printf 'load1_after_cpu_sample=%s\n' "$task_load" >> "$task_out/gate.txt"
awk -v l1="$task_load" 'BEGIN {exit !(l1+0 < 0.5)}' || { echo 'STOP: load1 >= 0.5' >&2; exit 1; }
task_available=$(awk '/^MemAvailable:/ {print $2*1024}' /proc/meminfo)
awk -v available="$task_available" 'BEGIN {exit !(available >= 5368709120)}' \
    || { echo 'STOP: need at least 5 GiB MemAvailable' >&2; exit 1; }
printf 'MemAvailable_bytes=%s\n' "$task_available" >> "$task_out/gate.txt"
gcc -O2 -std=c11 -Wall -Wextra -Werror "$task_dir/dram_lat.c" -o "$task_out/dram_lat"
"$task_out/dram_lat" --self-test > "$task_out/self_test.txt"
sha256sum "$task_out/dram_lat" > "$task_out/binary_sha256.txt"
task_args=()
if [[ -n ${CPU:-} ]]; then task_args+=(--cpu "$CPU"); fi
printf 'gcc -O2 -std=c11 -Wall -Wextra -Werror\nCPU=%s\nseed=20261012\n' "${CPU:-first allowed CPU}" > "$task_out/command.txt"
"$task_out/dram_lat" "${task_args[@]}" > "$task_out/results.partial.tsv" 2> "$task_out/run.log"
mv "$task_out/results.partial.tsv" "$task_out/results.tsv"
printf 'PASS\n' > "$task_out/STATUS"
task_ok=1
printf 'PASS: %s/results.tsv\n' "$task_out"
