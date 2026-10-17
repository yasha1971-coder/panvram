#!/usr/bin/env bash
# build.sh <aceapex checkout at e902934> [out dir] - builds rr3bench exactly as used for M6 (2026-10-07) and the
# mech558 counters (2026-10-09). The refrel sources here are aceapex research/refrel @ 5b6d5ce and rr3bench.cpp /
# bench_common.h are aceapex research/agc_vs_refrel3 @ 2cc43cc, unchanged. The aceapex library sources (src/) come from
# the public aceapex repository at e902934 (same src tree as 5b6d5ce): https://github.com/yasha1971-coder/aceapex
# With g++ 11.4.0 on ace-core (Ubuntu 22.04.5, AMD EPYC 4344P) this produces SHA-256
# f2f84fe9f396d3a67ede44d81eafb4959330baa7069584953c47d194a6e21673, the binary recorded in the M6 and mech558 inputs.
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd); AX=$(cd "${1:?aceapex checkout at e902934}" && pwd); O=${2:-$H}; mkdir -p "$O"; O=$(cd "$O" && pwd)
[ "$(grep -c '^int main(int argc, char\*\* argv) {$' "$H/refrel3v1.cpp")" = 1 ]
sed 's/^int main(int argc, char\*\* argv) {$/static int refrel3v1_main(int argc, char** argv) {/' "$H/refrel3v1.cpp" > "$O/refrel3v1_nomain.cpp"
g++ -std=c++17 -O3 -march=native -funroll-loops -I"$O" -I"$AX/src" -I"$H" "$H/rr3bench.cpp" "$AX/src/aceapex_api.cpp" \
    -lzstd -lpthread /usr/lib/x86_64-linux-gnu/libcrypto.so.3 -o "$O/rr3bench"
sha256sum "$O/rr3bench"
