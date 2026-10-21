#!/usr/bin/env bash
# build.sh <aceapex checkout at e902934> [out dir] [fault] - builds rv2dec, the decoder of the refrel v2 prototype
# (format unstable). rv2dec.cpp and rv2.h are here; the v1 container code (refrel3v1.cpp and its headers, unchanged) and
# bench_common.h are taken from tools/rr3bench, the aceapex library sources (src/) from the public aceapex repository
# at e902934: https://github.com/yasha1971-coder/aceapex . Same compiler flags as the prototype runs.
# With g++ 11.4.0 on ace-core (Ubuntu 22.04.5, AMD EPYC 4344P) this produces SHA-256
#   6e1f91b7a1fa1846e6d902d0ad33f052e948591d5fcaf22af3986165e011ad97
# [fault]: one of break_lift, break_abs, break_bt3, break_ref, break_self - builds a deliberately wrong decoder (negative
# control, see logs/refrel-v2-2026-10-21/negative_controls/); without it no fault code is compiled.
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd); RB=$(cd "$H/../../tools/rr3bench" && pwd); AX=$(cd "${1:?aceapex checkout at e902934}" && pwd)
O=${2:-$H}; mkdir -p "$O"; O=$(cd "$O" && pwd); F=${3:-}; D=
case $F in "") ;; break_lift|break_abs|break_bt3|break_ref|break_self) D=-DRV2_FAULT_$(echo "${F#break_}" | tr a-z A-Z);; *) echo "unknown fault $F" >&2; exit 1;; esac
[ "$(grep -c '^int main(int argc, char\*\* argv) {$' "$RB/refrel3v1.cpp")" = 1 ]
sed 's/^int main(int argc, char\*\* argv) {$/static int refrel3v1_main(int argc, char** argv) {/' "$RB/refrel3v1.cpp" > "$O/refrel3v1_nomain.cpp"
g++ -std=c++17 -O3 -march=native -funroll-loops $D -I"$O" -I"$AX/src" -I"$RB" -I"$H" "$H/rv2dec.cpp" "$AX/src/aceapex_api.cpp" \
    -lzstd -lpthread /usr/lib/x86_64-linux-gnu/libcrypto.so.3 -o "$O/rv2dec${F:+_$F}"
sha256sum "$O/rv2dec${F:+_$F}"
