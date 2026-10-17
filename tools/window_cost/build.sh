#!/usr/bin/env bash
# build.sh [out dir] - builds window_cost from this repository's csrc only (needs g++ >= 11, libzstd-dev).
# Flags are those of the v1 cohort encoder (FORMAT.md section 7); the measurement does not depend on them, the
# binary hash does: record `sha256sum window_cost` with every run.
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$H/../.." && pwd); O=${1:-$H}
g++ -std=c++17 -O3 -march=native -funroll-loops -I"$REPO/csrc" "$H/window_cost.cpp" -lzstd -lpthread -o "$O/window_cost"
sha256sum "$O/window_cost"
