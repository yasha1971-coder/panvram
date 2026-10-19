#!/usr/bin/env bash
# build.sh [out dir] - builds psearch and crispr from this repository's csrc only (needs g++ >= 11, libzstd-dev, zlib).
# The only project header they include is csrc/pv_v1.h (the public refrel3 v1 reader); nothing from an encoder or from
# a research prototype is used. The build is deterministic on one host: with g++ 11.4.0 on ace-core (AMD EPYC 4344P,
# Ubuntu 22.04.5) it reproduces the binaries whose SHA-256 are recorded in logs/psearch-2026-10-19/REPRODUCIBILITY.md.
set -euo pipefail
H=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$H/../.." && pwd); O=${1:-$H}; mkdir -p "$O"; O=$(cd "$O" && pwd)
for t in psearch crispr; do
  g++ -std=c++17 -O3 -march=native -I"$REPO/csrc" "$H/$t.cpp" -lzstd -lz -lpthread -o "$O/$t"
done
sha256sum "$O/psearch" "$O/crispr"
