#!/usr/bin/env bash
# make_dataset_tars.sh <v1 archive dir> <out dir> - the two Zenodo tars, deterministic: cohort558_<ds>.tar holds the 558
# <name>.<ds>.rr3 of MANIFEST.tsv and SHA256SUMS.<ds> as SHA256SUMS, at the tar root (the layout run_colab_cohort558.sh
# and quickstart_558.ipynb unpack); GNU tar --sort=name --owner=0 --group=0 --numeric-owner --mtime='2026-10-04 00:00Z'
# --format=gnu, files read through a staging directory of symlinks (-h); <tar>.sha256 next to each. Check: every member
# read back from the tar (streamed, nothing written) and its SHA-256 == MANIFEST.tsv; member list == manifest + SHA256SUMS.
set -euo pipefail
ulimit -c 0
H=$(cd "$(dirname "$0")/.." && pwd); V1=$(cd "${1:?}" && pwd); O=${2:?}; mkdir -p "$O"; O=$(cd "$O" && pwd)
for ds in q4k q16k; do
  S=$O/stage_$ds; mkdir -p "$S"
  col=$(head -1 "$H/MANIFEST.tsv" | tr '\t' '\n' | grep -n "^${ds}_sha256$" | cut -d: -f1)
  tail -n +2 "$H/MANIFEST.tsv" | cut -f1 | while read -r n; do ln -sfn "$V1/$n.$ds.rr3" "$S/$n.$ds.rr3"; done
  cp "$H/SHA256SUMS.$ds" "$S/SHA256SUMS"
  [ "$(ls "$S" | wc -l)" = 559 ]
  T=$O/cohort558_$ds.tar
  [ -e "$T" ] && { echo "exists: $T"; exit 1; }
  (cd "$S" && LC_ALL=C tar --sort=name --owner=0 --group=0 --numeric-owner --mtime='2026-10-04 00:00Z' --format=gnu -h -cf "$T" .)
  (cd "$O" && sha256sum "$(basename "$T")" > "$(basename "$T").sha256")
  python3 - "$T" "$H/MANIFEST.tsv" "$ds" "$col" <<'PY'
import hashlib, sys, tarfile
t, man, ds, col = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
want = {l.split("\t")[0] + "." + ds + ".rr3": l.rstrip("\n").split("\t")[col - 1] for l in list(open(man))[1:]}
seen, ok = set(), 0
with tarfile.open(t) as tf:
    for m in tf:
        if not m.isfile(): continue
        name = m.name[2:] if m.name.startswith("./") else m.name
        h = hashlib.sha256(); f = tf.extractfile(m)
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
        seen.add(name)
        if name in want: ok += h.hexdigest() == want[name]
print(f"{t}: {len(seen)} members, {ok}/{len(want)} archives SHA-256 == MANIFEST.tsv, SHA256SUMS present {'SHA256SUMS' in seen}, extra {sorted(seen - set(want) - {'SHA256SUMS'})}")
sys.exit(0 if ok == len(want) and seen == set(want) | {"SHA256SUMS"} else 1)
PY
  stat -c '%s %n' "$T"; cat "$T.sha256"
done
