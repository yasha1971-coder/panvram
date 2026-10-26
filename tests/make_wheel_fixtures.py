"""Generate five distinct synthetic v1 assemblies; never labelled HPRC subset5."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import synth


def main():
    encoder, target = sys.argv[1:]
    root = Path(target)
    root.mkdir(parents=True, exist_ok=False)
    (root / "q4k").mkdir()
    ref = root / "synth_ref.fa"
    records = synth.reference()
    ref.write_bytes(synth._fasta(records, [80] * len(records)))
    assemblies = []
    for seed in range(1, 6):
        name = f"wheel_synth{seed}"
        records = synth.assembly(seed)
        fasta = root / (name + ".fa")
        fasta.write_bytes(synth._fasta(records, [60 if i % 2 == 0 else 70 for i in range(len(records))]))
        archive = root / "q4k" / (name + ".q4k.rr3")
        subprocess.run([encoder, "encode", str(ref), "4096", "1", str(fasta), str(archive)], check=True)
        assemblies.append({"name": name, "seed": seed, "fasta": fasta.name,
                           "archive": archive.relative_to(root).as_posix()})
    hashes = {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(root.rglob("*")) if p.is_file()}
    manifest = {"kind": "synthetic-v1-five-assemblies", "reference": ref.name,
                "assemblies": assemblies, "sha256": hashes}
    (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
