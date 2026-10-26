"""Installed-wheel CPU check, outside the source import path; emits hash evidence."""
import hashlib
import json
from pathlib import Path
import random
import sys

import torch
import panvram
from fasta_cut import FastaIndex


def main():
    root = Path(sys.argv[1]).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    assert manifest["kind"] == "synthetic-v1-five-assemblies"
    assert len(manifest["assemblies"]) == 5
    assert len({a["name"] for a in manifest["assemblies"]}) == 5
    for name, expected in manifest["sha256"].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected, name
    assert not panvram.with_cuda
    assert "site-packages" in str(Path(panvram.__file__).resolve()), panvram.__file__
    c = panvram.Cohort.open(str(root / "q4k"), reference=str(root / manifest["reference"]),
                           dataset="q4k", device="cpu", threads=1)
    assert c.names == sorted(a["name"] for a in manifest["assemblies"])
    results = []
    rng = random.Random(731)
    for a in manifest["assemblies"]:
        source = (root / a["fasta"]).read_bytes()
        decoded = c.fasta(a["name"])
        assert decoded == source, a["name"]
        index = FastaIndex(source)
        windows = 0
        for contig, _, length, offset, _ in c.contigs(a["name"]):
            for width in (1, 13, 4096, 8192):
                if width > length:
                    continue
                for start in (0, length - width, rng.randrange(length - width + 1)):
                    truth = index.cut(contig, start, width)
                    got = c.windows(torch.tensor([c.names.index(a["name"])]),
                                    torch.tensor([offset + start]), width)
                    assert got.dtype == torch.uint8 and got.device.type == "cpu"
                    assert got.numpy().tobytes() == truth
                    windows += 1
        results.append({"assembly": a["name"], "fasta_sha256": hashlib.sha256(source).hexdigest(),
                        "decoded_sha256": hashlib.sha256(decoded).hexdigest(), "windows": windows})
    assert len({x["fasta_sha256"] for x in results}) == 5
    print(json.dumps({"status": "PASS", "kind": manifest["kind"], "assemblies": results,
                      "torch": torch.__version__, "panvram": panvram.__version__,
                      "module": panvram.__file__}, sort_keys=True))


if __name__ == "__main__":
    main()
