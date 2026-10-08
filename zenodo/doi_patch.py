# SPDX-License-Identifier: MIT
"""doi_patch.py <reserved DOI> <out.patch> - writes (does not apply) a git patch that puts the reserved dataset DOI and the
repository URL into CITATION.cff, .zenodo.json, README.md and notebooks/quickstart_558.ipynb. Run from the repository
root on a clean tree; the tree is restored after the diff is taken."""
import json, re, subprocess, sys
from pathlib import Path

DOI, OUT = sys.argv[1], sys.argv[2]
REPO_URL = "https://github.com/yasha1971-coder/panvram"
RECORD = DOI.rsplit(".", 1)[1]
if not re.fullmatch(r"10\.5281/zenodo\.\d+", DOI): sys.exit("expected a Zenodo DOI 10.5281/zenodo.<id>")
if subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True).stdout.strip(): sys.exit("tree not clean")
edits = {}
p = Path("CITATION.cff"); s = p.read_text()
s = s.replace('license: MIT\n', f'license: MIT\nrepository-code: "{REPO_URL}"\nidentifiers:\n  - type: doi\n    value: "{DOI}"\n    description: "The HPRC cohort dataset (refrel3 v1 archives, q4k and q16k) this software reads"\n')
edits[p] = s
p = Path(".zenodo.json"); z = json.loads(p.read_text())
for r in z.get("related_identifiers", []):
    if r.get("identifier", "").startswith("10.5281/zenodo."): r["identifier"] = DOI
z["related_identifiers"] = [r for r in z["related_identifiers"] if not r["identifier"].startswith(REPO_URL)] + [{"identifier": REPO_URL, "relation": "isSupplementTo", "scheme": "url", "resource_type": "software"}]
edits[p] = json.dumps(z, indent=2, ensure_ascii=False) + "\n"
p = Path("README.md"); s = p.read_text()
s = s.replace("Version 1.0.0. Numbers below come only from the logs named next to them.",
              f"Version 1.0.0. Repository: {REPO_URL}. Dataset (the 558 refrel3 v1 archives, q4k and q16k): DOI {DOI}.\nNumbers below come only from the logs named next to them.")
edits[p] = s
p = Path("notebooks/quickstart_558.ipynb"); nb = json.loads(p.read_text())
for c in nb["cells"]:
    src = "".join(c["source"])
    src = src.replace('DATASET_BASE_URL = "https://zenodo.org/records/REPLACE-RECORD-ID/files"   # placeholder until the DOI exists; a local/Drive directory path also works',
                      f'DATASET_BASE_URL = "https://zenodo.org/records/{RECORD}/files"   # the dataset record, DOI {DOI}; a local/Drive directory path also works')
    src = src.replace("**Before running:** set `DATASET_BASE_URL` (the dataset record - placeholder until the DOI exists - or a local / Drive directory holding the same files) and `PANVRAM_SOURCE` (the package release tag, or a local tarball) - the package source.",
                      f"**Before running:** nothing - `DATASET_BASE_URL` is the dataset record (DOI {DOI}) and `PANVRAM_SOURCE` the package release (tag v1.0.0); a local / Drive directory or tarball also works.")
    c["source"] = src.splitlines(True)
edits[p] = json.dumps(nb, indent=1, ensure_ascii=False)
originals = {p: p.read_bytes() for p in edits}
try:
    for p, s in edits.items(): p.write_text(s)
    diff = subprocess.run(["git", "diff"], capture_output=True, text=True, check=True).stdout
finally:
    for p, b in originals.items(): p.write_bytes(b)
Path(OUT).write_text(diff)
print("patch written:", OUT, len(diff.splitlines()), "lines; tree restored:", not subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True).stdout.strip())
