# SPDX-License-Identifier: MIT
"""colab_public_check.py <patched quickstart_558.ipynb> <source label> <out cell.txt> <out notebook.ipynb> - the final
external check of the release (release_phase2.sh step 5): the flow of the DOI-patched public notebook, run unattended
on Colab from the public URLs only (Zenodo record, GitHub archive; no token, no Drive input), dataset q4k, form 1,1.
Cells 2-4 of the notebook (package, T2T, the dataset check chain) are taken verbatim; the gate (cell 5) runs the same
command with its whole output in the log. Every step records STEP <name>: OK | FAIL <reason> | SKIP in the log; the log
is copied to Drive (ACEAPEX-OVH/quickstart_public_<UTC>.txt) after every line; the end always flushes Drive and
releases the runtime. Writes a one-cell text (paste into one Colab cell) and the same flow as a notebook."""
import json
import re
import sys

NB, SRC, OUT_TXT, OUT_NB = sys.argv[1:5]
RECID, DOI = "23232317", "10.5281/zenodo.23232317"
cells = json.load(open(NB))["cells"]
src = ["".join(c["source"]) for c in cells]


def bash_body(i, must):
    s = src[i]
    assert s.startswith("%%bash\n"), f"cell {i} is not a %%bash cell"
    lines = s.split("\n")[1:]
    lines = [l for l in lines if l.strip() != "set -euo pipefail"]   # the runner adds set -euo pipefail itself
    body = "\n".join(lines).strip("\n") + "\n"
    for m in must:
        assert m in body, f"cell {i}: expected text not found: {m}"
    assert '"""' not in body
    return body


settings = src[1]
base = re.search(r'^DATASET_BASE_URL = "([^"]+)"', settings, re.M).group(1)
pv_src = re.search(r'^PANVRAM_SOURCE = "([^"]+)"', settings, re.M).group(1)
assert base == f"https://zenodo.org/records/{RECID}/files", f"DATASET_BASE_URL not the public record: {base}"
assert pv_src == "https://github.com/yasha1971-coder/panvram/archive/refs/tags/v1.0.0.tar.gz", pv_src
assert 'SUFFIX = "?download=1" if DATASET_BASE_URL.startswith("http") else ""' in settings
PANVRAM_SH = bash_body(2, ["3b6941efeccb00a7b092573b0482281711f249e62465f771325704d9334d31d0  MANIFEST.tsv", "--strip-components=1"])
T2T_SH = bash_body(3, ["cd1e52ce400c027ed0b7ab4b9d613f5a  t2t.fa"])
ARCHIVES_SH = bash_body(4, ['grep "  $t.part" PARTS.sha256 | sha256sum -c --quiet', "sha256sum -c $t.sha256", "sha256sum -c --quiet want.sha256"])
assert "python scripts/cohort558.py $W/cohort MANIFEST.tsv $DATASET" in src[5]
GATE_SH = """ulimit -c 0
cd $W/panvram
for form in $FORMS; do
  PANVRAM_PACKED_REF=${form%,*} PANVRAM_COMPACT_BLOCKS=${form#*,} python scripts/cohort558.py $W/cohort MANIFEST.tsv $DATASET
done 2>&1
"""
PROV_SH = """cd $W
echo "panvram.tar.gz: $(stat -c %s panvram.tar.gz) B, sha256 $(sha256sum panvram.tar.gz | cut -d' ' -f1)"
c=$(gzip -dc panvram.tar.gz | git get-tar-commit-id 2>/dev/null || true); echo "commit in the tar header: ${c:-none}"
if [ -n "$EXPECTED_ARCHIVE_SHA256" ]; then echo "$EXPECTED_ARCHIVE_SHA256  panvram.tar.gz" | sha256sum -c -
else echo "EXPECTED_ARCHIVE_SHA256 empty: sha256 logged only"; fi
"""

MOUNT = """# 0. Drive: only the log destination (every input comes from the public URLs)
from google.colab import drive
drive.mount('/content/drive')
"""
SETTINGS = f'''# 1. Settings - the public URLs of the DOI-patched quickstart_558.ipynb ({SRC})
import os, time, json, shutil, subprocess, urllib.request
DATASET_BASE_URL = "{base}"   # the dataset record, DOI {DOI}
SUFFIX = "?download=1"
PANVRAM_SOURCE = "{pv_src}"  # the package (tag v1.0.0)
RECORD, DOI = "{RECID}", "{DOI}"
DATASET = "q4k"            # q4k or q16k
FORMS = "1,1"              # final external check: the compact form
W = "/content/panvram_run"
LOG_DIR = "/content/drive/MyDrive/ACEAPEX-OVH"   # the log goes here as quickstart_public_<UTC>.txt
EXPECTED_ARCHIVE_SHA256 = ""   # optional: the sha256 recorded by release_phase2.sh step 3 (github-archive-v1.0.0.tar.gz.sha256); empty = logged only
SKIP_STEPS = set()          # names of steps to skip
os.environ.update(dict(DATASET_BASE_URL=DATASET_BASE_URL, SUFFIX=SUFFIX, PANVRAM_SOURCE=PANVRAM_SOURCE, DATASET=DATASET, FORMS=FORMS, W=W,
                       EXPECTED_ARCHIVE_SHA256=EXPECTED_ARCHIVE_SHA256))

# 1b. Step runner for an unattended run: status per step in the log, FAIL -> later steps SKIP, log copied to Drive after every line
STATUS = {{}}; os.makedirs(W, exist_ok=True)
LOG = os.path.join(W, "quickstart_public_log.txt"); LOG_PATH = os.path.join(LOG_DIR, "quickstart_public_" + time.strftime("%Y-%m-%dT%H%MZ", time.gmtime()) + ".txt")
def save_log():
    try: os.makedirs(LOG_DIR, exist_ok=True); shutil.copy(LOG, LOG_PATH)
    except Exception as e: print("log copy to Drive failed:", e)
def log(line):
    print(line); open(LOG, "a").write(line + "\\n"); save_log()
def failed(): return any(v[0] == "FAIL" for v in STATUS.values())
def _skip(name):
    if name in SKIP_STEPS: STATUS[name] = ("SKIP", "SKIP_STEPS"); log(f"STEP {{name}}: SKIP (listed in SKIP_STEPS)"); return True
    if failed(): STATUS[name] = ("SKIP", "earlier FAIL"); log(f"STEP {{name}}: SKIP (an earlier step failed)"); return True
    return False
def run_step(name, script):
    """one bash step (set -euo pipefail); its whole output goes to the log, the tail to the screen"""
    if _skip(name): return
    t0 = time.time(); r = subprocess.run(["bash", "-c", "set -euo pipefail\\n" + script], capture_output=True, text=True, env=os.environ)
    with open(LOG, "a") as f: f.write(r.stdout + r.stderr)
    save_log(); print(r.stdout[-2000:] + r.stderr[-2000:])
    if r.returncode == 0: STATUS[name] = ("OK", ""); log(f"STEP {{name}}: OK ({{time.time() - t0:.0f}} s)")
    else:
        reason = ((r.stderr.strip().splitlines() or r.stdout.strip().splitlines() or ["no output"])[-1])[:300]
        STATUS[name] = ("FAIL", reason); log(f"STEP {{name}}: FAIL (exit {{r.returncode}}, {{time.time() - t0:.0f}} s): {{reason}}")
def run_py(name, fn):
    """one Python step; an exception is a FAIL with its message"""
    if _skip(name): return
    t0 = time.time()
    try: fn(); STATUS[name] = ("OK", ""); log(f"STEP {{name}}: OK ({{time.time() - t0:.0f}} s)")
    except BaseException as e: STATUS[name] = ("FAIL", f"{{type(e).__name__}}: {{e}}"[:300]); log(f"STEP {{name}}: FAIL ({{time.time() - t0:.0f}} s): {{STATUS[name][1]}}")
log(f"panvram public quickstart check: dataset {{DATASET}}, forms {{FORMS}}, record {{RECORD}} (DOI {{DOI}}), package {{PANVRAM_SOURCE}}; log {{LOG_PATH}}")
'''
RECORD_PY = '''# 2a. The public record, anonymously: published, 48 files, the DOI resolves to it
def _record():
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{RECORD}", timeout=120) as r: rec = json.load(r)
    files = rec.get("files", []); n = len(files.get("entries", files) if isinstance(files, dict) else files)
    log(f"record {RECORD}: status {rec.get('status')}, access {rec.get('access', {}).get('record')}, doi {rec.get('doi') or rec.get('pids', {}).get('doi', {}).get('identifier')}, {n} files")
    assert rec.get("access", {}).get("record", "public") == "public" and n == 48, "record not public or not 48 files"
    with urllib.request.urlopen(f"https://doi.org/{DOI}", timeout=120) as r: final = r.geturl()
    log(f"https://doi.org/{DOI} -> {final}"); assert RECORD in final, "the DOI does not resolve to the record"
run_py("public-record", _record)
'''
FINAL = '''# Always: final status per step and RESULT to the log (on Drive), flush Drive, release the runtime
try:
    for k, v in STATUS.items(): log(f"FINAL {k}: {v[0]} {v[1]}".rstrip())
    log("RESULT: " + ("PASS" if STATUS and all(v[0] == "OK" for v in STATUS.values()) else "FAIL"))
    save_log()
finally:
    from google.colab import drive, runtime
    try: drive.flush_and_unmount()
    finally: runtime.unassign()
'''


def step(name, body, title):
    return f'# {title}\nrun_step("{name}", r"""\n{body}""")\n'


STEPS = [
    ("panvram", PANVRAM_SH, "2. panvram from the public GitHub archive (cell 2 of the notebook, verbatim): MANIFEST.tsv SHA-256, CUDA build"),
    ("provenance", PROV_SH, "2b. The package tarball as downloaded: size, sha256, commit id in its header (and == EXPECTED_ARCHIVE_SHA256 if set)"),
    ("t2t", T2T_SH, "3. Reference: T2T-CHM13v2.0 from NCBI, md5 checked (cell 3, verbatim)"),
    ("archives", ARCHIVES_SH, "4. Archives from the public record (cell 4, verbatim): parts -> PARTS.sha256 -> cat -> tar .sha256 -> MANIFEST.tsv"),
    ("gate", GATE_SH, "5. On the GPU (cell 5, same command, whole output in the log): sample, fetch, full decode of every assembly == manifest XXH3"),
]
MD = f"""# panvram 1.0 - final external check from the public record (DOI {DOI})

The public quickstart flow (`notebooks/quickstart_558.ipynb` after the DOI patch; {SRC}) run unattended from the
public URLs only: the dataset from the Zenodo record {RECID} (anonymous, `?download=1`), the package from
`{pv_src}`. No token, nothing read from Drive; Drive only receives the log `ACEAPEX-OVH/quickstart_public_<UTC>.txt`.
Dataset q4k, resident form 1,1. Every step records `STEP <name>: OK | FAIL <reason> | SKIP`; a FAIL makes the later
steps SKIP; the last cell always flushes Drive and releases the runtime. Needs a CUDA GPU of sm_80 or newer and ~40 GB
of disk (q4k).
"""

nb_cells = [("markdown", MD), ("code", MOUNT), ("code", SETTINGS), ("code", RECORD_PY)] + \
           [("code", step(n, b, t)) for n, b, t in STEPS] + [("code", FINAL)]
nb = {"cells": [{"cell_type": t, "metadata": {}, "source": s.splitlines(True), **({"execution_count": None, "outputs": []} if t == "code" else {})}
                for t, s in nb_cells],
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python"},
                   "accelerator": "GPU", "colab": {"provenance": []}},
      "nbformat": 4, "nbformat_minor": 0}
open(OUT_NB, "w").write(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")

indent = lambda s: "".join(("    " + l if l.strip() else l) for l in s.splitlines(True))
txt = (f"# panvram 1.0 - final external check from the public record (DOI {DOI}); one Colab cell, GPU runtime, Run.\n"
       f"# Generated by zenodo/colab_public_check.py from quickstart_558.ipynb ({SRC}). No token; Drive only for the log.\n"
       + MOUNT + "\n" + SETTINGS + "\ntry:\n" + indent(RECORD_PY) + "".join(indent(step(n, b, t)) for n, b, t in STEPS)
       + "finally:\n" + indent(FINAL))
assert txt.rstrip().splitlines()[-1].strip().endswith("runtime.unassign()")
open(OUT_TXT, "w").write(txt.rstrip() + "\n")
print("written:", OUT_TXT, OUT_NB)
