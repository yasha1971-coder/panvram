# panvram 1.0 - release phase 2 (public): runbook

Script: `zenodo/release_phase2.sh` (bash, `set -euo pipefail`). Default mode is the dry run: every local and read-only
check runs for real (hashes, the GET of the Zenodo draft with the token read by Python only and masked, the history
scan, the DOI patch applied in a throwaway worktree `.wk-release-dryrun` + the CI tests), external commands are only
printed (`WOULD RUN [name]: ...`). Nothing is tagged, pushed, published or made public in a dry run.

```
cd ~/panvram
bash zenodo/release_phase2.sh                                   # dry run, all steps (= --dry-run)
bash zenodo/release_phase2.sh --dry-run --step 4                # one step
bash zenodo/release_phase2.sh --execute --step 1 --date 2026-10-09   # for real: one step at a time
```

`--execute`: an external command runs only if every check of its step before it PASSes, and only after its name is
typed back (`type 'push-v1.0.0' to run: ...`). Check states: PASS, FAIL, PEND (waits for an external action of an
earlier step that a dry run does not do). `--date` = the release day (default: today UTC) written to CHANGELOG.md and
CITATION.cff and compared with the Zenodo publication_date. Outputs: `~/outgoing/release_phase2/`.

Order (fixed; Zenodo publish is the last external action):

| step | external actions (names to type) | checks before | checks after |
|---|---|---|---|
| 1 aceapex sources tag | `tag-sources`, `push-sources-tag` | commit d6899513 exists, orphan, tree 04e259bd, research/{refrel,cleanroom,bench2610} == 5b6d5ce / 927a9dd / 927a9dd subtrees == stated hashes, 414 files, author e-mail, secret scan of the tree and message = 0, tag absent or == d6899513 (local, remote), check_source_links against d6899513 = 0 broken | `git ls-remote --tags origin panvram-1.0-sources` == d6899513, check_source_links against the tag = 0 broken |
| 2 panvram 1.0.0 | `release-commit`, `tag-v1.0.0`, `push-v1.0.0` | branch master, clean tree, origin/main ancestor of HEAD, tag v1.0.0 absent (local, remote), CHANGELOG placeholder line; then doi_patch.py -> git apply, release date, no placeholders left, DOI / URL in README, CITATION, .zenodo.json, notebook, CI tests (test_synth.py, test_gate.py, CPU path built in place) | tarball from the tag commit (`git archive --prefix=panvram-1.0.0/ <commit> \| gzip -n`), its commit id, .sha256, rebuild == file, MANIFEST.tsv sha256 3b6941ef..., check_source_links = 0 broken, origin main and v1.0.0^{} == tag commit |
| 3 GitHub | `make-public`, `gh-release` | gh auth, repository private, origin main / tag == tag commit, CI check runs on the tag commit all success, history scan (counts), release notes = CHANGELOG 1.0.0, assets present and of the tag commit, GitHub tarball files == git archive files | `.private` == false, anonymous api.github.com 200, release asset downloaded anonymously sha256 == local, archive URL sha256 recorded (`github-archive-v1.0.0.tar.gz.sha256`), its files == asset files, its MANIFEST.tsv sha256 |
| 4 Zenodo publish (LAST) | `zenodo-metadata` (only if needed), `zenodo-publish` | GET deposition: state unsubmitted, reserved DOI, 48 files, (name, size, md5) == FILES.md table, metadata == dataset.zenodo.json (title, version, license, upload_type, keywords, notes, description as text, creators, related_identifiers), access open, publication_date, the GitHub URL of related_identifiers public (anonymous 200); with --execute also the aceapex tag on origin and the GitHub release | anonymous api/records/23232317 public, record page 200, https://doi.org/10.5281/zenodo.23232317 redirects to Zenodo, anonymous MANIFEST.tsv sha256 == 3b6941ef... |
| 5 Colab final external check | none (generates files; the run is the user's) | the DOI-patched notebook (tag v1.0.0, else HEAD + the step-2 patch) | `colab_final_external_check.txt` compiles, no token / getpass / deposit API, last line `runtime.unassign()`; `quickstart_public_check.ipynb` valid |

## Step 1 by hand (the agent cannot create tags)

```
cd ~/pubrepo
git cat-file -t d68995137cd26689b266964828a7cd79a1588293                  # commit
git rev-parse d68995137cd26689b266964828a7cd79a1588293^{tree}             # 04e259bd56bff97117f2a654605113241e0e3cc7
git tag panvram-1.0-sources d68995137cd26689b266964828a7cd79a1588293
git push origin refs/tags/panvram-1.0-sources
git ls-remote --tags origin panvram-1.0-sources                           # d68995137cd26689b266964828a7cd79a1588293  refs/tags/panvram-1.0-sources
python3 ~/panvram/scripts/check_source_links.py ~/pubrepo panvram-1.0-sources   # links checked: 13 distinct paths; broken: 0
```

(lightweight tag, so `ls-remote` shows the commit itself; the same as `bash zenodo/release_phase2.sh --execute --step 1`.)

## Steps 2-4 for real

```
cd ~/panvram
bash zenodo/release_phase2.sh --execute --step 2 --date <YYYY-MM-DD>   # release-commit, tag-v1.0.0, push-v1.0.0
# wait for CI on the tag commit (gh run list -R yasha1971-coder/panvram)
bash zenodo/release_phase2.sh --execute --step 3 --date <YYYY-MM-DD>   # make-public, gh-release
bash zenodo/release_phase2.sh --execute --step 4 --date <YYYY-MM-DD>   # zenodo-metadata (if needed), zenodo-publish
bash zenodo/release_phase2.sh --dry-run --step 5                       # regenerates the Colab check from the tag
```

The external commands, as the script runs them:

```
git -C ~/panvram add CHANGELOG.md CITATION.cff .zenodo.json README.md notebooks/quickstart_558.ipynb
git -C ~/panvram -c user.name="Yakiv Shavidze" -c user.email=aceapex.contour@outlook.com commit -m "panvram 1.0.0: ..." -- <the same 5 paths>
git -C ~/panvram -c user.name="Yakiv Shavidze" -c user.email=aceapex.contour@outlook.com tag -a v1.0.0 -m "panvram 1.0.0" HEAD
git -C ~/panvram archive --format=tar --prefix=panvram-1.0.0/ <tag commit> | gzip -n > ~/outgoing/release_phase2/panvram-1.0.0.tar.gz
git -C ~/panvram push origin master:main refs/tags/v1.0.0
gh api -X PATCH repos/yasha1971-coder/panvram -F private=false
gh release create v1.0.0 ~/outgoing/release_phase2/panvram-1.0.0.tar.gz ~/outgoing/release_phase2/panvram-1.0.0.tar.gz.sha256 \
   --verify-tag --title "panvram 1.0.0" --notes-file ~/outgoing/release_phase2/release-notes-v1.0.0.md
PUT  https://zenodo.org/api/deposit/depositions/23232317   (only if the metadata differs: description as HTML text, publication_date)
POST https://zenodo.org/api/deposit/depositions/23232317/actions/publish
```

`-F private=false` (typed: JSON boolean); `-f private=false` would send the string "false".

## Decisions and findings of the dry run (before --execute)

1. `notebooks/gate558_colab.ipynb` still holds the placeholder `10.5281/zenodo.XXXXXXX` / `DATASET_RECORD = "XXXXXXX"` and
   expects files the record does not have (`cohort558_q4k.tar`, `manifest_v1.tsv`, `panvram-1.0.0.tar.gz`); doi_patch.py
   does not touch it. Step 2 FAILs on it until it is updated or removed (user's decision).
2. Zenodo sanitises the HTML description: `cohort558_<dataset>.tar.sha256` of dataset.zenodo.json is stored on the draft
   as `cohort558_.tar.sha256`. Step 4 offers `zenodo-metadata` (PUT with `<` `>` escaped as entities, publication_date =
   `--date`); the checks are repeated after it. Without it the published description reads `cohort558_.tar.sha256`.
3. The draft's publication_date is 2026-10-08; if the release day differs, the same `zenodo-metadata` action sets it.
4. GitHub archive vs release asset: the notebook downloads `PANVRAM_SOURCE =
   https://github.com/yasha1971-coder/panvram/archive/refs/tags/v1.0.0.tar.gz` (GitHub's archive), not the release asset
   `panvram-1.0.0.tar.gz` (`git archive | gzip -n`). Same files (checked on 989d415), different bytes, and GitHub does
   not promise stable archive bytes. The notebook checks only MANIFEST.tsv SHA-256 after extraction, so either source
   passes; the archive URL sha256 is recorded in step 3 for the record, and the Colab check logs the sha256 it got
   (`EXPECTED_ARCHIVE_SHA256` can pin it). A notebook cannot carry the SHA-256 of a tarball that contains it.
5. GitHub-Zenodo integration (a software DOI from `.zenodo.json` on each release): not enabled by the script; if wanted,
   switch the repository on in Zenodo -> GitHub before `gh-release` (that is what triggers it).
6. `RELEASE_CHECKLIST.md` becomes public as it is ("prepared locally, not pushed, not public", an open item about
   placeholders and links to branch refrel).
7. The aceapex sources tree has 202 lines with `/home/aeterna` (paths in research logs; not secrets); with the tag they
   are reachable from panvram links. panvram history: 0 hits for gmail addresses, tokens (exact Zenodo and gh token,
   ghp_/github_pat_/gh*_ forms), AKIA keys, private keys; `/home/aeterna` only in RELEASE_CHECKLIST.md and this script.
8. Local commits on master not on origin/main (the release scripts) are pushed together with the release commit.

## Step 5 - the Colab cell of the final external check (after step 4)

Generated by `zenodo/colab_public_check.py` from the DOI-patched `notebooks/quickstart_558.ipynb`: cells 2-4 verbatim
(package with the MANIFEST.tsv stop check, T2T md5, the dataset check chain parts -> PARTS.sha256 -> tar .sha256 ->
MANIFEST.tsv), the gate with form 1,1, dataset q4k, public URLs only (no token), every step OK / FAIL / SKIP in the log
`ACEAPEX-OVH/quickstart_public_<UTC>.txt` on Drive, last line `runtime.unassign()`. Files:
`~/outgoing/release_phase2/colab_final_external_check.txt` (one cell, below) and
`~/outgoing/release_phase2/quickstart_public_check.ipynb` (the same as a notebook). Colab: GPU runtime (A100), paste,
run; the runtime is released at the end.

```python
# panvram 1.0 - final external check from the public record (DOI 10.5281/zenodo.23232317); one Colab cell, GPU runtime, Run.
# Generated by zenodo/colab_public_check.py from quickstart_558.ipynb (HEAD 881e9a1 + dryrun-doi_23232317.patch). No token; Drive only for the log.
# 0. Drive: only the log destination (every input comes from the public URLs)
from google.colab import drive
drive.mount('/content/drive')

# 1. Settings - the public URLs of the DOI-patched quickstart_558.ipynb (HEAD 881e9a1 + dryrun-doi_23232317.patch)
import os, time, json, shutil, subprocess, urllib.request
DATASET_BASE_URL = "https://zenodo.org/records/23232317/files"   # the dataset record, DOI 10.5281/zenodo.23232317
SUFFIX = "?download=1"
PANVRAM_SOURCE = "https://github.com/yasha1971-coder/panvram/archive/refs/tags/v1.0.0.tar.gz"  # the package (tag v1.0.0)
RECORD, DOI = "23232317", "10.5281/zenodo.23232317"
DATASET = "q4k"            # q4k or q16k
FORMS = "1,1"              # final external check: the compact form
W = "/content/panvram_run"
LOG_DIR = "/content/drive/MyDrive/ACEAPEX-OVH"   # the log goes here as quickstart_public_<UTC>.txt
EXPECTED_ARCHIVE_SHA256 = ""   # optional: the sha256 recorded by release_phase2.sh step 3 (github-archive-v1.0.0.tar.gz.sha256); empty = logged only
SKIP_STEPS = set()          # names of steps to skip
os.environ.update(dict(DATASET_BASE_URL=DATASET_BASE_URL, SUFFIX=SUFFIX, PANVRAM_SOURCE=PANVRAM_SOURCE, DATASET=DATASET, FORMS=FORMS, W=W,
                       EXPECTED_ARCHIVE_SHA256=EXPECTED_ARCHIVE_SHA256))

# 1b. Step runner for an unattended run: status per step in the log, FAIL -> later steps SKIP, log copied to Drive after every line
STATUS = {}; os.makedirs(W, exist_ok=True)
LOG = os.path.join(W, "quickstart_public_log.txt"); LOG_PATH = os.path.join(LOG_DIR, "quickstart_public_" + time.strftime("%Y-%m-%dT%H%MZ", time.gmtime()) + ".txt")
def save_log():
    try: os.makedirs(LOG_DIR, exist_ok=True); shutil.copy(LOG, LOG_PATH)
    except Exception as e: print("log copy to Drive failed:", e)
def log(line):
    print(line); open(LOG, "a").write(line + "\n"); save_log()
def failed(): return any(v[0] == "FAIL" for v in STATUS.values())
def _skip(name):
    if name in SKIP_STEPS: STATUS[name] = ("SKIP", "SKIP_STEPS"); log(f"STEP {name}: SKIP (listed in SKIP_STEPS)"); return True
    if failed(): STATUS[name] = ("SKIP", "earlier FAIL"); log(f"STEP {name}: SKIP (an earlier step failed)"); return True
    return False
def run_step(name, script):
    """one bash step (set -euo pipefail); its whole output goes to the log, the tail to the screen"""
    if _skip(name): return
    t0 = time.time(); r = subprocess.run(["bash", "-c", "set -euo pipefail\n" + script], capture_output=True, text=True, env=os.environ)
    with open(LOG, "a") as f: f.write(r.stdout + r.stderr)
    save_log(); print(r.stdout[-2000:] + r.stderr[-2000:])
    if r.returncode == 0: STATUS[name] = ("OK", ""); log(f"STEP {name}: OK ({time.time() - t0:.0f} s)")
    else:
        reason = ((r.stderr.strip().splitlines() or r.stdout.strip().splitlines() or ["no output"])[-1])[:300]
        STATUS[name] = ("FAIL", reason); log(f"STEP {name}: FAIL (exit {r.returncode}, {time.time() - t0:.0f} s): {reason}")
def run_py(name, fn):
    """one Python step; an exception is a FAIL with its message"""
    if _skip(name): return
    t0 = time.time()
    try: fn(); STATUS[name] = ("OK", ""); log(f"STEP {name}: OK ({time.time() - t0:.0f} s)")
    except BaseException as e: STATUS[name] = ("FAIL", f"{type(e).__name__}: {e}"[:300]); log(f"STEP {name}: FAIL ({time.time() - t0:.0f} s): {STATUS[name][1]}")
log(f"panvram public quickstart check: dataset {DATASET}, forms {FORMS}, record {RECORD} (DOI {DOI}), package {PANVRAM_SOURCE}; log {LOG_PATH}")

try:
    # 2a. The public record, anonymously: published, 48 files, the DOI resolves to it
    def _record():
        with urllib.request.urlopen(f"https://zenodo.org/api/records/{RECORD}", timeout=120) as r: rec = json.load(r)
        files = rec.get("files", []); n = len(files.get("entries", files) if isinstance(files, dict) else files)
        log(f"record {RECORD}: status {rec.get('status')}, access {rec.get('access', {}).get('record')}, doi {rec.get('doi') or rec.get('pids', {}).get('doi', {}).get('identifier')}, {n} files")
        assert rec.get("access", {}).get("record", "public") == "public" and n == 48, "record not public or not 48 files"
        with urllib.request.urlopen(f"https://doi.org/{DOI}", timeout=120) as r: final = r.geturl()
        log(f"https://doi.org/{DOI} -> {final}"); assert RECORD in final, "the DOI does not resolve to the record"
    run_py("public-record", _record)
    # 2. panvram from the public GitHub archive (cell 2 of the notebook, verbatim): MANIFEST.tsv SHA-256, CUDA build
    run_step("panvram", r"""
    # 2. panvram: download the package (MANIFEST.tsv SHA-256 == the published value, else stop) and build the CUDA path
    mkdir -p $W && cd $W
    [ -d panvram ] || { if [[ "$PANVRAM_SOURCE" == http* ]]; then curl -fL --retry 3 -o panvram.tar.gz "$PANVRAM_SOURCE"; else cp "$PANVRAM_SOURCE" panvram.tar.gz; fi; mkdir -p panvram && tar -xzf panvram.tar.gz -C panvram --strip-components=1; }
    cd panvram && echo "3b6941efeccb00a7b092573b0482281711f249e62465f771325704d9334d31d0  MANIFEST.tsv" | sha256sum -c - && pip install -q -e . && python -c "import panvram; assert panvram.with_cuda, 'no CUDA path'; print('panvram', panvram.__file__)"
    """)
    # 2b. The package tarball as downloaded: size, sha256, commit id in its header (and == EXPECTED_ARCHIVE_SHA256 if set)
    run_step("provenance", r"""
    cd $W
    echo "panvram.tar.gz: $(stat -c %s panvram.tar.gz) B, sha256 $(sha256sum panvram.tar.gz | cut -d' ' -f1)"
    c=$(gzip -dc panvram.tar.gz | git get-tar-commit-id 2>/dev/null || true); echo "commit in the tar header: ${c:-none}"
    if [ -n "$EXPECTED_ARCHIVE_SHA256" ]; then echo "$EXPECTED_ARCHIVE_SHA256  panvram.tar.gz" | sha256sum -c -
    else echo "EXPECTED_ARCHIVE_SHA256 empty: sha256 logged only"; fi
    """)
    # 3. Reference: T2T-CHM13v2.0 from NCBI, md5 checked (cell 3, verbatim)
    run_step("t2t", r"""
    # 3. Reference: T2T-CHM13v2.0 from NCBI, md5 of the unpacked FASTA checked
    mkdir -p $W/cohort && cd $W
    if [ ! -s t2t.fa ]; then
      curl -fL --retry 3 -o t2t.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
      gunzip -c t2t.fa.gz > t2t.fa && rm -f t2t.fa.gz
    fi
    echo "cd1e52ce400c027ed0b7ab4b9d613f5a  t2t.fa" | md5sum -c -
    ln -sf $W/t2t.fa $W/cohort/t2t.fa
    """)
    # 4. Archives from the public record (cell 4, verbatim): parts -> PARTS.sha256 -> cat -> tar .sha256 -> MANIFEST.tsv
    run_step("archives", r"""
    # 4. Archives: the dataset tar is shipped as parts of 500 MiB. Parts -> SHA-256 of every part (PARTS.sha256) -> cat in
    #    name order -> SHA-256 of the whole tar (.sha256) -> unpacked -> SHA-256 of every archive against MANIFEST.tsv
    cd $W
    t=cohort558_$DATASET.tar
    get() { [ -s "$1" ] && return 0; if [[ "$DATASET_BASE_URL" == http* ]]; then curl -fL --retry 3 -o "$1" "$DATASET_BASE_URL/$1$SUFFIX"; else cp "$DATASET_BASE_URL/$1" "$1"; fi; }
    get PARTS.sha256; get $t.sha256
    parts=$(awk '{print $2}' PARTS.sha256 | grep "^$t.part" | sort)
    for p in $parts; do get $p; done
    grep "  $t.part" PARTS.sha256 | sha256sum -c --quiet   # a mismatch stops the cell (set -e; no && after a check)
    echo "$(echo $parts | wc -w) parts OK"
    [ -s $t ] || cat $parts > $t
    sha256sum -c $t.sha256
    rm -f $parts
    tar -xf $t -C $W/cohort
    rm -f $t
    cd $W/cohort
    col=$(head -1 $W/panvram/MANIFEST.tsv | tr '\t' '\n' | grep -n "^${DATASET}_sha256$" | cut -d: -f1)
    tail -n +2 $W/panvram/MANIFEST.tsv | awk -F'\t' -v c=$col -v d=$DATASET '{print $c "  " $1 "." d ".rr3"}' > want.sha256
    sha256sum -c --quiet want.sha256
    echo "$(wc -l < want.sha256) archives, SHA-256 == MANIFEST.tsv"
    """)
    # 5. On the GPU (cell 5, same command, whole output in the log): sample, fetch, full decode of every assembly == manifest XXH3
    run_step("gate", r"""
    ulimit -c 0
    cd $W/panvram
    for form in $FORMS; do
      PANVRAM_PACKED_REF=${form%,*} PANVRAM_COMPACT_BLOCKS=${form#*,} python scripts/cohort558.py $W/cohort MANIFEST.tsv $DATASET
    done 2>&1
    """)
finally:
    # Always: final status per step and RESULT to the log (on Drive), flush Drive, release the runtime
    try:
        for k, v in STATUS.items(): log(f"FINAL {k}: {v[0]} {v[1]}".rstrip())
        log("RESULT: " + ("PASS" if STATUS and all(v[0] == "OK" for v in STATUS.values()) else "FAIL"))
        save_log()
    finally:
        from google.colab import drive, runtime
        try: drive.flush_and_unmount()
        finally: runtime.unassign()
```
