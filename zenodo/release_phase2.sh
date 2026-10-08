#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# release_phase2.sh - panvram 1.0 public release ("phase 2"), five steps in a fixed order; runbook: zenodo/RELEASE_PHASE2.md
#   1 aceapex: tag panvram-1.0-sources on the orphan sources commit + push of that tag
#   2 panvram: DOI patch, release date, tests, commit, annotated tag v1.0.0, tarball, push of main + tag
#   3 GitHub: repository public, release v1.0.0 with the tarball and its .sha256
#   4 Zenodo: publish the dataset draft 23232317 (LAST external action)
#   5 the Colab cell / notebook of the final external check from the public URLs
# Every step: checks before (PASS / FAIL / PEND = depends on an external action of an earlier step not done yet),
# the command, checks after.
# Usage: release_phase2.sh [--dry-run | --execute] [--step N] [--date YYYY-MM-DD]
#   --dry-run (default): every local and read-only check runs for real (hashes, GET of the Zenodo draft, history scan,
#     the DOI patch applied in a throwaway worktree .wk-release-dryrun + tests); external commands are only printed;
#     no tag is created, nothing is pushed or published.
#   --execute: external commands run only when all checks before them PASS and the step name is typed back.
# The Zenodo token (~/.zenodo_token) is read by Python only, never put on a command line, never printed.
set -euo pipefail
umask 022

MODE=dry; ONLY=""; DATE=$(date -u +%F)
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) MODE=dry ;;
    --execute) MODE=exec ;;
    --step) ONLY=$2; shift ;;
    --date) DATE=$2; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac; shift
done
[[ "$DATE" =~ ^20[0-9]{2}-[01][0-9]-[0-3][0-9]$ ]] || { echo "--date: YYYY-MM-DD expected" >&2; exit 2; }
[[ -z "$ONLY" || "$ONLY" =~ ^[1-5]$ ]] || { echo "--step: 1..5" >&2; exit 2; }

# ------------------------------------------------------------------ constants
ACEAPEX=${ACEAPEX:-$HOME/pubrepo}
PANVRAM=${PANVRAM:-$HOME/panvram}
OUT=${OUT:-$HOME/outgoing/release_phase2}
PY=${PY:-$PANVRAM/.venv/bin/python}          # interpreter with torch / numpy / pytest (the CI test set)
TOKFILE=$HOME/.zenodo_token
SRC_COMMIT=d68995137cd26689b266964828a7cd79a1588293
SRC_TREE=04e259bd56bff97117f2a654605113241e0e3cc7
SRC_TAG=panvram-1.0-sources
declare -A SUB_REF=([refrel]=5b6d5ce [cleanroom]=927a9dd [bench2610]=927a9dd)
declare -A SUB_TREE=([refrel]=c96b8f0167c54ae6645aa18d7c30c697a18178a0 [cleanroom]=1208ed2af9b373779cefc05edb94cd5f7c5554d7 [bench2610]=3649681f1517406bf752e8bbdc8e727887a2143f)
DOI=10.5281/zenodo.23232317; RECID=23232317
GH=yasha1971-coder/panvram; GH_URL=https://github.com/$GH
TAG=v1.0.0; PREFIX=panvram-1.0.0; TARBALL=$OUT/$PREFIX.tar.gz
ARCHIVE_URL=$GH_URL/archive/refs/tags/$TAG.tar.gz
ASSET_URL=$GH_URL/releases/download/$TAG/$PREFIX.tar.gz
MANIFEST_SHA=3b6941efeccb00a7b092573b0482281711f249e62465f771325704d9334d31d0
WT=$PANVRAM/.wk-release-dryrun
RELEASE_PATHS=(CHANGELOG.md CITATION.cff .zenodo.json README.md notebooks/quickstart_558.ipynb)
GIT_ID=(-c user.name="Yakiv Shavidze" -c user.email="aceapex.contour@outlook.com")
mkdir -p "$OUT"
TMP=$(mktemp -d); WT_MADE=0
cleanup() { if [ "$WT_MADE" = 1 ] && [ -d "$WT" ]; then git -C "$PANVRAM" worktree remove --force "$WT" >/dev/null 2>&1 || true; git -C "$PANVRAM" worktree prune || true; fi; rm -rf "$TMP"; }
trap cleanup EXIT
mask() { sed -E 's/[A-Za-z0-9]{60}/<token>/g'; }

# ------------------------------------------------------------------ check helpers (never abort the run)
P=0; F=0; N=0; STEP_F=0; SUMMARY=()
ok()   { echo "  PASS  $*"; P=$((P + 1)); }
bad()  { echo "  FAIL  $*"; F=$((F + 1)); STEP_F=$((STEP_F + 1)); }
pend() { echo "  PEND  $*"; N=$((N + 1)); }
note() { echo "  NOTE  $*"; }
chk()  { local d=$1; shift; if "$@" >/dev/null 2>&1; then ok "$d"; else bad "$d"; fi; }
eq()   { if [ "$2" = "$3" ]; then ok "$1: $2"; else bad "$1: '$2' != '$3'"; fi; }
# eq_or_pend <desc> <have> <want>: in a dry run a mismatch is PEND (the external action it waits for was not done)
eq_or_pend() { if [ "$2" = "$3" ]; then ok "$1: $2"; elif [ $MODE = dry ]; then pend "$1: '$2' (expected '$3' after the external action)"; else bad "$1: '$2' != '$3'"; fi; }
pycheck() { local out; out=$("$@" 2>&1 | mask) || true; echo "$out"
  P=$((P + $(grep -c '^  PASS' <<<"$out" || true))); N=$((N + $(grep -c '^  PEND' <<<"$out" || true)))
  local f; f=$(grep -c '^  FAIL' <<<"$out" || true); F=$((F + f)); STEP_F=$((STEP_F + f)); }
begin() { STEP=$1; STEP_F=0; P0=$P; F0=$F; N0=$N; echo; echo "================ STEP $1: $2 ($MODE, $(date -u +%FT%TZ))"; }
end_step() { local s="STEP $STEP: PASS $((P - P0)) FAIL $((F - F0)) PEND $((N - N0))"; echo "---------------- $s"; SUMMARY+=("$s"); }
# external <name> <command...>: printed in a dry run; with --execute only after all checks of the step PASS and the
# name typed back
external() { local name=$1; shift
  if [ $MODE = dry ]; then echo "  WOULD RUN [$name]: $*"; return 0; fi
  if [ "$STEP_F" != 0 ]; then echo "  REFUSED [$name]: $STEP_F check(s) of step $STEP failed"; exit 1; fi
  local ans; read -r -p "  type '$name' to run: $* > " ans </dev/tty
  [ "$ans" = "$name" ] || { echo "  NOT CONFIRMED [$name] - stop"; exit 1; }
  echo "  RUN [$name]"; "$@"; }
want() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }

# secret / identity scan of a text stream on stdin: prints "<pattern>: <count>" (counts only, never the matches)
SCAN_PATTERNS=('gmail' 'yasha1971@' 'ghp_[A-Za-z0-9]{20}' 'github_pat_' 'gh[osur]_[A-Za-z0-9]{20}' 'AKIA[0-9A-Z]{16}' 'PRIVATE KEY' 'api[_-]?key *[:=]' 'password *[:=]' 'token *= *["'"'"'][A-Za-z0-9]{16}' '/home/aeterna')
scan_stream() { local f=$TMP/scan.txt; cat > "$f"
  for p in "${SCAN_PATTERNS[@]}"; do echo "$p: $(grep -c -E -e "$p" "$f" || true)"; done
  echo "zenodo token (exact): $(grep -c -F -f <(tr -d '\n\r ' < "$TOKFILE") "$f" || true)"
  echo "gh token (exact): $(grep -c -F -f <(gh auth token 2>/dev/null | tr -d '\n\r ') "$f" || true)"; rm -f "$f"; }
# judge the scan: every count 0 except /home/aeterna (local paths in logs: NOTE, not a secret)
judge_scan() { local what=$1 line pat cnt
  while IFS= read -r line; do pat=${line%: *}; cnt=${line##*: }
    if [ "$pat" = /home/aeterna ]; then note "$what: $line (local paths, not secrets)"
    elif [ "$cnt" = 0 ]; then ok "$what: $line"; else bad "$what: $line"; fi; done; }

# ================================================================== STEP 1
step1() {
  begin 1 "aceapex: tag $SRC_TAG on $SRC_COMMIT + push of the tag"
  local A=$ACEAPEX s have
  chk "commit $SRC_COMMIT exists in $A" git -C "$A" cat-file -e "$SRC_COMMIT^{commit}"
  eq "root tree" "$(git -C "$A" rev-parse "$SRC_COMMIT^{tree}")" "$SRC_TREE"
  eq "orphan (no parents)" "$(git -C "$A" rev-list --parents -n 1 "$SRC_COMMIT")" "$SRC_COMMIT"
  eq "root entries" "$(git -C "$A" ls-tree --name-only "$SRC_COMMIT" | tr '\n' ' ')" "research "
  eq "research/ entries" "$(git -C "$A" ls-tree --name-only "$SRC_COMMIT:research" | tr '\n' ' ')" "bench2610 cleanroom refrel "
  for s in refrel cleanroom bench2610; do
    have=$(git -C "$A" rev-parse "$SRC_COMMIT:research/$s")
    eq "research/$s == ${SUB_REF[$s]}:research/$s" "$have" "$(git -C "$A" rev-parse "${SUB_REF[$s]}:research/$s")"
    eq "research/$s == stated tree" "$have" "${SUB_TREE[$s]}"
  done
  eq "files in the tree" "$(git -C "$A" ls-tree -r "$SRC_COMMIT" | wc -l)" 414
  eq "author / committer" "$(git -C "$A" log -1 --format='%an <%ae> / %cn <%ce>' "$SRC_COMMIT")" "Yakiv Shavidze <aceapex.contour@outlook.com> / Yakiv Shavidze <aceapex.contour@outlook.com>"
  git -C "$A" grep -h -I -e '' "$SRC_COMMIT" | scan_stream | judge_scan "sources tree scan"
  git -C "$A" log -1 --format='%B' "$SRC_COMMIT" | scan_stream | judge_scan "sources commit message scan"
  local loc rem
  loc=$(git -C "$A" rev-parse -q --verify "refs/tags/$SRC_TAG^{commit}" || true)
  rem=$(git -C "$A" ls-remote --tags origin "refs/tags/$SRC_TAG" "refs/tags/$SRC_TAG^{}" | awk '{print $1}' | tail -1)
  if [ -z "$loc" ]; then ok "local tag $SRC_TAG absent (to be created)"; elif [ "$loc" = "$SRC_COMMIT" ]; then ok "local tag $SRC_TAG -> $SRC_COMMIT (already created)"; else bad "local tag $SRC_TAG -> $loc"; fi
  if [ -z "$rem" ]; then ok "remote tag $SRC_TAG absent (to be pushed)"; elif [ "$rem" = "$SRC_COMMIT" ]; then ok "remote tag $SRC_TAG -> $SRC_COMMIT (already pushed)"; else bad "remote tag $SRC_TAG -> $rem"; fi
  # links of README / EVIDENCE / DATASET against the commit itself (the tag does not exist before this step)
  if python3 "$PANVRAM/scripts/check_source_links.py" "$A" "$SRC_COMMIT"; then ok "check_source_links against $SRC_COMMIT: 0 broken"; else bad "check_source_links against $SRC_COMMIT"; fi
  [ -n "$loc" ] || external tag-sources git -C "$A" tag "$SRC_TAG" "$SRC_COMMIT"
  [ "$rem" = "$SRC_COMMIT" ] || external push-sources-tag git -C "$A" push origin "refs/tags/$SRC_TAG"
  echo "  after:"
  rem=$(git -C "$A" ls-remote --tags origin "refs/tags/$SRC_TAG" "refs/tags/$SRC_TAG^{}" | awk '{print $1}' | tail -1)
  eq_or_pend "git ls-remote --tags origin $SRC_TAG" "$rem" "$SRC_COMMIT"
  if git -C "$A" rev-parse -q --verify "refs/tags/$SRC_TAG" >/dev/null; then
    if python3 "$PANVRAM/scripts/check_source_links.py" "$A" "$SRC_TAG"; then ok "check_source_links against $SRC_TAG: 0 broken"; else bad "check_source_links against $SRC_TAG"; fi
  else eq_or_pend "check_source_links against $SRC_TAG" "tag absent" "0 broken"; fi
  end_step
}

# ================================================================== STEP 2
# release_edit <tree>: DOI patch + release date in CHANGELOG.md and CITATION.cff
release_edit() { local T=$1 patch=$2
  if (cd "$T" && python3 zenodo/doi_patch.py "$DOI" "$patch"); then ok "doi_patch.py -> $patch ($(wc -l < "$patch") lines, sha256 $(sha256sum "$patch" | cut -c1-16)...)"; else bad "doi_patch.py"; return; fi
  if git -C "$T" apply --check "$patch" && git -C "$T" apply "$patch"; then ok "git apply"; else bad "git apply"; return; fi
  if (cd "$T" && python3 - "$DATE" <<'PY'
import sys
from pathlib import Path
d = sys.argv[1]
for name, old, new in (("CHANGELOG.md", "## 1.0.0 - (date set at release)\n", f"## 1.0.0 - {d}\n"),
                       ("CITATION.cff", 'date-released: "2026-10-08"\n', f'date-released: "{d}"\n')):
    p = Path(name); s = p.read_text()
    if s.count(old) != 1: sys.exit(f"{name}: expected line not found once")
    p.write_text(s.replace(old, new))
PY
  ); then ok "CHANGELOG.md '## 1.0.0 - $DATE', CITATION.cff date-released $DATE"; else bad "release date edit"; fi
}
release_content_checks() { local T=$1 n
  n=$(git -C "$T" grep -n -I -E 'XXXXXXX|REPLACE-RECORD-ID|date set at release' -- . ':!zenodo/doi_patch.py' ':!zenodo/release_phase2.sh' ':!zenodo/RELEASE_PHASE2.md' | tee "$TMP/placeholders.txt" | wc -l)
  if [ "$n" = 0 ]; then ok "no placeholders left (XXXXXXX, REPLACE-RECORD-ID, date set at release)"; else bad "$n placeholder line(s) left:"; sed 's/^/          /' "$TMP/placeholders.txt"; fi
  chk "README.md names the DOI and the repository" grep -q "Repository: $GH_URL. Dataset .*DOI $DOI" "$T/README.md"
  chk "CITATION.cff: DOI identifier, repository-code, date-released $DATE" bash -c "grep -q 'value: \"$DOI\"' '$T/CITATION.cff' && grep -q 'repository-code: \"$GH_URL\"' '$T/CITATION.cff' && grep -q 'date-released: \"$DATE\"' '$T/CITATION.cff'"
  chk ".zenodo.json: valid, related DOI $DOI, repository URL" python3 -c "import json,sys; z=json.load(open('$T/.zenodo.json')); ids=[r['identifier'] for r in z['related_identifiers']]; sys.exit(not ('$DOI' in ids and '$GH_URL' in ids))"
  chk "quickstart_558.ipynb: valid JSON, DATASET_BASE_URL = records/$RECID, PANVRAM_SOURCE = $ARCHIVE_URL" python3 -c "import json,sys; s=''.join(''.join(c['source']) for c in json.load(open('$T/notebooks/quickstart_558.ipynb'))['cells']); sys.exit(not ('DATASET_BASE_URL = \"https://zenodo.org/records/$RECID/files\"' in s and 'PANVRAM_SOURCE = \"$ARCHIVE_URL\"' in s and '$MANIFEST_SHA' in s))"
  chk "CHANGELOG.md: '## 1.0.0 - $DATE'" grep -qx "## 1.0.0 - $DATE" "$T/CHANGELOG.md"
}
# run_tests <tree>: the CI test set (CPU path) against the extension built in place in <tree>
run_tests() { local T=$1
  if [ ! -x "$PY" ]; then bad "test interpreter $PY missing"; return; fi
  if (cd "$T" && PANVRAM_CUDA=0 "$PY" setup.py -q build_ext --inplace > "$TMP/build.txt" 2>&1); then ok "build_ext --inplace (CPU path) in $T"; else bad "build (tail below)"; tail -20 "$TMP/build.txt"; return; fi
  eq "panvram imported from" "$(cd "$T/tests" && PYTHONPATH=$T "$PY" -c 'import panvram, os; print(os.path.dirname(panvram.__file__))')" "$T/panvram"
  local t r
  for t in test_synth.py test_gate.py; do
    if r=$(cd "$T/tests" && ulimit -c 0 && PYTHONPATH=$T "$PY" -m pytest -q -rs -p no:cacheprovider "$t" 2>&1); then ok "pytest $t: $(tail -1 <<<"$r")"; else bad "pytest $t: $(tail -1 <<<"$r")"; tail -30 <<<"$r"; fi
  done
}
release_commit() { git -C "$PANVRAM" add -- "${RELEASE_PATHS[@]}" && git -C "$PANVRAM" "${GIT_ID[@]}" commit -q -m "$RELEASE_MSG" -- "${RELEASE_PATHS[@]}"; }
step2() {
  begin 2 "panvram: DOI patch, release date $DATE, tests, commit, tag $TAG, tarball, push"
  local T C rem_main
  chk "panvram on branch master" test "$(git -C "$PANVRAM" symbolic-ref --short HEAD)" = master
  eq "clean tree (git status --porcelain)" "$(git -C "$PANVRAM" status --porcelain | wc -l)" 0
  rem_main=$(git -C "$PANVRAM" ls-remote origin refs/heads/main | awk '{print $1}')
  if git -C "$PANVRAM" rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then
    note "local tag $TAG exists -> $(git -C "$PANVRAM" rev-parse "$TAG^{commit}"): step 2 was done; only the checks after run"
    C=$(git -C "$PANVRAM" rev-parse "$TAG^{commit}"); T=$PANVRAM
  else
    chk "origin/main $rem_main is an ancestor of HEAD (the push of step 2 is a fast-forward)" git -C "$PANVRAM" merge-base --is-ancestor "$rem_main" HEAD
    note "commits on master not on origin/main yet (pushed with the release commit): $(git -C "$PANVRAM" rev-list --count "$rem_main..HEAD")"
    eq "remote tag $TAG absent" "$(git -C "$PANVRAM" ls-remote --tags origin "refs/tags/$TAG" | wc -l)" 0
    chk "CHANGELOG.md has '## 1.0.0 - (date set at release)'" grep -qx '## 1.0.0 - (date set at release)' "$PANVRAM/CHANGELOG.md"
    if [ $MODE = dry ]; then
      [ ! -e "$WT" ] || { bad "$WT exists (not created by this run) - not touching it"; end_step; return; }
      git -C "$PANVRAM" worktree add --detach "$WT" HEAD >/dev/null 2>&1 && WT_MADE=1 && ok "throwaway worktree $WT at $(git -C "$WT" rev-parse --short HEAD)" || { bad "worktree add"; end_step; return; }
      T=$WT; release_edit "$T" "$OUT/dryrun-doi_23232317.patch"
    else
      T=$PANVRAM; release_edit "$T" "$OUT/doi_23232317.patch"
    fi
    release_content_checks "$T"
    echo "  diff (git diff --stat; full diff below):"; git -C "$T" diff --stat | sed 's/^/    /'
    git -C "$T" diff | sed 's/^/    | /'
    [ $MODE = dry ] && cp "$T/notebooks/quickstart_558.ipynb" "$OUT/dryrun-quickstart_558.patched.ipynb"
    run_tests "$T"
    local msg="panvram 1.0.0: dataset DOI $DOI and repository URL (CITATION.cff, .zenodo.json, README, quickstart notebook), release date $DATE

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
    if [ $MODE = dry ]; then
      # commit on the detached HEAD of the throwaway worktree: a preview object only, gone with the worktree
      git -C "$T" add -- "${RELEASE_PATHS[@]}" && git -C "$T" "${GIT_ID[@]}" commit -q -m "$msg" -- "${RELEASE_PATHS[@]}" && ok "preview commit in the worktree: $(git -C "$T" rev-parse --short HEAD) ($(git -C "$T" show --stat --format= HEAD | tail -1))" || bad "preview commit"
      eq "tree clean after the commit (only the 5 release paths changed)" "$(git -C "$T" status --porcelain --untracked-files=no | wc -l)" 0
      C=$(git -C "$T" rev-parse HEAD)
      echo "  WOULD RUN [release-commit]: git -C $PANVRAM add ${RELEASE_PATHS[*]} && git ${GIT_ID[*]} commit -m '<message above>' -- ${RELEASE_PATHS[*]}"
      external tag-v1.0.0 git -C "$PANVRAM" "${GIT_ID[@]}" tag -a "$TAG" -m "panvram 1.0.0" HEAD
      local s; s=$(git -C "$T" archive --format=tar --prefix="$PREFIX/" "$C" | gzip -n | sha256sum | cut -d' ' -f1)
      note "PREVIEW tarball sha256 (preview commit; the real one differs - other commit time and id): $s"
      eq "PREVIEW tarball: MANIFEST.tsv sha256" "$(git -C "$T" archive --format=tar --prefix="$PREFIX/" "$C" | tar -xOf - "$PREFIX/MANIFEST.tsv" | sha256sum | cut -d' ' -f1)" "$MANIFEST_SHA"
      if python3 "$T/scripts/check_source_links.py" "$ACEAPEX" "$SRC_COMMIT"; then ok "check_source_links (patched README/EVIDENCE/DATASET) against $SRC_COMMIT: 0 broken"; else bad "check_source_links (patched)"; fi
      external tarball bash -c "git -C $PANVRAM archive --format=tar --prefix=$PREFIX/ $TAG^{commit} | gzip -n > $TARBALL && (cd $OUT && sha256sum $PREFIX.tar.gz > $PREFIX.tar.gz.sha256)"
      external push-v1.0.0 git -C "$PANVRAM" push origin master:main "refs/tags/$TAG"
      git -C "$PANVRAM" worktree remove --force "$WT" && git -C "$PANVRAM" worktree prune && WT_MADE=0 && ok "throwaway worktree removed"
      eq "panvram tree unchanged by the dry run" "$(git -C "$PANVRAM" status --porcelain | wc -l)" 0
      end_step; return
    fi
    RELEASE_MSG=$msg
    external release-commit release_commit
    external tag-v1.0.0 git -C "$PANVRAM" "${GIT_ID[@]}" tag -a "$TAG" -m "panvram 1.0.0" HEAD
    C=$(git -C "$PANVRAM" rev-parse "$TAG^{commit}")
  fi
  # tarball from the tag commit (local)
  if [ ! -s "$TARBALL" ]; then
    git -C "$PANVRAM" archive --format=tar --prefix="$PREFIX/" "$C" | gzip -n > "$TARBALL"
    (cd "$OUT" && sha256sum "$PREFIX.tar.gz" > "$PREFIX.tar.gz.sha256")
  fi
  eq "tarball commit id (git get-tar-commit-id)" "$(gzip -dc "$TARBALL" | git get-tar-commit-id)" "$C"
  chk "tarball .sha256" bash -c "cd '$OUT' && sha256sum -c --quiet '$PREFIX.tar.gz.sha256'"
  eq "tarball rebuilt from $C == $TARBALL" "$(git -C "$PANVRAM" archive --format=tar --prefix="$PREFIX/" "$C" | gzip -n | sha256sum | cut -d' ' -f1)" "$(sha256sum "$TARBALL" | cut -d' ' -f1)"
  eq "tarball MANIFEST.tsv sha256" "$(tar -xzOf "$TARBALL" "$PREFIX/MANIFEST.tsv" | sha256sum | cut -d' ' -f1)" "$MANIFEST_SHA"
  note "tarball $TARBALL sha256 $(sha256sum "$TARBALL" | cut -d' ' -f1), $(stat -c %s "$TARBALL") B"
  release_content_checks "$PANVRAM"
  local ref=$SRC_TAG; git -C "$ACEAPEX" rev-parse -q --verify "refs/tags/$SRC_TAG" >/dev/null || ref=$SRC_COMMIT
  if python3 "$PANVRAM/scripts/check_source_links.py" "$ACEAPEX" "$ref"; then ok "check_source_links against $ref: 0 broken"; else bad "check_source_links against $ref"; fi
  [ "$(git -C "$PANVRAM" rev-parse "$TAG^{commit}")" = "$(git -C "$PANVRAM" rev-parse HEAD)" ] || bad "HEAD is not the tag commit"
  [ "$(git -C "$PANVRAM" ls-remote origin refs/heads/main | awk '{print $1}')" = "$C" ] || external push-v1.0.0 git -C "$PANVRAM" push origin master:main "refs/tags/$TAG"
  echo "  after:"
  eq "origin main == tag commit" "$(git -C "$PANVRAM" ls-remote origin refs/heads/main | awk '{print $1}')" "$C"
  eq "origin $TAG^{} == tag commit" "$(git -C "$PANVRAM" ls-remote --tags origin "refs/tags/$TAG^{}" | awk '{print $1}')" "$C"
  end_step
}

# ================================================================== STEP 3
step3() {
  begin 3 "GitHub: $GH public, release $TAG"
  local C priv notes runs
  chk "gh authenticated" gh auth status
  priv=$(gh api "repos/$GH" --jq .private 2>/dev/null || echo "?")
  if [ "$priv" = true ]; then ok "repository private now (to be made public)"; elif [ "$priv" = false ]; then note "repository already public"; else bad "GET repos/$GH: $priv"; fi
  C=$(git -C "$PANVRAM" rev-parse -q --verify "refs/tags/$TAG^{commit}" || true)
  if [ -n "$C" ]; then
    eq "origin main == $TAG commit" "$(git -C "$PANVRAM" ls-remote origin refs/heads/main | awk '{print $1}')" "$C"
    eq "origin $TAG^{} == $TAG commit" "$(git -C "$PANVRAM" ls-remote --tags origin "refs/tags/$TAG^{}" | awk '{print $1}')" "$C"
  else
    eq_or_pend "local tag $TAG" "absent" "present"; C=$(git -C "$PANVRAM" rev-parse HEAD)
    note "CI below is checked on HEAD $C (the tag commit does not exist yet)"
  fi
  runs=$(gh api "repos/$GH/commits/$C/check-runs" --jq '.check_runs[] | "\(.name)=\(.status)/\(.conclusion)"' 2>/dev/null | sort | tr '\n' ' ' || true)
  if [ -n "$runs" ] && ! grep -qv 'completed/success' <<<"$(tr ' ' '\n' <<<"$runs" | sed '/^$/d')"; then ok "CI on $C: $runs"; else bad "CI on $C: '${runs:-no check runs}'"; fi
  echo "  history scan of panvram (git log --all -p: every commit, diff and message of every local and remote-tracking ref):"
  echo "    refs scanned: $(git -C "$PANVRAM" for-each-ref --format='%(refname:short)' | tr '\n' ' ')"
  echo "    refs on GitHub (become public): $(git -C "$PANVRAM" ls-remote origin | awk '{print $2}' | tr '\n' ' ')"
  echo "    commits: $(git -C "$PANVRAM" rev-list --all | wc -l); author/committer e-mails: $(git -C "$PANVRAM" log --all --format='%ae%n%ce' | sort | uniq -c | tr -s ' ' | tr '\n' ';')"
  git -C "$PANVRAM" log --all -p --format='%H%n%an <%ae>%n%cn <%ce>%n%B' --no-ext-diff | scan_stream | judge_scan "panvram history"
  git -C "$PANVRAM" ls-files -z | (cd "$PANVRAM" && xargs -0 cat 2>/dev/null) | scan_stream | judge_scan "panvram tracked files (HEAD tree)"
  # release notes = CHANGELOG section 1.0.0
  notes=$OUT/release-notes-$TAG.md; local cl=$PANVRAM/CHANGELOG.md
  [ $MODE = dry ] && [ -z "$(git -C "$PANVRAM" rev-parse -q --verify "refs/tags/$TAG" || true)" ] && notes=$OUT/dryrun-release-notes-$TAG.md
  if [ "$notes" = "$OUT/release-notes-$TAG.md" ]; then git -C "$PANVRAM" show "$TAG:CHANGELOG.md" > "$TMP/cl.md"; else sed "s/^## 1.0.0 - (date set at release)\$/## 1.0.0 - $DATE/" "$cl" > "$TMP/cl.md"; fi
  awk '/^## 1\.0\.0 - /{on=1} on&&/^## /&&!/^## 1\.0\.0 - /{exit} on' "$TMP/cl.md" | sed '1s/^## /# panvram /' > "$notes"
  if [ -s "$notes" ] && head -1 "$notes" | grep -q "^# panvram 1.0.0 - $DATE"; then ok "release notes $notes ($(wc -l < "$notes") lines, from CHANGELOG 1.0.0)"; else bad "release notes $notes"; fi
  if [ -s "$TARBALL" ] && [ -s "$TARBALL.sha256" ] && (cd "$OUT" && sha256sum -c --quiet "$PREFIX.tar.gz.sha256") && [ "$(gzip -dc "$TARBALL" | git get-tar-commit-id)" = "$C" ]; then ok "release assets $TARBALL(.sha256) of $C"
  else eq_or_pend "release assets $TARBALL(.sha256) of the tag commit" "absent or other commit" "present"; fi
  # GitHub archive vs git archive: same files, different bytes (GitHub's gzip / tar header); checked on $C
  local ga="$TMP/gh.tar.gz"
  if gh api "repos/$GH/tarball/$C" > "$ga" 2>/dev/null; then
    mkdir -p "$TMP/a" "$TMP/b"; tar -xzf "$ga" -C "$TMP/a" --strip-components=1; git -C "$PANVRAM" archive --format=tar "$C" | tar -xf - -C "$TMP/b"
    if diff -r "$TMP/a" "$TMP/b" >/dev/null; then ok "GitHub tarball of $C has the same files as git archive ($(find "$TMP/b" -type f | wc -l) files); bytes: GitHub $(sha256sum "$ga" | cut -c1-16)... vs gzip -n $(git -C "$PANVRAM" archive --format=tar --prefix="$PREFIX/" "$C" | gzip -n | sha256sum | cut -c1-16)... (differ by design)"
    else bad "GitHub tarball content != git archive of $C"; fi
    rm -rf "$TMP/a" "$TMP/b"
  else bad "GET repos/$GH/tarball/$C"; fi
  note "the notebook downloads PANVRAM_SOURCE = $ARCHIVE_URL (GitHub archive, not the release asset); it checks MANIFEST.tsv SHA-256 after extraction, not the tarball bytes"
  note "decision: GitHub-Zenodo integration (software DOI from .zenodo.json) is NOT enabled by this script; enable it in Zenodo -> GitHub before 'gh-release' if a software record is wanted"
  external make-public gh api -X PATCH "repos/$GH" -F private=false
  external gh-release gh release create "$TAG" "$TARBALL" "$TARBALL.sha256" --verify-tag --title "panvram 1.0.0" --notes-file "$notes"
  echo "  after:"
  eq_or_pend "gh api repos/$GH --jq .private" "$(gh api "repos/$GH" --jq .private 2>/dev/null || echo '?')" false
  eq_or_pend "anonymous GET api.github.com/repos/$GH" "$(curl -s -o /dev/null -w '%{http_code}' "https://api.github.com/repos/$GH")" 200
  local code
  code=$(curl -sL -o "$TMP/asset.tar.gz" -w '%{http_code}' "$ASSET_URL" || true)
  if [ "$code" = 200 ] && [ -s "$TARBALL" ]; then eq "release asset sha256 (anonymous download) == local" "$(sha256sum < "$TMP/asset.tar.gz" | cut -d' ' -f1)" "$(sha256sum < "$TARBALL" | cut -d' ' -f1)"
  else eq_or_pend "release asset $ASSET_URL" "HTTP $code" "HTTP 200"; fi
  code=$(curl -sL -o "$TMP/archive.tar.gz" -w '%{http_code}' "$ARCHIVE_URL" || true)
  if [ "$code" = 200 ]; then
    local s; s=$(sha256sum < "$TMP/archive.tar.gz" | cut -d' ' -f1); echo "$s  $TAG.tar.gz ($ARCHIVE_URL, $(date -u +%FT%TZ))" > "$OUT/github-archive-$TAG.tar.gz.sha256"
    ok "archive URL sha256 recorded: $s -> $OUT/github-archive-$TAG.tar.gz.sha256"
    mkdir -p "$TMP/a" "$TMP/b"; tar -xzf "$TMP/archive.tar.gz" -C "$TMP/a" --strip-components=1; tar -xzf "$TARBALL" -C "$TMP/b" --strip-components=1
    if diff -r "$TMP/a" "$TMP/b" >/dev/null; then ok "archive URL files == release asset files"; else bad "archive URL files != release asset files"; fi
    eq "archive URL MANIFEST.tsv sha256" "$(sha256sum < "$TMP/a/MANIFEST.tsv" | cut -d' ' -f1)" "$MANIFEST_SHA"
    rm -rf "$TMP/a" "$TMP/b"
  else eq_or_pend "archive URL $ARCHIVE_URL" "HTTP $code" "HTTP 200"; fi
  end_step
}

# ================================================================== STEP 4
zenodo_py() { python3 - "$@" <<'PY'
import hashlib, html, json, os, re, sys, urllib.request, urllib.error
from pathlib import Path
cmd, recid, doi, files_md, meta_json, mode = sys.argv[1:7]
TOK = Path("~/.zenodo_token").expanduser().read_text().strip()
mask = lambda s: str(s).replace(TOK, "<token>") if TOK else str(s)
DEP = f"https://zenodo.org/api/deposit/depositions/{recid}"
def call(method, url, auth=True):
    r = urllib.request.Request(url, method=method)
    if auth: r.add_header("Authorization", "Bearer " + TOK)
    try:
        with urllib.request.urlopen(r, timeout=180) as f: return f.status, f.read()
    except urllib.error.HTTPError as e: return e.code, e.read()
    except Exception as e: return -1, mask(f"{type(e).__name__}: {e}").encode()
def res(ok, msg, pend=False):
    print(("  PASS  " if ok else ("  PEND  " if pend and mode == "dry" else "  FAIL  ")) + mask(msg))
if cmd == "before":
    st, body = call("GET", DEP)
    if st != 200: res(False, f"GET deposition {recid}: HTTP {st} {mask(body[:200])}"); sys.exit()
    d = json.loads(body); m = d["metadata"]
    res(d["state"] == "unsubmitted" and d["submitted"] is False, f"state {d['state']}, submitted {d['submitted']}")
    res(m.get("prereserve_doi", {}).get("doi") == doi, f"reserved DOI {m.get('prereserve_doi', {}).get('doi')}")
    server = {(f["filename"], int(f["filesize"]), f["checksum"]) for f in d["files"]}
    table = set()
    for line in Path(files_md).read_text().splitlines():
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if len(c) == 5 and c[0].isdigit(): table.add((c[1], int(c[2].replace(" ", "")), c[3]))
    res(len(d["files"]) == 48 and len(server) == 48, f"files on the draft: {len(d['files'])} (48 expected)")
    res(server == table, f"file list (name, size, md5) == FILES.md table ({len(table)} rows; only on server {len(server - table)}, only in table {len(table - server)})")
    print(f"  NOTE  total {sum(s for _, s, _ in server)} B")
    want = json.loads(Path(meta_json).read_text())
    for k in ("title", "version", "license", "upload_type", "keywords", "notes"):
        res(m.get(k) == want.get(k), f"metadata {k} == dataset.zenodo.json")
    norm = lambda s: re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", s or ""))).strip()
    res(norm(m.get("description")) == norm(want.get("description")), "metadata description == dataset.zenodo.json (HTML tags and whitespace normalised)")
    cr = lambda L: [(c.get("name"), c.get("orcid")) for c in L]
    res(cr(m.get("creators", [])) == cr(want["creators"]), f"metadata creators == dataset.zenodo.json {cr(want['creators'])}")
    ri = lambda L: sorted((r["identifier"], r["relation"], r.get("scheme"), r.get("resource_type")) for r in L)
    res(ri(m.get("related_identifiers", [])) == ri(want["related_identifiers"]), f"metadata related_identifiers == dataset.zenodo.json ({len(want['related_identifiers'])})")
    res(m.get("access_right") == "open", f"access_right {m.get('access_right')}")
    print(f"  NOTE  publication_date on the draft: {m.get('publication_date')} (kept as is by publish; edit before publishing if the release day differs)")
    gh = [r["identifier"] for r in want["related_identifiers"] if r["identifier"].startswith("https://github.com/")]
    for u in gh:
        st, _ = call("GET", u, auth=False)
        res(st == 200, f"anonymous GET {u}: HTTP {st} (public after step 3)", pend=True)
elif cmd == "publish":
    st, body = call("POST", DEP + "/actions/publish")
    print(f"  publish: HTTP {st}")
    if st not in (200, 202): print("  " + mask(body[:500].decode(errors="replace"))); sys.exit(1)
    d = json.loads(body); print(f"  state {d.get('state')}, submitted {d.get('submitted')}, doi {d.get('doi')}, record {d.get('record_id')}")
elif cmd == "after":
    st, body = call("GET", f"https://zenodo.org/api/records/{recid}", auth=False)
    pub = st == 200 and json.loads(body).get("status", "published") in ("published",) and json.loads(body).get("access", {}).get("record", "public") == "public"
    res(pub, f"anonymous GET api/records/{recid}: HTTP {st}" + (f", files {len(json.loads(body).get('files', []))}" if st == 200 else ""), pend=True)
    st, _ = call("GET", f"https://zenodo.org/records/{recid}", auth=False)
    res(st == 200, f"anonymous GET https://zenodo.org/records/{recid}: HTTP {st}", pend=True)
    r = urllib.request.Request(f"https://doi.org/{doi}", method="HEAD")
    class NoRedir(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, *a, **k): return None
    try: st, loc = urllib.request.build_opener(NoRedir).open(r, timeout=60).status, ""
    except urllib.error.HTTPError as e: st, loc = e.code, e.headers.get("Location", "")
    except Exception as e: st, loc = -1, str(e)
    res(st in (301, 302, 303, 307, 308) and "zenodo" in loc, f"https://doi.org/{doi}: HTTP {st} -> {loc or '-'}", pend=True)
    st, body = call("GET", f"https://zenodo.org/records/{recid}/files/MANIFEST.tsv?download=1", auth=False)
    have = hashlib.sha256(body).hexdigest() if st == 200 else f"HTTP {st}"
    res(have == "3b6941efeccb00a7b092573b0482281711f249e62465f771325704d9334d31d0", f"anonymous MANIFEST.tsv sha256: {have}", pend=True)
PY
}
step4() {
  begin 4 "Zenodo: publish draft $RECID (DOI $DOI) - the LAST external action"
  chk "token file $TOKFILE present (mode $(stat -c %a "$TOKFILE" 2>/dev/null))" test -s "$TOKFILE"
  pycheck zenodo_py before "$RECID" "$DOI" "$PANVRAM/zenodo/FILES.md" "$PANVRAM/zenodo/dataset.zenodo.json" "$MODE"
  if [ $MODE = exec ]; then
    eq "aceapex $SRC_TAG on origin" "$(git -C "$ACEAPEX" ls-remote --tags origin "refs/tags/$SRC_TAG" "refs/tags/$SRC_TAG^{}" | awk '{print $1}' | tail -1)" "$SRC_COMMIT"
    eq "GitHub release $TAG exists" "$(gh release view "$TAG" -R "$GH" --json tagName --jq .tagName 2>/dev/null || echo absent)" "$TAG"
  fi
  external zenodo-publish zenodo_py publish "$RECID" "$DOI" - - "$MODE"
  echo "  after:"
  pycheck zenodo_py after "$RECID" "$DOI" - - "$MODE"
  end_step
}

# ================================================================== STEP 5
# Colab final external check, generated from the DOI-patched public notebook (tag v1.0.0 if it exists, else HEAD +
# the DOI patch of step 2): cells 2-4 verbatim through a step runner; the gate with forms 1,1; log on Drive.
step5() {
  begin 5 "Colab: final external check from the public URLs (cell text + notebook)"
  local NB=$TMP/quickstart_558.ipynb src
  if git -C "$PANVRAM" rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then git -C "$PANVRAM" show "$TAG:notebooks/quickstart_558.ipynb" > "$NB"; src="$TAG"
  else
    local patch=$OUT/doi_23232317.patch; [ -s "$patch" ] || patch=$OUT/dryrun-doi_23232317.patch
    if [ ! -s "$patch" ]; then bad "no DOI patch ($OUT/[dryrun-]doi_23232317.patch): run step 2 first"; end_step; return; fi
    mkdir -p "$TMP/nb/notebooks"; git -C "$PANVRAM" show HEAD:notebooks/quickstart_558.ipynb > "$TMP/nb/notebooks/quickstart_558.ipynb"
    if (cd "$TMP/nb" && git apply --include=notebooks/quickstart_558.ipynb "$patch"); then cp "$TMP/nb/notebooks/quickstart_558.ipynb" "$NB"; src="HEAD $(git -C "$PANVRAM" rev-parse --short HEAD) + $(basename "$patch")"; ok "patched notebook from $src"
    else bad "applying the notebook hunk of $patch"; end_step; return; fi
  fi
  mkdir -p "$OUT"
  if python3 "$PANVRAM/zenodo/colab_public_check.py" "$NB" "$src" "$OUT/colab_final_external_check.txt" "$OUT/quickstart_public_check.ipynb"; then
    ok "generated $OUT/colab_final_external_check.txt ($(sha256sum < "$OUT/colab_final_external_check.txt" | cut -c1-16)...) and $OUT/quickstart_public_check.ipynb ($(sha256sum < "$OUT/quickstart_public_check.ipynb" | cut -c1-16)...)"
  else bad "colab_public_check.py"; fi
  chk "cell text compiles (python)" python3 -c "import ast,sys; ast.parse(open('$OUT/colab_final_external_check.txt').read())"
  chk "cell text: no token, no getpass, no deposit API" bash -c "! grep -qE 'getpass|api/deposit|zenodo_token|Bearer' '$OUT/colab_final_external_check.txt'"
  chk "cell text ends with runtime.unassign()" bash -c "tail -1 '$OUT/colab_final_external_check.txt' | grep -q 'runtime.unassign()'"
  chk "notebook: valid JSON, last cell ends with runtime.unassign()" python3 -c "import json; c=json.load(open('$OUT/quickstart_public_check.ipynb'))['cells']; assert ''.join(c[-1]['source']).rstrip().splitlines()[-1].strip().endswith('runtime.unassign()')"
  note "run only after step 4: Colab GPU runtime, paste the cell (or upload the notebook), Run all; log ACEAPEX-OVH/quickstart_public_<UTC>.txt"
  end_step
}

echo "release_phase2.sh mode=$MODE step=${ONLY:-all} date=$DATE  panvram $(git -C "$PANVRAM" rev-parse --short HEAD)  aceapex $(git -C "$ACEAPEX" rev-parse --short HEAD)  $(date -u +%FT%TZ)"
want 1 && step1
want 2 && step2
want 3 && step3
want 4 && step4
want 5 && step5
echo; echo "================ SUMMARY ($MODE)"; printf '%s\n' "${SUMMARY[@]}"
echo "TOTAL: PASS $P FAIL $F PEND $N"
[ "$F" = 0 ]
