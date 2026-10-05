#!/usr/bin/env bash
# run_colab_cohort558.sh - the whole HPRC cohort (558 assemblies, refrel3 v1 q4k) resident on one Colab GPU, measured and
# checked without source files. From Drive (uploaded from ace-core, nothing fetched but T2T if absent):
#   MyDrive/panvram/panvram.tar.gz          the package (git archive, the commit inside)
#   MyDrive/panvram/cohort558/              558 x <name>.q4k.rr3 + SHA256SUMS, or instead
#   MyDrive/panvram/cohort558_q4k.tar       the same in one uncompressed tar + cohort558_q4k.tar.sha256 next to it
#                                           (the tar's sha256 checked first, then unpacked to the VM disk)
#   MyDrive/panvram/manifest_v1.tsv         names, source hashes, fasta_xxh3, archive sha256 (aceapex research/refrel)
#   MyDrive/aceapex_corpus/t2t.fa[.gz]      T2T-CHM13v2.0 (md5 checked; else NCBI)
# Stages (fail-fast, the failing stage in ..._FAILED.txt):
#   S0 build: pip install -e . (CUDA, the GPU's arch)
#   S1 data: archives copied to the VM disk, sha256sum -c, SHA256SUMS == manifest's q4k_sha256 for all 558
#   S2 scripts/cohort558.py: open time, resident bytes, torch allocated, nvidia-smi before / after; sample 1024 x 8192
#      (== CPU decoder); fetch at 1000 coordinates (== CPU decoder); every assembly decoded in full on the GPU, FASTA
#      rebuilt from the contig table, XXH3 == fasta_xxh3 of the manifest (all 558, time)
# Result: MyDrive/aceapex_logs/panvram_cohort558_<date>.txt; failure: panvram_cohort558_<date>_FAILED.txt.
# Env: DRIVE, W, PV_TAR, COHORT_SRC (folder), COHORT_TAR, MANIFEST
set -Eeuo pipefail
ulimit -c 0   # no core dumps from test runs
shopt -s nullglob
DRIVE=${DRIVE:-/content/drive/MyDrive}
STORE=$DRIVE/aceapex_corpus; LOGS=$DRIVE/aceapex_logs
PV_TAR=${PV_TAR:-$DRIVE/panvram/panvram.tar.gz}
COHORT_SRC=${COHORT_SRC:-$DRIVE/panvram/cohort558}
COHORT_TAR=${COHORT_TAR:-$DRIVE/panvram/cohort558_q4k.tar}
MANIFEST=${MANIFEST:-$DRIVE/panvram/manifest_v1.tsv}
W=${W:-/content/cohort558_run}; PV=$W/panvram; C=$W/cohort
[ -d $DRIVE ] || { echo "Drive not mounted: from google.colab import drive; drive.mount('/content/drive')"; echo "DONE — выключи runtime"; exit 1; }
for f in $PV_TAR $MANIFEST; do [ -s $f ] || { echo "missing $f (upload from ace-core)"; echo "DONE — выключи runtime"; exit 1; }; done
if [ -s $COHORT_SRC/SHA256SUMS ]; then SRC=dir
elif [ -s $COHORT_TAR ] && [ -s $COHORT_TAR.sha256 ]; then SRC=tar
else echo "missing the cohort: $COHORT_SRC/ (with SHA256SUMS) or $COHORT_TAR + $COHORT_TAR.sha256"; echo "DONE — выключи runtime"; exit 1; fi
mkdir -p $W $C $LOGS
DAY=$(date -u +%Y-%m-%d); RUNLOG=$W/run-$DAY.log; OUT=$LOGS/panvram_cohort558_$DAY.txt; STAGE=setup
exec > >(tee -a $RUNLOG) 2>&1
on_err() { local rc=$? ln=$1 cmd=$2
  echo "FAILED in $STAGE: exit $rc at line $ln: $cmd"
  { echo "# panvram cohort558 FAILED $(date -u +%FT%TZ) in stage $STAGE: exit $rc at line $ln: $cmd"; tail -n 150 $RUNLOG; } > $LOGS/panvram_cohort558_${DAY}_FAILED.txt 2>/dev/null || true
  echo "details: $LOGS/panvram_cohort558_${DAY}_FAILED.txt"; echo "DONE — выключи runtime"; }
trap 'on_err $LINENO "$BASH_COMMAND"' ERR
T=$(nproc)
nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv,noheader | tee $W/gpu.txt
command -v nvcc >/dev/null || export PATH=/usr/local/cuda/bin:$PATH
nvcc --version | tail -n 2 | tee -a $W/gpu.txt
python3 -c "import torch; print('torch', torch.__version__, 'cuda', torch.version.cuda)" | tee -a $W/gpu.txt
echo "host: $T threads, $(free -g | awk '/Mem:/{print $2}') GB RAM, $(df -BG --output=avail /content | tail -1 | tr -d ' ') free on /content" | tee -a $W/gpu.txt
dpkg -s libzstd-dev >/dev/null 2>&1 || apt-get install -y -qq libzstd-dev >/dev/null

# ---------------------------------------------------------------- S0 build
STAGE=S0
rm -rf $PV; tar -xzf $PV_TAR -C $W
PV_COMMIT=$( (gzip -dc $PV_TAR 2>/dev/null || true) | git get-tar-commit-id 2>/dev/null || true); PV_COMMIT=${PV_COMMIT:-unknown}; echo "panvram $PV_COMMIT"
cd $PV; pip uninstall -y -q panvram 2>/dev/null || true
t0=$(date +%s); MAX_JOBS=$T pip install -e . > $W/build.txt 2>&1 || { tail -40 $W/build.txt; false; }
echo "S0 pip install -e . : $(( $(date +%s) - t0 )) s"
python3 -c "import panvram, torch; assert panvram.with_cuda, 'built without CUDA'; print('S0 panvram', panvram.__version__, 'device', torch.cuda.get_device_name(0), 'sm_%d%d' % torch.cuda.get_device_capability(0))"

# ---------------------------------------------------------------- S1 data
STAGE=S1
FA=$W/t2t.fa
if [ ! -s $FA ]; then
  if [ -s $STORE/t2t.fa ]; then cp $STORE/t2t.fa $FA
  elif [ -s $STORE/t2t.fa.gz ]; then gunzip -c $STORE/t2t.fa.gz > $FA
  else curl -fsSL --retry 3 -o $W/t2t.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
    gunzip -c $W/t2t.fa.gz > $FA; rm -f $W/t2t.fa.gz; fi
fi
echo "cd1e52ce400c027ed0b7ab4b9d613f5a  $FA" | md5sum -c -
ln -sf $FA $C/t2t.fa
t0=$(date +%s)
if [ $SRC = tar ]; then
  cp $COHORT_TAR $W/cohort.tar; echo "S1 tar copied from Drive: $(stat -c%s $W/cohort.tar) B in $(( $(date +%s) - t0 )) s"
  want=$(cut -d' ' -f1 $COHORT_TAR.sha256); have=$(sha256sum $W/cohort.tar | cut -d' ' -f1)
  [ "$want" = "$have" ] || { echo "tar sha256 $have != $want"; false; }
  echo "S1 tar sha256 == $(basename $COHORT_TAR).sha256 ($have)"
  tar -xf $W/cohort.tar -C $C; rm -f $W/cohort.tar
else cp $COHORT_SRC/SHA256SUMS $C/; cp $COHORT_SRC/*.q4k.rr3 $C/; fi
echo "S1 cohort from $SRC: $(ls $C/*.q4k.rr3 | wc -l) archives, $(du -cb $C/*.q4k.rr3 | tail -1 | cut -f1) B on the VM disk, $(( $(date +%s) - t0 )) s"
t0=$(date +%s); ( cd $C && sha256sum -c --quiet SHA256SUMS ); echo "S1 sha256sum -c: all OK ($(( $(date +%s) - t0 )) s)"
awk -F'\t' 'NR>1{print $8"  "$1".q4k.rr3"}' $MANIFEST | sort > $W/man.sums; sort $C/SHA256SUMS > $W/drive.sums
cmp -s $W/man.sums $W/drive.sums || { echo "SHA256SUMS differ from the manifest's q4k_sha256"; diff $W/man.sums $W/drive.sums | head; false; }
echo "S1 SHA256SUMS == manifest ($(wc -l < $W/man.sums) archives)"

# ---------------------------------------------------------------- S2 cohort
STAGE=S2
python3 scripts/cohort558.py $C $MANIFEST q4k 2>&1 | tee $W/c558.txt
grep -q "C558 RESULT.PASS" $W/c558.txt

# ---------------------------------------------------------------- report
STAGE=report
{ echo "# panvram cohort558 (q4k) on the GPU - $(head -n 1 $W/gpu.txt); panvram $PV_COMMIT; $(date -u +%FT%TZ)"
  sed -n '2,$p' $W/gpu.txt
  echo; grep -h "^S0\|^S1" $RUNLOG
  echo; grep -vP "^C558 full\t" $W/c558.txt | grep "^C558"
  echo; echo "## per assembly (GPU full decode, FASTA XXH3 vs manifest)"; grep -P "^C558 full\t" $W/c558.txt
} > $OUT
echo "written: $OUT"
echo "DONE — выключи runtime"
