#!/usr/bin/env bash
# run_colab_panvram.sh - panvram on a Colab GPU (sm_80 or newer): build, fast tests, the public gate, throughput.
# panvram has no remote repository: upload the tarball made on ace-core (git archive, the commit inside) to
#   MyDrive/panvram/panvram.tar.gz
# The cohort for the gate: T2T and four HPRC year-1 haplotypes (HG00438.1/.2, HG00621.1/.2) from the Drive store
# run_colab_r3.sh uses (MyDrive/aceapex_corpus: t2t.fa[.gz], hprc/<name>.fa.gz; else fetched: T2T NCBI + md5, HPRC S3
# + sha256 of the year-1 index). The refrel3 v1 archives (q4k, q16k) come from ace-core via Drive:
#   MyDrive/panvram/cohort/<name>.<q4k|q16k>.rr3 + SHA256SUMS (== MANIFEST.tsv), checked with sha256sum -c;
# if that directory is absent they are encoded on the VM by aceapex research/refrel/refrel3v1.cpp at the pinned commit,
# built with the flags of MANIFEST.tsv (gcc -O3 -march=x86-64-v3 -funroll-loops: FMA contraction on; the encoder's
# match scores are doubles, so other flags / compilers can give other archive bytes - every archive decodes == FASTA
# either way); sha256 of each archive printed next to MANIFEST.tsv's.
#   S0 build: pip install -e . (torch cpp_extension, the GPU's arch), panvram.with_cuda must be true
#   S1 tests/test_synth.py: synthetic fixtures, CPU and CUDA paths, CUDA == CPU, refusals
#   S2 gate (tests/test_gate.py) for q4k and q16k: sample(1024, 8192) on the GPU == the windows of the FASTA rebuilt by
#      the CPU decoder (== the source files), == the CPU decoder's windows; fetch at 1000 coordinates == FASTA
#   S3 scripts/bench.py q4k and q16k: windows/s by W and n, fetch latency, resident bytes
# Result: MyDrive/aceapex_logs/panvram_<date>.txt; any failure: ..._FAILED.txt with the stage, line and command.
# Env: DRIVE, W, PV_TAR, COHORT_SRC, ACEAPEX_COMMIT, DATASETS ("q4k q16k")
set -Eeuo pipefail
ulimit -c 0   # no core dumps from test runs
shopt -s nullglob
DRIVE=${DRIVE:-/content/drive/MyDrive}
STORE=$DRIVE/aceapex_corpus; LOGS=$DRIVE/aceapex_logs
PV_TAR=${PV_TAR:-$DRIVE/panvram/panvram.tar.gz}
COHORT_SRC=${COHORT_SRC:-$DRIVE/panvram/cohort}
ACEAPEX_COMMIT=${ACEAPEX_COMMIT:-5b6d5cec0f5962a561ac48822a1b5c48793a5b47}
DATASETS=(${DATASETS:-q4k q16k})
W=${W:-/content/panvram_run}; PV=$W/panvram; ACE=$W/aceapex; D=$W/fa; C=$W/cohort; B=$W/bin
NAMES=(HG00438.1 HG00438.2 HG00621.1 HG00621.2)
[ -d $DRIVE ] || { echo "Drive not mounted: from google.colab import drive; drive.mount('/content/drive')"; echo "DONE — выключи runtime"; exit 1; }
[ -s $PV_TAR ] || { echo "no $PV_TAR: upload the panvram tarball from ace-core there"; echo "DONE — выключи runtime"; exit 1; }
mkdir -p $W $D $C $B $STORE/hprc $LOGS
DAY=$(date -u +%Y-%m-%d); RUNLOG=$W/run-$DAY.log; OUT=$LOGS/panvram_$DAY.txt; STAGE=setup
exec > >(tee -a $RUNLOG) 2>&1
on_err() { local rc=$? ln=$1 cmd=$2
  echo "FAILED in $STAGE: exit $rc at line $ln: $cmd"
  { echo "# panvram FAILED $(date -u +%FT%TZ) in stage $STAGE: exit $rc at line $ln: $cmd"; tail -n 120 $RUNLOG; } > $LOGS/panvram_${DAY}_FAILED.txt 2>/dev/null || true
  echo "details: $LOGS/panvram_${DAY}_FAILED.txt"; echo "DONE — выключи runtime"; }
trap 'on_err $LINENO "$BASH_COMMAND"' ERR
T=$(nproc)
nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv,noheader | tee $W/gpu.txt
command -v nvcc >/dev/null || export PATH=/usr/local/cuda/bin:$PATH
nvcc --version | tail -n 2 | tee -a $W/gpu.txt
python3 -c "import torch; print('torch', torch.__version__, 'cuda', torch.version.cuda, 'arch list', torch.cuda.get_arch_list())" | tee -a $W/gpu.txt
dpkg -s libzstd-dev >/dev/null 2>&1 || apt-get install -y -qq libzstd-dev >/dev/null

# ---------------------------------------------------------------- code
rm -rf $PV; mkdir -p $PV; tar -xzf $PV_TAR -C $W
PV_COMMIT=$( (gzip -dc $PV_TAR 2>/dev/null || true) | git get-tar-commit-id 2>/dev/null || true); PV_COMMIT=${PV_COMMIT:-unknown}; echo "panvram $PV_COMMIT"
if [ -d $ACE/.git ]; then git -C $ACE fetch -q origin refrel; else git clone -q -b refrel https://github.com/yasha1971-coder/aceapex.git $ACE; fi
git -C $ACE checkout -q $ACEAPEX_COMMIT; echo "aceapex $(git -C $ACE rev-parse HEAD) (encoder)"
( cd $ACE && g++ -std=c++17 -O3 -march=x86-64-v3 -funroll-loops -Isrc -Iresearch/refrel research/refrel/refrel3v1.cpp src/aceapex_api.cpp -lzstd -lpthread -o $B/refrel3v1 2> $W/enc_build.txt )

# ---------------------------------------------------------------- data
FA=$W/t2t.fa
if [ ! -s $FA ]; then
  if [ -s $STORE/t2t.fa ]; then cp $STORE/t2t.fa $FA
  elif [ -s $STORE/t2t.fa.gz ]; then gunzip -c $STORE/t2t.fa.gz > $FA
  else curl -fsSL --retry 3 -o $W/t2t.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz
    gunzip -c $W/t2t.fa.gz > $FA; cp $W/t2t.fa.gz $STORE/t2t.fa.gz; rm -f $W/t2t.fa.gz; fi
fi
echo "cd1e52ce400c027ed0b7ab4b9d613f5a  $FA" | md5sum -c -
HD=$STORE/hprc; IDX=$HD/index.tsv
[ -s $IDX ] || curl -fsSL --retry 3 -o $IDX https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index
for nm in "${NAMES[@]}"; do
  smp=${nm%.*}; hap=${nm##*.}
  read -r url sha < <(awk -F'\t' -v s=$smp -v h=$hap '$1==s{ if (h==1) print $2, $6; else print $3, $7 }' $IDX)
  [ -n "$url" ] || { echo "$nm not in the year-1 index"; false; }
  u=${url/s3:\/\/human-pangenomics/https:\/\/s3-us-west-2.amazonaws.com\/human-pangenomics}
  if [ ! -s $HD/$nm.fa.gz ]; then curl -fsSL --retry 3 -o $HD/$nm.fa.gz.part "$u"; mv $HD/$nm.fa.gz.part $HD/$nm.fa.gz; fi
  [ "$(sha256sum $HD/$nm.fa.gz | cut -d' ' -f1)" = "$sha" ] || { echo "$nm: sha256 mismatch"; false; }
  [ -s $D/$nm.fa ] || gunzip -c $HD/$nm.fa.gz > $D/$nm.fa
done
echo "HPRC: ${NAMES[*]} (sha256 == the year-1 index)"
ln -sf $FA $C/t2t.fa
STAGE=archives
gcc --version | head -n 1 | tee -a $W/gpu.txt
if [ -s $COHORT_SRC/SHA256SUMS ]; then
  for ds in "${DATASETS[@]}"; do for nm in "${NAMES[@]}"; do cp $COHORT_SRC/$nm.$ds.rr3 $C/; done; done
  ( cd $C && grep -E "\.($(echo ${DATASETS[@]} | tr ' ' '|'))\.rr3$" $COHORT_SRC/SHA256SUMS | sha256sum -c - ) | tee -a $W/encode.txt
  echo "archives from $COHORT_SRC (sha256 checked)" | tee -a $W/encode.txt
fi
for ds in "${DATASETS[@]}"; do q=$([ $ds = q4k ] && echo 4096 || echo 16384)
  for nm in "${NAMES[@]}"; do [ -s $C/$nm.$ds.rr3 ] || $B/refrel3v1 encode $FA $q $T $D/$nm.fa $C/$nm.$ds.rr3 | tee -a $W/encode.txt
    m=$(awk -F'\t' -v n="y1_$nm" -v d=$ds '$1==n{ print (d=="q4k") ? $8 : $10 }' $PV/MANIFEST.tsv 2>/dev/null || true)
    echo "ARCHIVE $ds $nm $(stat -c%s $C/$nm.$ds.rr3) B sha256 $(sha256sum $C/$nm.$ds.rr3 | cut -d' ' -f1) MANIFEST ${m:--}" | tee -a $W/encode.txt; done; done

# ---------------------------------------------------------------- S0 build
STAGE=S0; cd $PV
pip uninstall -y -q panvram 2>/dev/null || true
t0=$(date +%s); MAX_JOBS=$T pip install -e . > $W/build.txt 2>&1 || { tail -40 $W/build.txt; false; }
echo "S0 pip install -e . : $(( $(date +%s) - t0 )) s"
python3 -c "import panvram, torch; assert panvram.with_cuda, 'built without CUDA'; print('S0 panvram', panvram.__version__, 'with_cuda', panvram.with_cuda, 'device', torch.cuda.get_device_name(0), 'sm_%d%d' % torch.cuda.get_device_capability(0))"
grep -E "warning" $W/build.txt | grep -v -E "torch/include|pybind11" | head -20 || true

# ---------------------------------------------------------------- S1 fast tests
STAGE=S1
python3 -m pytest -q -rs tests/test_synth.py 2>&1 | tee $W/s1.txt
grep -qE "passed" $W/s1.txt
if grep -qE "failed|error|skipped" $W/s1.txt; then echo "S1: not every test passed"; false; fi

# ---------------------------------------------------------------- S2 gate
for ds in "${DATASETS[@]}"; do STAGE="S2 gate $ds"
  PANVRAM_COHORT=$C PANVRAM_DATASET=$ds PANVRAM_DEVICE=cuda PANVRAM_SOURCES=$D python3 -m pytest -q -s tests/test_gate.py 2>&1 | tee $W/gate_$ds.txt
  grep -q "GATE RESULT.PASS" $W/gate_$ds.txt; done

# ---------------------------------------------------------------- S3 throughput
for ds in "${DATASETS[@]}"; do STAGE="S3 bench $ds"
  python3 scripts/bench.py $C $ds cuda 2>&1 | tee $W/bench_$ds.txt; grep -q "BENCH done" $W/bench_$ds.txt; done

# ---------------------------------------------------------------- report
STAGE=report
{ echo "# panvram on the GPU - $(head -n 1 $W/gpu.txt); panvram $PV_COMMIT; encoder aceapex $ACEAPEX_COMMIT; $(date -u +%FT%TZ); $T host threads"
  sed -n '2,$p' $W/gpu.txt
  echo; echo "## archives (encoded on the VM)"; grep -h "^ARCHIVE" $W/encode.txt
  echo; echo "## S0 build"; grep -h "^S0" $RUNLOG
  echo; echo "## S1 tests/test_synth.py"; tail -n 3 $W/s1.txt
  echo; echo "## S2 gate"; for ds in "${DATASETS[@]}"; do echo "$ds:"; grep -h "^GATE" $W/gate_$ds.txt; done
  echo; echo "## S3 throughput"; for ds in "${DATASETS[@]}"; do grep -h "^BENCH" $W/bench_$ds.txt; done
} > $OUT
echo "written: $OUT"
echo "DONE — выключи runtime"
