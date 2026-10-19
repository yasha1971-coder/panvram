#!/usr/bin/env bash
# Cas-OFFinder over cohort assemblies, one at a time: the source FASTA is decompressed, searched and deleted, because
# 50 assemblies uncompressed would need about 150 GB and only 22 GB are free. Times per assembly to casoff_times.txt.
set -euo pipefail
cd "$HOME/pubrepo/.wk/gs"
export OCL_ICD_VENDORS=$PWD/icd LD_LIBRARY_PATH=$PWD/OpenCL-ICD-Loader/build
NAMES=${1:?}; OUT=${2:?}; T=asmrun
mkdir -p $T/genome; : > "$OUT"; : > casoff_times.txt
printf '%s\n' "$T/genome" 'NNNNNNNNNNNNNNNNNNNNNGG' 'GAGTCCGAGCAGAAGAAGAANNN 6' 'GGGTGGGGGGAGTTTGCTCCNNN 6' 'GACCCCCTCCACCCCGCCTCNNN 6' 'GGTGAGTGAGTGTGTGCGTGNNN 6' > $T/in.txt
while read -r nm; do
  base=${nm#y1_}
  rm -f $T/genome/*.fa
  s=$SECONDS; zcat "$HOME/pubrepo/.wk/hprc50/$base.fa.gz" > "$T/genome/$base.fa"; gz=$((SECONDS-s))
  s=$SECONDS; ./cas-offinder/build/cas-offinder $T/in.txt C $T/out.txt > /dev/null 2>&1; se=$((SECONDS-s))
  awk -v a="$nm" 'BEGIN{FS=OFS="\t"} !/^#/{print a,$0}' $T/out.txt >> "$OUT"
  n=$(grep -vc '^#' $T/out.txt || true)
  echo -e "$nm\tgunzip ${gz}s\tsearch ${se}s\tsites $n" | tee -a casoff_times.txt
  rm -f "$T/genome/$base.fa" $T/out.txt
done < "$NAMES"
rm -rf $T
