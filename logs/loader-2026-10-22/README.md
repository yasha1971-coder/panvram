# A window loader for DNA model training: panvram against GenVarLoader and pyfaidx, CPU (2026-10-22)

Training a sequence model on a pangenome means pulling batches of random fixed-length windows. This directory measures
how fast three loaders deliver such batches on one machine, and what each of them actually returns.

**The three are not interchangeable.** panvram returns **the assembly's own bases**, decoded from its refrel3 v1
archive. GenVarLoader returns **the reference with that sample's variants applied**, at reference coordinates. pyfaidx
returns the assembly's bases from a plain indexed FASTA. So **only the speed is compared here**; the contents of A and
B differ by construction and no agreement between them is claimed or measured. A's contents are checked against the
source FASTA (below); B's are not compared with A at all.

Run by the author on ace-core (AMD EPYC 4344P, 8 cores / 16 threads, 125 GB, Ubuntu 22.04.5, kernel 5.15.0-163,
g++ 11.4.0, Python 3.10.12), behind a silence gate (one-minute load average below 0.5) before every row. No GPU is
present on this host (`nvcc` absent), so every row is the CPU path; the GPU path of A is the notebook in this
directory, **not run**.

Protocol, identical for all three: batches of **1 024 windows**, W = **8 192** and **32 768** bases, 3 warm-up batches
then **9 timed batches**, median reported, output a `uint8 [1024, W]` array. A and C answer **the same requests**: the
coordinates of A's first timed batch are written to `requests/` and replayed by C. B cannot answer those coordinates
(it works in reference space), so it draws its own windows of the same size and count.

## Result

| loader | scope | 8 KiB windows/s | 8 KiB GB/s | 32 KiB windows/s | 32 KiB GB/s | threads | open | peak RSS |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **A** panvram v1 q16k | 5 assemblies | 438 468 | 3.592 | 155 929 | 5.109 | 16 | 8.1 s | 3.5 GB |
| A | 5 assemblies | 94 099 | 0.771 | 31 495 | 1.032 | 1 | 8.0 s | - |
| **A** panvram v1 q16k | **558 assemblies** | 411 668 | 3.372 | 180 609 | 5.918 | 16 | 10.8 s | 20.7 GB |
| A | 558 assemblies | 80 011 | 0.655 | 31 847 | 1.044 | 1 | 16.2 s | - |
| **B** GenVarLoader 0.43.1 | chr1, 3 samples | **2 337 735** | **19.151** | 211 398 | 6.927 | 16 | 0.2 s | 1.4 GB |
| B | chr1, 3 samples | 1 259 337 | 10.316 | 77 266 | 2.532 | 1 | 0.2 s | - |
| **C** pyfaidx 0.9.1, bgzip FASTA | 5 assemblies | 9 073 | 0.074 | 3 951 | 0.129 | 1 | 0.2 s | 0.36 GB |

**GenVarLoader is 5.3 times faster than panvram at 8 KiB** (2 337 735 against 438 468 windows/s, both with every
thread available). At 32 KiB the two are close: 211 398 against 180 609, a factor of 1.17. Against pyfaidx, panvram is
48 times faster at 8 KiB and 39 times at 32 KiB.

Read those numbers with the scopes in the table: B covers **one chromosome (248 Mb) and 3 samples**, a working set that
stays in RAM and is re-read every batch; A covers 15 Gb of bases over 5 assemblies, or 1.68 Tb over 558. Raw logs and
the per-row JSON are in `cpu/`.

## The GPU path of A: measured, not yet publishable

One Colab run of `panvram_loader_A100_v3.ipynb`, 2026-10-25, **NVIDIA A100-SXM4-40GB**, 558 assemblies resident in
**11.567 GB** of VRAM (`nvidia-smi` 11 466 MiB), opened in **24.499 s**; same batch of 1 024, 9 timed runs, median;
`CHECK W 8192 windows 256 identical 256 differ 0`. Raw output `gpu/v3_run.log`, with derived figures
`gpu/v3_result.json`.

| W | ms per batch | min / max | windows/s | GB/s | blocks decoded | us per decoded block |
|---:|---:|---|---:|---:|---:|---:|
| 8 192 | 3.872 | 3.693 / 4.007 | 264 475 | 2.167 | 1 536 | **2.521** |
| 32 768 | 2.575 | 2.501 / 2.766 | 397 729 | 13.033 | 3 072 | **0.838** |

**Open question, and the reason no GPU row joins the table above.** Each warp decodes exactly one 16 KiB block, so
microseconds per decoded block must not depend on W; here it differs by 3.0x. The block counts come from the launch
formula `bpw = (W + Q - 1) / Q + 1` (`csrc/pv_cuda.cu:128`) and from how many blocks a window touches on average,
`1 + (W-1)/Q`: 1.50 at 8 KiB and 3.00 at 32 KiB. So the 8 KiB batch decodes **half** as many blocks and must be the
cheaper one - with an equal cost per block the 32 KiB batch should have taken 7.744 ms, and it took 2.575 ms. The
standing hypothesis is the order of the sweep: 8 KiB was measured first, right after ~14 GB were pushed to the card
and after three 8 MB warm-up batches, far too little to bring the clocks up. The run recorded no clocks, so the cause
is **NOT MEASURED**. `gpu/README.md` states it in full, and `panvram_loader_A100_v4.ipynb` shuffles the W order,
repeats the sweep three times, ramps the card first and prints microseconds per decoded block and the clocks.

A property that is **not** an artefact, and shows in the CPU rows too: at q16k a 8 KiB window decodes 24 575 bases to
deliver 8 192 (3.00x) while a 32 KiB window decodes 49 151 to deliver 32 768 (1.50x), so the small window costs twice
the decode per delivered base. That is why GB/s rises with W everywhere in this directory.

## What each loader does per window, and where its time goes

- **A, panvram v1 q16k.** Decodes the rANS token stream of every 16 KiB block the window touches and executes the copy
  operations against the decoded reference. It is the only one of the three that decompresses a per-assembly
  representation. 16 threads give 4.7 times one thread (438 468 against 94 099 at 8 KiB), so it scales with cores and is
  bound by work, not by memory. Going from 8 to 32 KiB raises its rate from 3.59 to 5.11 GB/s: the fixed per-request
  cost is amortised over more bases. Growing the cohort from 5 to 558 assemblies costs 6 % at 8 KiB and nothing at
  32 KiB, while peak RSS goes from 3.5 to 20.7 GB.
- **B, GenVarLoader.** Copies the reference slice and applies the few variants inside it. The normalised chr1 VCF holds
  590 710 records for the three samples, so an 8 KiB window meets on the order of ten; the work is memcpy from an
  uncompressed in-RAM reference plus a handful of edits, which is why 19 GB/s is reachable. That the variants really
  are applied was checked, not assumed: over 64 windows, haplotype 1 differs from the reference track in 110 965 bases
  and haplotype 2 in 132 757 of 524 288, and 50 of the 64 windows carry at least one difference.
- **C, pyfaidx on bgzip FASTA.** Inflates the bgzip blocks of the window on every request, single-threaded, and builds
  a Python string per window. 0.074 GB/s. This is the cost of the usual tooling around an indexed FASTA, not a limit of
  the file format.

## Correctness of A

`scripts/check_a_vs_fasta.py`: every window of both request sets, decoded by panvram, compared byte for byte with the
same coordinates read from the source FASTA (bgzip copies of the very files the archives were built from), case kept.

| W | windows | identical | differ |
|---:|---:|---:|---:|
| 8 192 | 1 024 | **1 024** | 0 |
| 32 768 | 1 024 | **1 024** | 0 |

B's contents are **not** checked against A or against any assembly: reference-plus-variants and the assembly's own
bases are different objects. Structural variants, large insertions and the 14.9 Gb of non-ACGT runs of these
assemblies are not represented in B at all.

## How B was scoped, and why

GenVarLoader requires left-aligned, bi-allelic, atomized variants. The HPRC v1.1 minigraph-cactus VCF against CHM13
(`hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz`, 3 021 114 234 B, from the HPRC release of BioProject PRJNA730823) is
multi-allelic, so it was subset to **chr1 and the three samples whose haplotypes make up the five assemblies of A**
(HG00438, HG00621, HG00673), then `bcftools norm -m -any -f <reference>` and `bcftools sort`: 19 441 645 B,
**590 710 records**. `bcftools` is not installed on this host; 1.21 was built from source for this run. The whole
genome did not fit the time budget of this measurement. The reference for B is the same T2T FASTA with its contigs
renamed to `chr1..chrM` to match the VCF (`t2t.chr.fa.gz`, SHA-256 in `inputs.sha256`).

C needed its own preparation: the assembly FASTA on this host are plain gzip without an index, which pyfaidx cannot
read, so the five were re-compressed with `bgzip` and indexed with `samtools faidx`.

## Files

- `cpu/` - one JSON and one log per row: `A5`, `A558`, `A5_1t`, `A558_1t` (panvram), `B5`, `B5_1cpu` (GenVarLoader),
  `C5` (pyfaidx). The logs keep the result lines, the wall time and the peak RSS of `/usr/bin/time -v`.
- `requests/A5.coords.W8192.json`, `requests/A5.coords.W32768.json` - the 1 024 windows of each request set as
  (assembly index, contig index, start in the contig); the assembly order is the one `cpu/A5.json` records.
- `scripts/bench_a.py`, `bench_b.py`, `bench_c.py`, `check_a_vs_fasta.py` - exactly what was run. `PV_THREADS` selects
  A's thread count.
- `inputs.sha256` - the reference, the renamed reference, both VCF, the five bgzip FASTA and the five v1 archives.
- `panvram_loader_A100_v4.ipynb` - the GPU path of A for Colab (A100): same protocol, same batch and window sizes,
  the result tensor already on the card, plus a check of 256 windows against the CPU decoder. **The only GPU numbers
  so far are the two times above, from v3; see `gpu/README.md`.** `N_ASM` is 558 by default; the cohort is fetched and verified by
  `scripts/colab_prepare.py`, the same file the notebook embeds, and any 0-byte archive or wrong SHA-256 stops the run.
  **Its five-assembly mode is not the five assemblies of the CPU rows:** the record ships the archives only inside a
  tar ordered by name, so the first five complete ones are release-2 haplotypes (`r2_HG00097_hap1`, `r2_HG00097_hap2`,
  `r2_HG00099_hap1`, `r2_HG00099_hap2`, `r2_HG00126_hap1`), while the CPU rows used five year-1 haplotypes. That mode
  controls cohort size, not the names.
- `scripts/colab_prepare.py` - the fetch and verification used by the notebook, tested on this host in both modes
  (below).
- `gpu/README.md` - what the GPU run reported, why it is not a window-size result, and what v4 changes.
- `panvram_lane0_A100_v2.ipynb` - a separate study: what share of a block's time lane 0 spends decoding the rANS
  stream while lanes 1..31 wait, by Nsight Compute on the shipped kernel and by `clock64()` in an instrumented build
  (`scripts/lane0_patch.py`). **No numbers yet.**
- The superseded notebooks are kept as a record, each with the reason on its first cell:
  `panvram_loader_A100_OLD_do_not_use.ipynb` fetched the record's individual `*.q16k.rr3` files, which are not
  addressable there (HTTP 404), so it left 0-byte archives; `..._v2_OLD_...` checked the reference against the SHA-256
  of a different rendering of the same genome; `..._v3_OLD_...` swept the window sizes in a fixed order;
  `panvram_lane0_A100_OLD_do_not_use.ipynb` rebuilt the extension and tried to `importlib.reload` it in a kernel that
  had already loaded the old shared object.

### What was tested on this host before the notebook was published

`scripts/colab_prepare.py` was run in both modes, up to the point where a GPU is required:

| mode | fetch | verification | cohort opened (CPU) |
|---|---|---|---|
| `N_ASM = 5` | 104 857 600 B of `part00` by HTTP range; 7 tar members seen, 5 archives complete | 5 / 5 SHA-256 equal `SHA256SUMS.q16k`, no empty file | 5 assemblies, 8.05 s; batches of 1 024 drawn at both W |
| `N_ASM = 558` | all 16 parts, each against `PARTS.sha256`; joined tar 8 098 467 840 B, SHA-256 equals the record's | 558 / 558 SHA-256 equal the `SHA256SUMS` inside the tar, no empty file | 558 assemblies, 10.74 s; batches of 1 024 drawn at both W |

The two stop conditions were tested by breaking a verified cohort on purpose: an archive truncated to 0 bytes gives
`r2_HG00097_hap1_...q16k.rr3 is 0 bytes - the download did not deliver the file`, and one extra byte appended gives
`SHA-256 6ba0a939... does not match the record's 43ec9ec6...`. Both raise before anything is opened.

Versions: panvram 1.0.0 (CPU extension), torch 2.14.1+cpu, genvarloader 0.43.1, genoray 6.0.2, polars 1.44.2,
pyfaidx 0.9.1, biopython 1.88, htslib/samtools 1.13, bcftools 1.21.

HPRC assemblies and variants: BioProject PRJNA730823. Their use here implies no endorsement by the HPRC.
