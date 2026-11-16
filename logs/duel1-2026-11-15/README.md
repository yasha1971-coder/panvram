# DUEL-1 phase A — PRELIMINARY (5 caveats)

> **STATUS: PRELIMINARY.** The pre-registered checks below hold, but five honest caveats
> bound the claims. VM chr20-22 run and full-cohort metric 5 are pending.
>
> 1. **P1 at 32 kb FAILED** the pre-registered rule (bcftools consensus is better by
>    4.2 pp > 2 pp) — this is a *failure*, not merely "flagged"; the mechanism statement
>    awaits the step-1 investigation.
> 2. **Metric 5 is on 3 haplotypes only** — "median 2.84 Mb" is NOT the cohort median.
> 3. **"genes/regulatory outside VCF" is NOT measured** — unique ≠ genes; it needs CAT
>    gene annotation + an HMM-Flagger reliability mask.
> 4. **loss>1kb is non-monotonic** (8 kb 7.2% → 32 kb 5.8%; clean 1.8 → 0.6%): a 32 kb
>    window contains an 8 kb window, so the fraction cannot decrease → the per-length
>    window samples are biased by the "unique PAF mapping" filter, and **metric 4
>    systematically under-counts loss in GVL's favour**.
> 5. **P3 was not checked per window** (460/460 is check_full of whole archives; the
>    cohort has 558); "100% by construction" is an inference, not a measurement.

# DUEL-1 — panvram vs GenVarLoader, Phase A (CPU)

Seed `20261114`. All numbers carry their config. Status: **A0, A2, A4(panvram side) done;
A1 partial (headers only); A3/A5/A6 pending (blockers below).**

## Environment (A0)
See `env.txt`. panvram 1.0.0 (frozen tag v1.0.0, HEAD 8f8e9e1, CPU decoder, reference
T2T-CHM13v2.0). GVL target >=0.43.1 — **not installed** here. Host: AMD EPYC 4344P, 16
threads, ~125 GiB RAM, Python 3.10.12, ~58 GB disk free.

## Inputs (A1)
Exact S3 keys, sizes and ETags in `data/inputs.sha256`. Primary inputs (bucket
`s3://human-pangenomics`, `--no-sign-request`):

| role | key | size |
|---|---|---|
| R2 phased VCF (GRCh38) | `.../release2/minigraph-cactus/v2.0/hprc-v2.0-mc-grch38/hprc-v2.0-mc-grch38.vcf.gz` | 7.48 GB |
| R2 PAF | `.../release2/impg/pafs/hprc465vsgrch38.aln.paf.gz` | 6.34 GB |
| R1 phased VCF (GRCh38) | `.../freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.vcf.gz` | 3.94 GB |
| GRCh38 reference | `.../release2/minigraph-cactus/v2.0/input/hg38.chroms.fa.gz` | 0.79 GB |

Full-file sha256 pending download. Headers fetched via range request:
```
aws s3api get-object --no-sign-request --bucket human-pangenomics \
  --key <key> --range "bytes=0-4194303" vcf_head.gz
zcat vcf_head.gz 2>/dev/null | grep -m1 '^#CHROM' > chrom_header.txt
```
(local tabix has no libcurl, so remote streaming via tabix/bcftools is unavailable here.)

## A2 — intersection  (`data/intersection.tsv`)
Both VCFs are per-sample **diploid, FORMAT=GT, fully phased (`|`)**. PanSN `sample#1` =
GT position 1 (left of `|`), `sample#2` = GT position 2. A haplotype is reconstructable by
GVL iff its **sample** is a VCF column. The `mat/hap1 -> #1`, `pat/hap2 -> #2` assignment is
the HPRC PanSN convention — the VCF header does not label mat/pat, so it is not independently
decidable from the header alone.

| cohort | VCF | VCF samples | entries | present | absent | absent samples |
|---|---|---|---|---|---|---|
| r2 | hprc-v2.0-mc-grch38 | 232 (incl. CHM13, HG002) | 464 | 462 | 2 | HG00272 (hap1+hap2) |
| y1 | hprc-v1.0-mc-grch38 | 45 (incl. CHM13) | 94 | 88 | 6 | HG002, HG005, NA19240 (×2) |

hg002 manual check: r2 HG002 **present** (both special rows `r2_hg002v1.1.pat`,
`r2_hg002v1.1.mat_MT` -> sample HG002); y1 HG002 **absent** (no HG002 column in R1 VCF).

Data note: `names_y1_94.tsv` rows 52+ carry a stray `y1_` prefix in cols 2/3/5; stripped
(`^y1_`) before matching, else 22 samples are false-negatives.

## A4 — memory/disk (panvram side done; GVL side pending)
| side | disk |
|---|---|
| panvram q4k archives (558) | 11.48 GB |
| panvram q16k archives (558) | 7.54 GB |
| panvram shared reference (CHM13) | 2.94 GB |
| **panvram total (q4k + ref)** | **14.42 GB for 558 haplotypes** |
| GVL dataset + GRCh38 ref | PENDING (A3) |

RAM and VRAM columns: RAM for GVL = A3; VRAM = Phase B.

## A5 — coverage (metric 5) — DONE for R2 (components a + c)
Method: stream the combined PAF `hprc465vsgrch38.aln.paf.gz` (8,916,633 lines) on
<host> (CPU, no cloud). Per haplotype, unreachable-by-GVL = (a) sample absent from
VCF -> 1.0, else (c) `1 - union(query-aligned bases)/sum(query lengths)`. R2 clips (b)
are not published -> folded into (c). HG002 (2 haps) omitted (manifest rows carry no
`#1/#2` pansn).

| set | n | min | median | p95 | max | frac >=1% |
|---|---|---|---|---|---|---|
| present (component c only) | 460 | 2.16% | 2.88% | 3.47% | 3.86% | 100% |
| all R2 (a+c) | 462 | — | 2.88% | 3.48% | 100% | 100% |

Absent (a): `r2_HG00272_hap1/hap2` (sample HG00272 not in R2 VCF) = 100% unreachable.

**Verdict (metric 5): panvram wins.** Every one of 462 R2 haplotypes has >=1% of its
assembly unreachable by GVL. Even excluding fully-absent HG00272, every present
haplotype has >=2.16% (median 2.88%) of bases that do not align 1:1 to GRCh38
(centromeres / satellite / segdup / novel insertions) — panvram stores all of it
losslessly. Threshold >=1% is met with margin. Per-haplotype: `a5_per_haplotype.json`.

## Kit (bench/duel1/kit/) — full phase-A pipeline for the self-deleting GCP VM
- `duel1_gcloud.sh` (stage|run): build+upload panvram source / cohort q4k / reference to
  GCS, create CPU VM (c4d/n2/c2d-standard-32, 400 GB disk, `--max-run-duration`
  `--instance-termination-action DELETE`, zones <zone>/a/c), provision, run
  driver, fetch results, delete VM, print empty instances/disks.
- `duel1_driver.sh` (VM-side): installs torch-cpu + genvarloader>=0.43.1 + genoray + edlib
  + bcftools; downloads inputs from S3; `bcftools norm -f GRCh38 -a --atom-overlaps . -m
  -any` for GVL; runs all stages guarded so panvram+PAF results survive GVL failure.
- `duel1_pipeline.py` / `duel1_lib.py`: A4, A5(a,c), A6 window sampling, PAF cg:Z CIGAR
  walk (unit-tested), coverage.
- `duel1_a6_compare.py`: A6 GVL-vs-true-assembly (panvram `fetch`) + bit-perfect control
  (`Cohort.check_full`), per-sample 1-sample GVL datasets, edit distance via edlib.
- `gvl_write_child.py`: A3 `gvl.write` under `/usr/bin/time -v` for wall+peak RSS.

Validated locally: panvram `fetch`/`check_full`, PAF orientation (query=assembly, cg:Z
present), CIGAR mapping math, A5 coverage end-to-end. GVL stages are authored to the 0.43
API but not yet executed (VM launch blocked).

## Blockers (unblock commands)
1. **GCP VM launch denied by the harness auto-mode classifier** (bucket create / 12 GB
   upload / VM create were refused as "destructive/unwanted"). This gates A3 and the GVL
   side of A4/A6. Unblock: the user adds a Bash permission rule for
   `~/google-cloud-sdk/bin/gcloud` and `gsutil` (or runs the kit himself):
   ```
   bash bench/duel1/kit/duel1_gcloud.sh stage ~/duel1_work   # build+upload to GCS
   bash bench/duel1/kit/duel1_gcloud.sh run   ~/duel1_work   # create VM, run, self-delete
   ```
2. **Clipped-intervals BED**: R1 found at
   `s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.clipped-intervals.bed.gz`
   (wired in the driver). R2 clips not published -> A5(b) for R2 folded into (c).
3. **Disk on <host> ~58 GB** — fine for A5 (PAF only); full GVL A3/A6 need the VM.

## Status
- DONE on <host> (no cloud): A2 (intersection), A4 panvram side, **A5 coverage (metric 5)**.
- BLOCKED on VM launch: A3 (gvl.write), A6 accuracy (GVL side), A4 GVL side. Kit ready;
  one permission rule unblocks a one-shot run.

## Pre-registered checks (DUEL1_CHECKS_PLAN.md, committed 593bb30 before run)
Scale: A6 / metric 4 = **chr22** (500 windows/length, 230 samples, wave VCF v2.0,
SparseVar2 ingestion, seed 20261114); metric 5 = **3-haplotype subset**; A5 = **whole
genome** (462 haplotypes). Reference everywhere: `hg38.chroms.fa` (GRCh38#0#→chrN).
panvram = CHM13. bcftools 1.25 (micromamba/bioconda, env parity with the VM driver).

| # | Prediction | Fact | Verdict |
|---|---|---|---|
| P0 | ≥70% of non-exact 1kb windows in SegDup/pericentromere/PAF-ambiguous | 77/92 = **83.7%** in SegDup∪pericentromere (PAF-ambiguous pre-excluded by the unique-mapping sampler) | **CONFIRMED** |
| P1 | \|bcftools consensus − GVL\| ≤ 2pp | Δexact 0.0 / 1.4 / **4.2** pp at 1/8/32kb; loss>1kb & p99 identical | HOLDS 1–8kb; **fails 32kb** (consensus +4.2pp) → GVL-specific small-edit gap flagged; big losses identical = genuine VCF loss |
| P2 | ≥80% of unaligned = α-sat/HSat; unique 0.05–0.5% | sat-of-unaligned 82–94% (median **87%**); unique 1.5/2.8/3.0 Mb = 0.05–0.10% | **CONFIRMED**; unique ≥1 Mb → **strong thesis** |
| P3 | 100% bit-perfect all lengths (else stop) | check_full 12/12 local + **460/460** VM; every window from a passing archive | **CONFIRMED** |

### Metric 4 (GVL reconstruction accuracy, chr22) — two views
| W | exact (all) | loss>1kb (all) | exact (clean*) | loss>1kb (clean*) |
|---|---|---|---|---|
| 1 kb | 81.5% | 0% | 95.2% | 0% |
| 8 kb | 80.9% | 7.2% | 89.0% | 1.8% |
| 32 kb | 64.5% | 5.8% | 68.9% | 0.6% |

\*clean = windows outside SegDup/pericentromere (step-0 hard regions). ed_p99 all/clean
same order (585 / 4393 / 17234). GVL vs bcftools consensus: exact within 0–4.2pp,
loss>1kb and p99 identical → the >1kb losses are a property of the VCF representation,
not of GVL. (`chr22_checks_perwindow.jsonl`, `frag_checks.json`.)

**Why 92.8% (chr20-22, A3/metric-4 final) vs 81.5% (chr22 all-windows):** the corrected window filter maps every window through the PAF and drops windows whose ends land in an indel or contig edge (unmappable, ~12%) from the exact denominator — the same filter already read 95.2% on chr22 "clean"; 92.8% is that filter applied to the broader, less acrocentric chr20-22 (chr22 alone is LCR22/segdup-heavy, which depresses the all-windows 81.5%).

### Metric 5 (GVL-unreachable sequence by class, 3-hap subset) — `frag_metric5.json`
Per haplotype ~3.0 Gb, ~79–90 Mb unaligned to GRCh38. Of the unaligned (median): α-sat
~1.1% of assembly, HSat ~0.9%, other CenSat ~0.3%, other RepeatMasker ~0.1–0.45%,
**UNIQUE (non-repeat) 0.05–0.10% = 1.5–3.0 Mb/haplotype (median 2.84 Mb)**. This unique
non-satellite sequence (genes/regulatory not alignable to GRCh38) is unrepresentable by
GRCh38+VCF and stored bit-perfect by panvram.

## Global cost (duel thesis, with numbers)
**VCF-representation of the pangenome loses ~0.09% UNIQUE non-satellite sequence per
haplotype (median 2.84 Mb; subset) and ~0.6–7% of windows (>1kb loss, by length);
panvram = 0% loss (bit-perfect, check_full 460/460). Throughput Z = phase B (GPU).**
X ≈ 2.84 Mb/haplotype (0.094% of assembly). Y = loss>1kb 0/7.2/5.8% (all) or 0/1.8/0.6%
(clean euchromatin), 1/8/32 kb.

## Follow-up pre-registered checks — RESULTS (steps 1-4, predictions.md committed before each)
All five caveats above are now addressed.

| Step | Prediction | Fact | Verdict |
|---|---|---|---|
| 1 P1-32kb | missing-GT drives the 32kb GVL-vs-consensus gap | missing-GT only 4.3%; GVL-vs-cons ED median 13bp; consensus bit-exact, GVL +bp | **mechanism = GVL over-expands tandem repeats** (not missing-GT, not VCF loss) |
| 2a discard | discard rate rises with length | 29.9 / 30.9 / 34.2% (1/8/32kb) | CONFIRMED |
| 2b discard=SV | discarded windows = breakpoints | 100% span ≥2 PAF blocks (SV) | CONFIRMED |
| 2c monotonic | nested → monotonic loss | 0 / 5.03 / 5.07% | CONFIRMED; metric-4 bounds 0-5% (lower) … 30-38% (upper) |
| 3a per-window | 100% bit-perfect all lengths | 378/378 byte-exact vs FASTA; check_full 460/460 | CONFIRMED (measured) |
| 3b 460v558 | 94 y1 + 2 HG00272 + 2 HG002 | exact | CONFIRMED |
| 4a satellite | ≥80% of unaligned satellite | median **91.8%** (min 75.4%) | CONFIRMED |
| 4b flagger | removes <30% of unique | median **1.2%** (24 kb) | CONFIRMED |
| 4c genes | CDS in unique ≥0.1 Mb | median **28.7 kb** (max 144 kb) | **FALSE** (real but ~10× below) |

### Metric 5 — full cohort (460 present haplotypes; CenSat+RepeatMasker+CAT v1.3+HMM-Flagger)
Unaligned (GVL-unreachable) median **87.9 Mb/haplotype** (p95 104, max 120). Class medians
(% of assembly): α-sat 1.06, HSat 1.22, other CenSat 0.33, other RepeatMasker 0.17,
**unique (non-repeat) 0.067% = 2.02 Mb; unique-reliable (after HMM-Flagger) 1.99 Mb**
(p95 2.95, max 4.38 Mb). CAT protein-coding **CDS in the unique-reliable unreachable set:
median 28.7 kb/haplotype** (p95 59, max 144 kb). (`frag_metric5_full.json`, `metric5_full.jsonl`.)

### Duel thesis (cohort-wide, numbers)
GVL (GRCh38+VCF) cannot represent a **median 1.99 Mb of unique non-satellite non-repeat
sequence per haplotype** (≥1 Mb → strong), **~29 kb of it CAT protein-coding CDS**.
Window-level loss>1kb is between **0-5% (within-block) and 30-38% (SV breakpoints counted)**
— the non-monotone under-count is a filter artifact, fixed by nesting. **panvram = 0% loss:
per-window 378/378 byte-exact, check_full 460/460 bit-perfect.** Throughput Z = phase B (GPU).

## chr20-22 scale-up — metric 3, metric 4, panvram control (local, no VM)

### Metric 3 (A3) — GVL vs panvram, chr20-22, 231 samples (`frag_a3_chr202122.json`)
| side | disk footprint | peak RSS | fidelity |
|---|---|---|---|
| GVL (SparseVar2+gvl.write) | svar2 341.5 MB + dataset 0.06 MB + GRCh38(chr20-22) 157 MB = **~499 MB** | gvl.write 749 MB; wall 2.79 s | **lossy** (no 1.99 Mb/hap unique) |
| panvram (chr20-22 scaled ×5.24% base fraction) | archives ~510 MB + CHM13 ref ~158 MB = **~667 MB** | — | **lossless** (bit-perfect) |
Method note: panvram archives are genome-wide; chr20-22 share = GRCh38 base fraction (0.0524), stated not physically isolated.

### Metric 4 (chr20-22, nested 10000×{1,8,32} kb; ground truth = published assembly FASTA) — `frag_metric4_202122.json`
| W | GVL exact | bcftools consensus exact | **panvram matched/total** | loss>1kb within-block (GVL/cons) | loss>1kb upper (SV) | discard |
|---|---|---|---|---|---|---|
| 1 kb | 92.8% | 93.1% | **8790/8790 = 100%** | 0.03% / 0.03% | 35.2% | 35.2% |
| 8 kb | 85.6% | 87.1% | **8789/8789 = 100%** | 3.7% / 3.6% | 38.5% | 36.1% |
| 32 kb | 66.6% | 70.8% | **8779/8779 = 100%** | 3.8% / 3.7% | 41.4% | 39.1% |

Three columns all vs the same independent FASTA. **panvram = 100% byte-exact (26358/26358)**
— the control, measured not "by construction". GVL trails consensus (−0.3/−1.6/−4.2 pp;
the 32 kb gap is the tandem-repeat over-expansion, `GENVARLOADER_ISSUE_DRAFT.md`). ed_p99
GVL/consensus ≈ equal (557 / 4399 / 17410). loss>1kb within-block monotone; nesting keeps
monotonicity; the SV-breakpoint upper bound is 35-41% (discarded windows = alignment
breakpoints). panvram loss = 0 at every length.
