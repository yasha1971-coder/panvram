# EVIDENCE - the five statements about panvram / refrel3 v1, and the container (R2)

Rule: every number below is copied from a file under `logs/` of this repository (`logs/SHA256SUMS`) named in the same row;
a statement without such a file is marked **NO LOG HERE** and is not supported by this file until the log is added. SHA-256 of each artifact
as of 2026-10-08. Every log cited below is a file of this repository under `logs/` (`logs/SHA256SUMS`); the encoder and the
research directories are referenced by tag of the aceapex repository (links `https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/` + path; tag `panvram-1.0-sources` of aceapex = one orphan commit whose tree holds exactly three subtrees of the aceapex history: `research/refrel` = `5b6d5ce:research/refrel` (tree c96b8f0167c54ae6645aa18d7c30c697a18178a0), `research/cleanroom` = `927a9dd:research/cleanroom` (tree 1208ed2af9b373779cefc05edb94cd5f7c5554d7), `research/bench2610` = `927a9dd:research/bench2610` (tree 3649681f1517406bf752e8bbdc8e727887a2143f); check in an aceapex clone: `git rev-parse panvram-1.0-sources:research/refrel panvram-1.0-sources:research/cleanroom panvram-1.0-sources:research/bench2610` == `git rev-parse 5b6d5ce:research/refrel 927a9dd:research/cleanroom 927a9dd:research/bench2610`), hw-apex paths by PR #63 of github.com/yasha1971-coder/hw-apex-bench.

## 1. 558 assemblies resident on one device

| | |
|---|---|
| number (log) | q4k, 558 assemblies + T2T resident, default form: 16 642 219 335 B (reference 3 117 292 070, payload 8 448 112 345, block table 4 923 167 388, model tables 153 647 532), open 12.9 s; q16k 11 563 752 056 B, open 10.6 s |
| where it was measured | **CPU path on ace-core** (`resident_bytes()` of the pools that the CUDA path copies to the device) - not a GPU run |
| full decode | check_full of all 1 116 archives (558 x q4k / q16k): 1116/1116 XXH3 == header == MANIFEST, CPU path (reader 087c638) |
| commit / log | panvram cf4efd0 `logs/resident-558-2026-10-06.log` (sha256 6e46b65e...); panvram 087c638 `logs/check-full-1116-2026-10-07.log` (e0f0f539...) |
| artifacts | MANIFEST.tsv (3b6941ef...) = SHA-256 of the 1 116 archives (SHA256SUMS.q4k / .q16k) |
| reproduce | `python3 scripts/resident.py <v1 dir> q4k` (resident bytes); `python3 scripts/check_full_all.py <v1 dir> <manifest_v1.tsv> <t2t.fa>` |
| GPU part of the statement (log) | Colab **NVIDIA RTX PRO 6000 Blackwell Server Edition** (97 887 MiB, driver 580.82.07, sm_120, CUDA 13.0, torch 2.11.0+cu130), panvram 02b23d0, q4k, 558 assemblies (0 missing): open 12.9 s, resident 16.642 GB (reference 3.117, payload 8.448, block table 4.923, model tables 0.154), torch allocated 16.646 GB, nvidia-smi used 0 -> 16 441 MiB; sample 1024 x 8192: median 0.73 ms per call = 1 401 185 windows/s (11.48 GB/s, 5 calls, 474 assemblies hit) == CPU decoder; fetch 1000/1000 == CPU decoder (median 0.10 ms); full decode on the GPU 558/558 XXH3 == manifest (1.680 T bases, 3.5 min); RESULT PASS |
| GPU commit / log | panvram `logs/colab-cohort558-blackwell-2026-10-05.txt` (sha256 2097320b...) |
| GPU, four resident forms (C1, log) | Colab **NVIDIA A100-SXM4-80GB** (81 920 MiB, driver 580.82.07, sm_80, CUDA 13.0), panvram 56c32dc (reader 087c638), q4k, 558 assemblies, every form PASS (sample == CPU decoder, 1000/1000 fetches, full decode 558/558 XXH3 == manifest): default 16.642 GB (nvidia-smi 0 -> 16 310 MiB, 743 271 windows/s); 2-bit reference 14.304 GB (14 080 MiB, 674 822 /s); compact block table 13.972 GB (13 764 MiB, 710 235 /s); **both compact forms 11.634 GB resident (nvidia-smi 0 -> 11 534 MiB), 736 762 windows/s (sample 1024 x 8192, median of 5 calls), 558/558 ==**. Cohort-region stage S3 not run (no request files on Drive) |
| C1 log | panvram `logs/colab-c1-558-forms-A100-2026-10-07.log` (sha256 7766b708...) |
| limitation | q4k only; windows/s at one batch size (1 024 x 8 192, median of 5 calls); A100 and Blackwell figures are different runs, not compared. The compact runs used an 80 GB card: "fits in 16 GB" (11.634 GB / 11 534 MiB) is from the measured memory, not from a run on a 16 GB card. |

## 2. Window law

| | |
|---|---|
| statement | t(W) = c0 + (W + Q - 1) / D_Q within ±20 %, model B, W = 1 outside the verdict |
| D_Q inputs (log) | one-thread full decode, HPRC N = 4, ace-core: refrel3 q4k 0.876 GB/s, q16k 0.915 GB/s (3 warm-ups + 9 runs, 48/48 == source) |
| commit / log | `logs/hwapex-dq-2026-10-05-SUMMARY.md` (copy of hw-apex PR #63 `review/axis3/results/ace-core-2026-10-05-407088c-dq/SUMMARY.md`); protocol v1.1 = hw-apex `review/axis3/PROTOCOL_AXIS3.md` at 6f573bb (f0bb442a...) |
| D_Q, two numbers | 0.876 GB/s (axis 3, 2026-10-05): one-thread decode **to the FASTA file** (bases + headers + line breaks rebuilt); 2.068 GB/s (verdict job, 2026-10-08): the same decoder's **canonical sequence bytes** only, no FASTA formatting, no block-hash FASTA XXH3 inside the timer - the two timed boundaries differ, neither number is wrong |
| verdict (log) | hw-apex B job (`tools.verdict_refrel3`, protocol v1.1 4a0027fe..., harness 0ace54d) on ace-core 2026-10-08, HPRC N = 4, one thread, 10 000 requests per W, D_Q 3 + 9: **the window law is confirmed on refrel3 and BGZF** - refrel3 q4k PASS (Q 4 096, D_Q 2.068 GB/s canonical, c0 7.25 us; errors at 1 / 8 / 64 KiB +1.5 / +1.3 / +3.2 %), q16k PASS (Q 16 384, D_Q 2.085 GB/s, c0 6.79 us; -0.1 / -2.7 / +3.2 %), BGZF default PASS (-1.2 / -10.0 / +1.8 %), BGZF matched-g PASS (-9.5 / -0.1 / +3.2 %). **For the large-granule formats model B is conservative (c0 < 0 -> FAIL by the rule): see model A** - zstd seekable 16 KiB frames, LZ4 4 MiB blocks, OZSEG 64 KiB / 1 MiB (B predictions within 2.4 % for LZ4 and 12 % for OZSEG, -43 / -41 % for zstd at 8 / 64 KiB; model A overestimates by 14-66 %); table: `logs/window-law-B-2026-10-08/MODEL_A_LARGE_GRANULE.md`. AGC FAILED (no canonical-granule evidence). verify: VERDICT_EVIDENCE_PASS |
| verdict log | `logs/window-law-B-2026-10-08/` (RESULTS_WINDOW_LAW_B.md, table.md, results.json, MODEL_A_LARGE_GRANULE.md; the raw window logs, 1.1 GB, are retained off-repository) |
| limitation | N = 4 cohort, one host, one thread. The law is confirmed on refrel3 and BGZF; for formats whose granule is >= the window (LZ4 4 MiB, OZSEG 1 MiB) or whose sequential D_Q includes output materialization, c0 comes out negative and model B is conservative - model A is then the upper bound. AGC not measured (no granule evidence). |

## 3. Archive bytes independent of threads and CPU, given the build recipe

| | |
|---|---|
| number (log) | encoder `refrel3v1` @ 5b6d5ce, g++ 11.4.0 `-O3 -march=x86-64-v3 -funroll-loops`: y1_HG00438.1 q4k and y1_HG00621.1 q16k == manifest SHA-256; -march native / x86-64-v4 / haswell / skylake-avx512 / icelake-server / sapphirerapids / znver3 same bytes; threads 1 / 4 / 8 / 12 / 16 same bytes; without FMA (x86-64-v2 or -ffp-contract=off) y1_HG00438.1 q4k 21 765 975 B (8c4a13f8...) instead of 21 765 741 B (7a0b3cc9...), both decode == FASTA |
| commit / log | `logs/aceapex-cohort-v1-BUILD-2026-10-04.md` (copy of aceapex refrel@927a9dd `research/refrel/logs/cohort-v1-BUILD.md`, written after 5b6d5ce and therefore not in the tag; the copy under `logs/` is the record); encoder source: https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/research/refrel/refrel3v1.cpp |
| clean-room rebuild | aceapex `research/cleanroom/repro.sh` (https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/research/cleanroom/repro.sh): frozen tool rebuilt from 5b6d5ce, binary b7d1842d..., 24 fixtures byte-identical, package zip 75ff3f60... reproduced |
| reproduce | the g++ line of cohort-v1-BUILD.md; `bash research/cleanroom/repro.sh <dir>` |
| limitation | another compiler can give other bytes. gcc 13.3.0 (Colab) gave 21 757 194 B for y1_HG00438.1 q4k (all 8 VM-encoded archives differ from MANIFEST, each decodes == FASTA; panvram `logs/colab-gate-blackwell-2026-10-04.txt`) - **the cause is not established** (ledger C-101: hypothesis; FMA is shown to change bytes on gcc 11.4, not shown to be the gcc 13.3 difference). |

## 4. One archive, several readers

| reader | number (log) | commit / log |
|---|---|---|
| CPU (panvram) | 1116/1116 check_full; gate 8 assemblies q4k + q16k: 1024/1024 windows, 1000/1000 fetches == FASTA | panvram 087c638 `logs/check-full-1116-2026-10-07.log`; f5feecd `logs/gate-cpu-2026-10-04.log` (cb0e2c08...) |
| independent implementation from the specification (Python), clean-room round 2 | a new agent with no sources, only the package `refrel3_cleanroom_v2.zip` (SHA-256 291ae64d...): full decodes 26/26, fetches 1 300/1 300, refusals 9/9 at the expected stage, test vectors 372/372; gaps reported 9, blocking 0 | `logs/cleanroom-round2-2026-10-07/` (RESULTS_R2.md, results.json, SPEC_GAPS_R2.md 22d493f6...) from cleanroom_round2.tar.gz (72187a3b...); refusals: stage 9/9, exact reason string 4/9 (§15 defines no decode-stage strings: gap G2) |
| our check decoder from the specification | Python (not C99): 26/26 full decodes, 1 300/1 300 fetches, 9/9 refusals on the clean-room package; ops == frozen code on every block (92 938 + 1 208); rANS trace == frozen code (9 358 events) | aceapex `research/cleanroom/v2/REPORT_V2.md` (https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/research/cleanroom/v2/REPORT_V2.md), `research/cleanroom/v2/specdec.py` (https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/research/cleanroom/v2/specdec.py) |
| GPU, reader 087c638 (with the two stricter checks) | Colab **NVIDIA A100-SXM4-40GB** (sm_80, driver 580.82.07, CUDA 13.0, torch 2.11.0+cu130), panvram e768115, 4 assemblies: q4k PASS and q16k PASS - 1024/1024 windows of 8 192 == FASTA and == the CPU decoder, 1000/1000 fetches == FASTA, 4 assemblies rebuilt, XXH3 == source and == source files; fast tests 30 passed. A100 numbers are not compared with the Blackwell run of 05.10 | panvram `logs/colab-gate-087c638-A100-2026-10-07.txt` (sha256 d0695bf8...) |
| GPU, earlier gate | Colab RTX PRO 6000 Blackwell (sm_120), panvram ed986f0: archives from ace-core 8/8 sha256 == MANIFEST; q4k and q16k gate PASS (4 assemblies, 1024/1024 windows, 1000/1000 fetches, FASTA XXH3 == source) | panvram `logs/colab-gate-blackwell-2026-10-04.txt` (sha256 e0e647e6...) |
| GPU, quickstart rehearsal from the Zenodo draft (the bytes of record 23232317), Colab A100, panvram 989d415, unattended run | 24/24 parts == PARTS.sha256, cohort558_q4k.tar == .sha256, 558/558 archives == MANIFEST.tsv; compact form (packed reference + compact blocks) 11.634 GB, nvidia-smi 11 534 MiB, open 55.1 s, 1024 x 8192 windows 736 335 /s == CPU decoder, 1000/1000 fetches == CPU decoder, C558 RESULT PASS | panvram 989d415 `logs/quickstart-dryrun-A100-2026-10-08.txt` (sha256 49df845e...) |
| region, haplotype coordinates | samtools 1.24 faidx on the 8 source FASTA vs panvram fetch, protocol v1.1 coordinates [start0,end0): 8000/8000 q4k, 8000/8000 q16k | panvram 0a5d8ff `logs/samtools-truth-2026-10-07/RESULT.txt` (919de928...), requests.tsv (359caa23...) |
| reproduce | `scripts/samtools_truth.py <work> .data/cohort sources.tsv 1000`; `scripts/check_full_all.py ...` | |
| limitation | reference (T2T) coordinates: not done (aceapex `research/bench2610/f_liftover`, https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/research/bench2610/f_liftover). A C99 decoder of refrel3 v1 does not exist; the independent implementation is Python. |

## 5. Corruption

| | |
|---|---|
| number (log) | 6 kinds x 10 000 cases per configuration, one reference archive: q4k with block XXH3: SILENT 0 of 60 000; q16k with block XXH3: SILENT 0 of 60 000; q4k without block hashes: SILENT 1 798 of 60 000; q16k without: SILENT 2 145 of 60 000; hang 0, crash 0 in all 240 000 |
| commit / log | `logs/aceapex-v1-corrupt-2026-10-04.log` (copy of aceapex refrel@927a9dd `research/refrel/logs/v1-corrupt-2026-10-04.log`, written after 5b6d5ce and not in the tag; sha256 740295ef...) |
| reproduce | `refrel3v1 corrupt <ref.fa> <archive> 10000 <verify 0/1> <seed> <parallel>` |
| limitation | "1 798 of 2 145" is not a row of the log: 1 798 (q4k) and 2 145 (q16k) are two configurations, each of 60 000. Without block hashes a fetch is protected only by the rANS checks (FORMAT_V1_SPEC section 14). |

## R2. CPU container

| | |
|---|---|
| number (log) | 8 assemblies x q4k / q16k inside the container, `--network none`: check_full XXH3 == header == MANIFEST and full decode SHA-256 == source, 16/16, 151 s |
| image | base `ubuntu:22.04@sha256:5ec03bb3441e8b0bf3b4f9cd4629a1ae763010dc3035bb8da3ae6cf026486401`; Python packages `pip --require-hashes` (requirements.lock); image id sha256:cf4044401953aef19f7eb901c04166ae1a4063904b8fa75930b39af19f5e10d1 (built under rootless dockerd 29.1.2) |
| commit / log | panvram af2ec56 `logs/r2-gate-2026-10-07.log` (sha256 f2a43dd7...), Dockerfile, requirements.lock |
| reproduce | `docker build -t panvram-cpu .` then the `docker run --network none ...` line in the log |
| limitation | apt packages inside the image are not pinned by version (only the base image and the Python wheels); the gate ran with reader af2ec56, before the two checks of 087c638. |

## M6. AGC against the cohort size N (request order matters)

| | |
|---|---|
| number (log) | AGC 3.2.4 random windows of 4 096 bases, one thread: 868 / 38 / 28 / 26 requests per s at N = 50 / 100 / 200 / 558 (16 threads 3 331 / 139 / 109 / 97); refrel3 q4k 200 884 / 178 591 / 186 487 / 145 222 (one thread); 1 Mb regions AGC 294 -> 26 /s; whole samples flat (AGC 7.6 s, refrel3 1.3 s per sample). Every answer == truth |
| caveat | **the AGC numbers hold for requests in random sample order.** Counters on the same 1 000 windows: in random order every window at N >= 100 reloads one metadata batch (3.8 MB zstd in -> 15.6 MB out, 105 MB of allocations); **sorted by sample** the same windows cost the same at every N (0 batch loads, 17 KB -> 15 KB). A workload that groups requests by sample does not see the N dependence; a workload that does not, does. refrel3 has no such order effect (counters not taken on refrel3). |
| commit / log | `logs/m6-2026-10-07/` (RESULTS_M6.md, curves_counters.csv), `logs/m6-curves-timing-2026-10-07.csv`; scripts and raw logs in aceapex `research/bench2610/m6` (https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/research/bench2610/m6) |
| reproduce | aceapex `research/bench2610/m6/run_all.sh` (https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/research/bench2610/m6/run_all.sh; PROTOCOL_M6 v1.1: https://github.com/yasha1971-coder/aceapex/tree/panvram-1.0-sources/research/bench2610/m6/PROTOCOL_M6.md) |
| limitation | one AGC build and create setting (-b 50 batch, -s 60000); the N effect is the metadata batch of the archive layout, measured, not a statement about AGC with other settings. |
