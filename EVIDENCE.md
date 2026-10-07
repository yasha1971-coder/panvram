# EVIDENCE - the five statements about panvram / refrel3 v1, and the container (R2)

Rule: every number below is copied from a log named in the same row; a statement whose log is not in a repository
on ace-core is marked **NO LOG HERE** and is not supported by this file until the log is added. SHA-256 of each artifact
as of 2026-10-07; aceapex paths are on branch `refrel` of aceapex (local commits), hw-apex paths on PR #63.

## 1. 558 assemblies resident on one device

| | |
|---|---|
| number (log) | q4k, 558 assemblies + T2T resident, default form: 16 642 219 335 B (reference 3 117 292 070, payload 8 448 112 345, block table 4 923 167 388, model tables 153 647 532), open 12.9 s; q16k 11 563 752 056 B, open 10.6 s |
| where it was measured | **CPU path on ace-core** (`resident_bytes()` of the pools that the CUDA path copies to the device) - not a GPU run |
| full decode | check_full of all 1 116 archives (558 x q4k / q16k): 1116/1116 XXH3 == header == MANIFEST, CPU path (reader 087c638) |
| commit / log | panvram cf4efd0 `logs/resident-558-2026-10-06.log` (sha256 6e46b65e...); panvram 087c638 `logs/check-full-1116-2026-10-07.log` (e0f0f539...) |
| artifacts | MANIFEST.tsv (3b6941ef...) = SHA-256 of the 1 116 archives (SHA256SUMS.q4k / .q16k) |
| reproduce | `python3 scripts/resident.py <v1 dir> q4k` (resident bytes); `python3 scripts/check_full_all.py <v1 dir> <manifest_v1.tsv> <t2t.fa>` |
| GPU part of the statement | 16.642 GB **VRAM**, 558/558 XXH3 **on the GPU**, 1.40 M windows/s, GPU type: **NO LOG HERE** (Colab run of `run_colab_cohort558.sh`; its log is not on ace-core) |
| limitation | until the Colab log is added, the 558-on-one-GPU statement rests on CPU-path byte counts only. |

## 2. Window law

| | |
|---|---|
| statement | t(W) = c0 + (W + Q - 1) / D_Q within ±20 %, model B, W = 1 outside the verdict |
| D_Q inputs (log) | one-thread full decode, HPRC N = 4, ace-core: refrel3 q4k 0.876 GB/s, q16k 0.915 GB/s (3 warm-ups + 9 runs, 48/48 == source) |
| commit / log | hw-apex PR #63 `review/axis3/results/ace-core-2026-10-05-407088c-dq/SUMMARY.md` (sha256 b34333a8...); protocol v1.1 = `review/axis3/PROTOCOL_AXIS3.md` at 6f573bb (f0bb442a...) |
| verdict | **NO LOG HERE** for a model-B verdict (c0, ±20 %) of refrel3 on ace-core. The only window-law check on record is the user's Colab Blackwell run (aceapex `research/refrel/RESULTS.md` section 3c: error -16 ... +7 % at 262 144 windows per batch; saturation model; log not in the repository) |
| other formats | 4 or more other formats: after S2 native integration (not done) |
| limitation | our format: D_Q measured; the verdict row needs the window run under protocol v1.1. |

## 3. Archive bytes independent of threads and CPU, given the build recipe

| | |
|---|---|
| number (log) | encoder `refrel3v1` @ 5b6d5ce, g++ 11.4.0 `-O3 -march=x86-64-v3 -funroll-loops`: y1_HG00438.1 q4k and y1_HG00621.1 q16k == manifest SHA-256; -march native / x86-64-v4 / haswell / skylake-avx512 / icelake-server / sapphirerapids / znver3 same bytes; threads 1 / 4 / 8 / 12 / 16 same bytes; without FMA (x86-64-v2 or -ffp-contract=off) y1_HG00438.1 q4k 21 765 975 B (8c4a13f8...) instead of 21 765 741 B (7a0b3cc9...), both decode == FASTA |
| commit / log | aceapex 7e4387a `research/refrel/logs/cohort-v1-BUILD.md` (sha256 3e068c12...) |
| clean-room rebuild | aceapex 05477c2 `research/cleanroom/repro.sh`: frozen tool rebuilt from 5b6d5ce, binary b7d1842d..., 24 fixtures byte-identical, package zip 75ff3f60... reproduced |
| reproduce | the g++ line of cohort-v1-BUILD.md; `bash research/cleanroom/repro.sh <dir>` |
| limitation | another compiler can give other bytes. gcc 13.3.0 (Colab) gave 21 757 194 B for y1_HG00438.1 q4k - **the cause is not established** (ledger C-101: hypothesis; the log above shows FMA changes bytes on gcc 11.4, not that FMA is the gcc 13.3 difference); Colab log NO LOG HERE. |

## 4. One archive, several readers

| reader | number (log) | commit / log |
|---|---|---|
| CPU (panvram) | 1116/1116 check_full; gate 8 assemblies q4k + q16k: 1024/1024 windows, 1000/1000 fetches == FASTA | panvram 087c638 `logs/check-full-1116-2026-10-07.log`; f5feecd `logs/gate-cpu-2026-10-04.log` (cb0e2c08...) |
| independent implementation from the specification (Python), clean-room round 2 | a new agent with no sources, only the package `refrel3_cleanroom_v2.zip` (SHA-256 291ae64d1ca43745f4063b49d5c2b310dc34bae3e701ec2147b4737ffc87bfb6): full decodes 26/26, fetches 1 300/1 300, refusals 9/9 at the expected stage, test vectors 372/372; gaps reported 9, blocking 0 | round-2 report of the user, 2026-10-07 - **NO LOG HERE** until the report file is added; package on Drive gdrive:ACEAPEX-OVH/outgoing (check 0 differences) |
| our check decoder from the specification | Python (not C99): 26/26 full decodes, 1 300/1 300 fetches, 9/9 refusals on the clean-room package; ops == frozen code on every block (92 938 + 1 208); rANS trace == frozen code (9 358 events) | aceapex 83ebfd9 `research/cleanroom/v2/REPORT_V2.md` (ebb65bd6...), `specdec.py` |
| GPU | Colab gate, 8 archives (sha256 8/8 == MANIFEST), q4k and q16k PASS | **NO LOG HERE** (MyDrive/aceapex_logs/panvram_<date>.txt); reader 087c638 not yet run on a GPU |
| region, haplotype coordinates | samtools 1.24 faidx on the 8 source FASTA vs panvram fetch, protocol v1.1 coordinates [start0,end0): 8000/8000 q4k, 8000/8000 q16k | panvram 0a5d8ff `logs/samtools-truth-2026-10-07/RESULT.txt` (919de928...), requests.tsv (359caa23...) |
| reproduce | `scripts/samtools_truth.py <work> .data/cohort sources.tsv 1000`; `scripts/check_full_all.py ...` | |
| limitation | reference (T2T) coordinates: not done (L1; aceapex research/bench2610/f_liftover). A C99 decoder of refrel3 v1 does not exist; the independent implementation is Python. |

## 5. Corruption

| | |
|---|---|
| number (log) | 6 kinds x 10 000 cases per configuration, one reference archive: q4k with block XXH3: SILENT 0 of 60 000; q16k with block XXH3: SILENT 0 of 60 000; q4k without block hashes: SILENT 1 798 of 60 000; q16k without: SILENT 2 145 of 60 000; hang 0, crash 0 in all 240 000 |
| commit / log | aceapex 7e4387a `research/refrel/logs/v1-corrupt-2026-10-04.log` (sha256 740295ef...) |
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
