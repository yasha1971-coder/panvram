# panvram

A pangenome cohort resident on the GPU: every assembly stored as refrel3 v1 (edits against one decoded reference,
T2T-CHM13v2.0 for HPRC; [FORMAT.md](FORMAT.md)), random training windows and coordinate slices decoded on the card by
a queue kernel, straight into a PyTorch tensor. A CPU decoder serves the same API without a GPU.

Status: local, not published. Numbers below come only from the logs named next to them.

```python
import torch, panvram
c = panvram.Cohort.open("cohort/", device="cuda")      # reference + every *.q4k.rr3, resident on the card
g = torch.Generator(device="cuda").manual_seed(0)
x = c.sample(1024, 8192, generator=g)                  # uint8 [1024, 8192] on the GPU, the FASTA's bytes
t = c.sample(1024, 8192, generator=g, tokens=True,     # A/a 0, C/c 1, G/g 2, T/t 3, anything else 4
             reverse_complement=0.5)                   # each window reverse-complemented with probability 0.5
y = c.fetch("y1_HG00438.1", "HG00438#1#JAHBCB010000001.1", 1_000_000, 5000)   # 0-based start, length
```

## Install

```
pip install torch            # first; the extension is built against it
pip install -e .             # CUDA when nvcc and a CUDA build of torch are present (sm_80 or newer), else CPU only
python -c "import panvram; print(panvram.with_cuda)"
```
Needs libzstd (`libzstd-dev`). `PANVRAM_CUDA=0|1` forces the path; `TORCH_CUDA_ARCH_LIST` overrides the GPU
architectures (default: the visible GPUs of sm_80 or newer, else `8.0;8.6;8.9;9.0+PTX`).

## API

`Cohort.open(path, device="cuda", dataset="q4k", reference=None, threads=0, assemblies=None)` - every
`*.<dataset>.rr3` in `path` (sorted by name; the assembly name is the file name without `.<dataset>.rr3`) and the
reference: `reference`, else `path/reference.fa`, else the file named in the archives' header in `path`. Every archive
is checked before any decode (FORMAT.md section 5: header XXH3, reference SHA-256 and size, section sizes, tables); a
wrong reference or a damaged header raises `panvram.FormatError`. One block size per cohort. The reference and the
pooled archives (payloads, block offsets and start states, model tables, lower-case runs) go to `device` once.

- `sample(n, W, generator=None, reverse_complement=False, tokens=False, return_coords=False, check=True)` - n windows
  uniform over every (assembly, contig, start) with start + W <= contig length; a window never crosses contigs.
  `return_coords=True` adds int64 [n, 4]: assembly, contig index, start in the contig, reverse-complemented.
- `fetch(assembly, contig, start, length, reverse_complement=False, tokens=False)` - a slice by coordinates; assembly
  and contig by name (header up to the first space or tab) or index.
- `windows(asm, start, W, reverse_complement=None, tokens=False, check=True)` - windows at stream offsets.
- `fasta(assembly, verify=True)` - full decode on the CPU -> the FASTA file's bytes; every block XXH3 and the FASTA
  XXH3 against the header (the source file's).
- `contigs(assembly)`, `info(assembly)`, `resident_bytes()`, `block_size`, `reference_sha256`, `names`.

Output bytes are the FASTA's (case kept: lower-case runs applied after decode); reverse complement keeps case and
leaves anything but ACGT/acgt as is. `check=True` waits for the kernel and raises if a block failed to decode (status
bits: 1 decode error, 4 bounded wait expired, 8 window outside an assembly). The GPU path does not check block XXH3
(the CPU `fasta` and `windows(..., verify)` do).

## How the GPU path works

The queue kernel (from aceapex `research/refrel/refrel3_gpu.cu`): one warp per refrel3 block a window touches; lane
0 decodes the block's rANS stream (`r3_decode_stream`, the code the CPU path runs), literals straight into the block
buffer in shared memory and copies into a ring of 64 ops; lanes 1..31 execute the copies as they arrive (reference,
reverse complement, in-block); the warp writes the window's part. Shared memory per warp: Q + 768 B. Then two small
kernels: lower-case runs, and reverse complement / tokens over pairs (i, W-1-i). Sampling (randint over the valid
starts, searchsorted over the contig table) runs on the device; no host round trip per batch except `check`.

## Tests

- `tests/test_synth.py` (seconds): synthetic reference and two assemblies (`tests/synth.py`, deterministic; contigs
  forward, reverse complement, inversion, novel + tandem repeat, N runs, IUPAC, lower case, lengths around block
  sizes), archives in `tests/data` encoded by aceapex `refrel3v1` @ 5b6d5ce for Q = 1024 / 4096 / 16384 (each decoded
  back by it == the FASTA). Every path (CPU, CUDA when present) against windows cut from the FASTA bytes
  (`tests/fasta_cut.py`), RC and tokens, whole contigs, CUDA == CPU, refusals (wrong reference, header byte, header
  XXH3, truncation, version, mixed block sizes), payload bit flips caught or harmless under full verification.
- `tests/test_gate.py` - **the public gate**: `import panvram` -> `Cohort.open` -> `sample(1024, 8192)` on the GPU ->
  the same windows cut from the FASTA the CPU decoder rebuilds (every block XXH3 and the FASTA XXH3 == the source's;
  `PANVRAM_SOURCES` also compares with the source files) byte for byte, and == the CPU decoder's own windows; fetch at
  1000 random coordinates (contig by length, length 1 .. 2^20) == the FASTA. `PANVRAM_COHORT=dir pytest -s tests/test_gate.py`.
- `run_colab_panvram.sh` - CUDA build, both tests and the gate for q4k and q16k on four HPRC haplotypes, throughput
  (`scripts/bench.py`). From Drive: the tarball (`MyDrive/panvram/panvram.tar.gz`, `git archive`) and the eight
  archives of MANIFEST.tsv (`MyDrive/panvram/cohort/*.rr3` + `SHA256SUMS`, checked); without them it encodes on the VM.

## Measured

| what | number | log |
|---|---|---|
| fast tests, CPU (ace-core) | 5 passed, 1 skipped (no CUDA) | `logs/tests-cpu-2026-10-04.log` |
| gate on the CPU path, 8 HPRC assemblies + T2T, q4k and q16k | 1024 / 1024 windows == FASTA, 1000 / 1000 fetches == FASTA, 8 FASTA rebuilt == source (XXH3), each dataset | `logs/gate-cpu-2026-10-04.log` |
| resident, 8 assemblies, q4k | reference 3.117 GB + payload 0.115 + block table 0.070 + model tables 0.002 GB | same |
| resident, 8 assemblies, q16k | reference 3.117 GB + payload 0.095 + block table 0.018 + model tables 0.002 GB | same |
| v1 archive per assembly (with block XXH3), 558 of MANIFEST.tsv, each decoded back == source (XXH3) | q4k 22.10 MB, q16k 14.51 MB mean (total 12.33 / 8.10 GB); FASTA 1.707 TB | `MANIFEST.tsv` |
| cohort, 558 HPRC assemblies, refrel3 before v1 (Q 16384, no block hashes), all decoded == FASTA | 13.04 MB per assembly | aceapex `research/refrel/logs/cohort-table-2026-10-04.txt` @ 6de4666 |
| GPU path | not measured yet (`run_colab_panvram.sh`) | - |

Arithmetic on the rows above (not a measurement): resident per assembly (payload + block table + tables) q4k
23.4 MB, q16k 14.4 MB; 558 assemblies + T2T: q4k about 16.2 GB, q16k about 11.1 GB.

## MANIFEST.tsv

The 558 HPRC haplotype assemblies (year-1 94, release 2 464; CHM13 v1.1 and GRCh38 rows of the indexes excluded):
name, source URL, source hash type and hash (from the HPRC indexes, checked at download), FASTA bytes, FASTA XXH3,
q4k and q16k archive bytes and SHA-256, status. The v1 archives are re-encoded from the cohort's refrel3 archives
(FASTA in memory, never re-downloaded); each is decoded back in full and compared with the source's XXH3 before its row
is written (aceapex `refrel3v1 recode`); all 558 rows `ok`.

Archive bytes depend on the encoder build (FORMAT.md section 7): the manifest's sha256 are those of gcc 11.4.0
`-O3 -march=native -funroll-loops` on ace-core, reproduced by `-O3 -march=x86-64-v3 -funroll-loops` (FMA contraction
on; checked on two archives). Without FMA (`-march=x86-64-v2`, `-ffp-contract=off`) or with another compiler the
bytes can differ while decoding to the same FASTA - compare archives by decoding, or build the encoder with these flags.

## Provenance

Ported from aceapex (github.com/yasha1971-coder/aceapex, branch `refrel`) `research/refrel` @
5b6d5cec0f5962a561ac48822a1b5c48793a5b47:

| panvram | source |
|---|---|
| `csrc/refrel3.h`, `csrc/refrel_format.h`, `csrc/rr_sha256.h` | verbatim |
| `csrc/xxhash.h` | aceapex `src/xxhash.h` (xxHash 0.8, BSD-2-Clause) verbatim |
| `csrc/pv_v1.h` | `open_v1`, `block_v1`, `full_v1` of `refrel3v1.cpp`; `build_tab` of `refrel3.cpp`; `exec_fast` of `refrel.cpp`; FASTA parsing of `refrel_io.h` - exceptions instead of exit, arrays in pooled memory |
| `csrc/pv_cuda.cu` | `r3q_kernel` of `refrel3_gpu.cu` - runtime Q, many assemblies resident, window = (assembly, start); new: case / RC / token kernels |
| `FORMAT.md` | `research/refrel/FORMAT.md` (format unchanged) |

## Limits

Reference under 2^32 bases (32-bit copy positions in the kernel's ring); Q in {1024, 2048, 4096, 16384}, one per
cohort; an assembly's payload under 4 GiB; sm_80 or newer for the CUDA path. Encoding is not part of panvram (aceapex
`refrel3v1 encode`).

License: MIT (see LICENSE); `csrc/xxhash.h` BSD-2-Clause.
