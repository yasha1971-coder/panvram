# panvram

A pangenome cohort resident on the GPU: every assembly stored as refrel3 v1 (edits against one decoded reference,
T2T-CHM13v2.0 for HPRC; [FORMAT.md](FORMAT.md)), random training windows and coordinate slices decoded on the card by
a queue kernel, straight into a PyTorch tensor. A CPU decoder serves the same API without a GPU.

Version 1.0.0 (prepared locally, not published). Numbers below come only from the logs named next to them.

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
- `decode(assembly)` - the whole base stream decoded on the device in chunks -> host tensor; `check_full(assembly)` -
  that, the FASTA rebuilt from the contig table and hashed: XXH3 == the header's source XXH3 (no source file needed).
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
- `notebooks/quickstart_558.ipynb` - public quickstart: downloads the package, T2T (md5) and the archives listed in
  MANIFEST.tsv (SHA-256 checked), then the 558-assembly run of `scripts/cohort558.py` per resident form (dataset location:
  placeholder until published).
- `run_colab_panvram.sh` - CUDA build, both tests and the gate for q4k and q16k on four HPRC haplotypes, throughput
  (`scripts/bench.py`). From Drive: the tarball (`MyDrive/panvram/panvram.tar.gz`, `git archive`) and the eight
  archives of MANIFEST.tsv (`MyDrive/panvram/cohort/*.rr3` + `SHA256SUMS`, checked); without them it encodes on the VM.

## Measured (each number from a file in `logs/`)

| what | number | log |
|---|---|---|
| 558 assemblies resident on one GPU, q4k, four resident forms (Colab NVIDIA A100-SXM4-80GB, sm_80, panvram 56c32dc) | default 16.642 GB (nvidia-smi 16 310 MiB), 743 271 windows/s; 2-bit reference 14.304 GB, 674 822 /s; compact block table 13.972 GB, 710 235 /s; **both compact 11.634 GB (nvidia-smi 11 534 MiB), 736 762 windows/s**; every form: sample == CPU decoder, 1000 / 1000 fetches, full decode on the GPU 558 / 558 XXH3 == manifest, PASS | `logs/colab-c1-558-forms-A100-2026-10-07.log` |
| 558 assemblies resident, default form (Colab NVIDIA RTX PRO 6000 Blackwell Server Edition, sm_120, panvram 02b23d0) | 16.642 GB (nvidia-smi 0 -> 16 441 MiB), open 12.9 s, 1 401 185 windows/s (1 024 x 8 192, median of 5 calls), 1000 / 1000 fetches, 558 / 558 XXH3 == manifest, PASS | `logs/colab-cohort558-blackwell-2026-10-05.txt` |
| public gate on the GPU, reader with the two stricter checks (Colab NVIDIA A100-SXM4-40GB, panvram e768115), 4 HPRC assemblies, q4k and q16k | 1024 / 1024 windows of 8 192 == FASTA and == CPU decoder, 1000 / 1000 fetches, 4 FASTA rebuilt XXH3 == source files; fast tests 30 passed; PASS | `logs/colab-gate-087c638-A100-2026-10-07.txt` |
| check_full of every manifest archive, CPU path (ace-core) | 1116 / 1116 (558 x q4k / q16k) XXH3 == header == manifest | `logs/check-full-1116-2026-10-07.log` |
| region fetch against samtools 1.24 on the source FASTA, 8 assemblies x 1 000 requests, [start0, end0) | 8000 / 8000 q4k, 8000 / 8000 q16k | `logs/samtools-truth-2026-10-07/RESULT.txt` |
| CPU container, `--network none`, base image pinned by digest, wheels by hash | 8 assemblies x q4k / q16k: check_full and full decode == source SHA-256, 16 / 16 | `logs/r2-gate-2026-10-07.log` |
| independent implementation from the specification (Python, clean-room round 2: no sources, the specification package only) | full decodes 26 / 26, fetches 1 300 / 1 300, refusals 9 / 9 at the expected stage, test vectors 372 / 372 | `logs/cleanroom-round2-2026-10-07/RESULTS_R2.md` |
| archive sizes, 558 assemblies | q4k 16.22 MB, q16k 13.04 MB per assembly without block hashes; AGC 3.2.4 (one create, T2T inside, T2T alone subtracted) 5.40 MB | `logs/sizes-558-2026-10-06.txt` |

Statement by statement, with commits, SHA-256 of every artifact, reproduction commands and one limitation each:
[EVIDENCE.md](EVIDENCE.md).

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

## Limitations

- **Size against AGC.** On the 558 HPRC assemblies refrel3 takes 13.04 (q16k) - 16.22 (q4k) MB per assembly without
  block hashes against 5.40 MB for AGC 3.2.4 (`logs/sizes-558-2026-10-06.txt`): AGC also compresses against the other
  assemblies, refrel3 encodes every assembly alone against the reference.
- **"Fits in 16 GB" is a memory measurement, not a run on a 16 GB card.** The compact forms hold 558 assemblies in
  11.634 GB (nvidia-smi 11 534 MiB) on an 80 GB A100; no run was made on a 16 GB card.
- **Archive bytes are reproducible with the encoder recipe only.** gcc 11.4.0 `-O3 -march=x86-64-v3 -funroll-loops`
  (any thread count, any x86-64-v3+ CPU) gives the manifest's bytes; other compilers or flags can give other bytes
  (gcc 13.3 on Colab did) that decode to the same FASTA (FORMAT.md section 7).
- **Coordinates are the haplotype's own.** `fetch` takes an assembly's contig and 0-based offsets; there is no
  translation from reference (T2T) coordinates to haplotype coordinates.
- **The reference is resident decoded.** T2T-CHM13v2.0 sits on the device as 3.117 G bases (0.779 GB in the 2-bit
  form) for every cohort, on top of the archives.
- **GPU: CUDA only, sm_80 or newer.** Other devices use the CPU path. The GPU path does not check block XXH3 (the CPU
  `fasta` and `windows(..., verify)` do).
- One block size per cohort; reference under 2^32 bases; an assembly's payload under 4 GiB.
- Encoding is not part of panvram (aceapex `refrel3v1 encode`).

## Data use (HPRC)

The cohort is the Human Pangenome Reference Consortium's assemblies (public domain / CC0 on the AWS Registry of Open
Data; full terms quoted in [DATASET.md](DATASET.md)). The HPRC asks users not to sell all or part of the Pangenome
Reference, not to seek intellectual property protection that would block access to it, and not to take any action that
could reasonably be expected to result in participant re-identification, including contacting participants. When
publishing results that use it, acknowledge: "We would like to acknowledge the Human Pangenome Reference Consortium
(BioProject ID: PRJNA730823) and its funder, the National Human Genome Research Institute (NHGRI)."

License: MIT (see LICENSE), Yakiv Shavidze; `csrc/xxhash.h` is BSD-2-Clause (xxHash, its own notice kept).
Cite: [CITATION.cff](CITATION.cff) (ORCID 0009-0008-3622-3448).
