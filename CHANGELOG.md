# Changelog

## 1.0.0 - 2026-10-09

- Version 1.0.0 of the API of 0.1.0 (`Cohort.open`, `sample`, `fetch`, `windows`, `decode`, `check_full`, `fasta`),
  format refrel3 v1 unchanged (FORMAT.md).
- README: measured numbers only with their evidence; Limitations (size against AGC, resident T2T, CUDA sm_80+).
- DATASET.md and SHA256SUMS of the HPRC v1 archives (q4k, q16k).
- CI workflow (CPU path, synthetic fixtures).
- Colab notebook: the 558-assembly gate on the dataset.
- CITATION.cff.
- Opt-in compact resident forms, format unchanged, default path unchanged: `Cohort.open(..., packed_reference=True)`
  (reference in 2 bits per base + table of non-ACGT runs) and `compact_blocks=True` (two-level block table: per 32
  blocks offset + state, per block 16-bit length + 16-bit state code, exception table). CPU path tested bit-exact
  (synthetic, and the 8-assembly gate in all four forms: `logs/gate-cpu-compact-2026-10-05.log`); bytes on the CPU path
  for the 558 cohort: `logs/resident-558-2026-10-06.log`. On the GPU: 558 x 4 forms PASS on A100-SXM4-80GB
  (`logs/colab-c1-558-forms-A100-2026-10-07.log`).
- Reader: two checks of the byte-level specification (a symbol from an all-zero context fails the block; the meta is
  exactly one zstd frame); check_full 1116 / 1116 (`logs/check-full-1116-2026-10-07.log`).
- EVIDENCE.md: five statements with numbers, commits, logs, SHA-256, reproduction and limitations.
- CPU container: Dockerfile with the base image pinned by digest, requirements.lock with wheel SHA-256.
- LICENSE: full MIT text, copyright Yakiv Shavidze.
- Removed before release: `notebooks/gate558_colab.ipynb` (per-archive download with record placeholders; the public
  quickstart covers the gate) and the local RELEASE_CHECKLIST.md.
- Dataset: the 558 refrel3 v1 archives (q4k, q16k) as one tar per dataset shipped in 500 MiB parts on Zenodo,
  DOI 10.5281/zenodo.23232317 (`zenodo/FILES.md`: 48 files, every part with its SHA-256 in `PARTS.sha256`, the whole tar
  with `cohort558_<dataset>.tar.sha256`; DATASET.md: how to rebuild the tar from its parts).
- Window law (EVIDENCE statement 2): the hw-apex verdict job on HPRC N = 4 - refrel3 q4k / q16k PASS (|error| <= 3.3 % at
  1 / 8 / 64 KiB), model A for zstd / LZ4 / OZSEG (`logs/window-law-B-2026-10-08/`).
- M6 (EVIDENCE statement 4): AGC windows 868 -> 26 /s as N grows 50 -> 558 (one thread), refrel3 q4k 200 884 -> 145 222 /s,
  counters and curves (`logs/m6-2026-10-07/`).
- `logs/`: every log cited by EVIDENCE.md is a file of this repository (`logs/SHA256SUMS`).
- Public notebook: the dataset tar from its parts (`PARTS.sha256` -> cat -> tar SHA-256 -> MANIFEST.tsv), MANIFEST.tsv of
  the package checked by SHA-256, no Google Drive.

## 0.1.0 - 2026-10-04 (local)

- refrel3 v1 reader with every check, CPU decoder, CUDA queue kernel ported from aceapex research/refrel @ 5b6d5ce.
- `Cohort.open` / `sample` / `fetch` / `fasta`; reverse complement, tokens ACGTN -> 0..4.
- Public gate test, synthetic tests, Colab scripts for the gate and the 558 cohort.
