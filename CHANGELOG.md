# Changelog

## 1.0.0 - unreleased (prepared locally)

- Version 1.0.0 of the API of 0.1.0 (`Cohort.open`, `sample`, `fetch`, `windows`, `decode`, `check_full`, `fasta`),
  format refrel3 v1 unchanged (FORMAT.md).
- README: measured numbers only with their evidence; Limitations (size against AGC, resident T2T, CUDA sm_80+).
- DATASET.md and SHA256SUMS of the HPRC v1 archives (q4k, q16k).
- CI workflow (CPU path, synthetic fixtures).
- Colab notebook: the 558-assembly gate on the dataset.
- CITATION.cff.

## 0.1.0 - 2026-10-04 (local)

- refrel3 v1 reader with every check, CPU decoder, CUDA queue kernel ported from aceapex research/refrel @ 5b6d5ce.
- `Cohort.open` / `sample` / `fetch` / `fasta`; reverse complement, tokens ACGTN -> 0..4.
- Public gate test, synthetic tests, Colab scripts for the gate and the 558 cohort.
