# Contributing

- Build: `pip install torch` first, then `pip install -e .` (CPU path without nvcc; CUDA path with nvcc and a CUDA
  torch, sm_80 or newer).
- Tests that must pass on the CPU path: `cd tests && python -m pytest -q test_synth.py` (synthetic fixtures; no data
  download). The public gate on a real cohort: `PANVRAM_COHORT=<dir> python -m pytest -s tests/test_gate.py`.
- The archive format is refrel3 v1 (FORMAT.md) and is frozen: changes that alter what a v1 reader accepts are a new
  format version, not a patch. Resident forms (`packed_reference`, `compact_blocks`) must stay bit-exact with the
  default path (the tests compare them).
- Numbers in README / BENCHMARKS come only from logs committed with the change (or named precisely); no estimates.
- Every source file starts with `SPDX-License-Identifier: MIT`; contributions are under the MIT license.
