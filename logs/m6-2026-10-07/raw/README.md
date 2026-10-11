# M6 raw material (AGC 3.2.4 vs refrel3 at N = 5 / 50 / 100 / 200 / 558)

Copy of aceapex `research/bench2610/m6` (tag `panvram-1.0-sources`) plus the query harness `harness/agcbench.cpp`
(aceapex `research/agc_vs_refrel3/agcbench.cpp` @ a6a27e2, not in the tag).

- Protocol: `PROTOCOL_M6.md` (v1.1; sha256 in `PROTOCOL_M6.sha256`).
- AGC: 3.2.4, commit e67e3fc, built on ace-core with g++ 11.4 `-O3 -march=native -std=c++20` (platform `native`).
- Archives: per N one `agc create -d -k 31 -l 20 -s 60000 -b 50 -t 16 -o agcN.agc <T2T> <N FIFOs>`
  (`run_n.sh`; `logs/N*/create.time`, `create.err`, `archive.txt`, `integrity.txt`).
- Queries: `agcbench req <archive> <requests> <truth> <threads> <warm-ups> <runs>` (windows of 4 096, 1 Mb regions) and
  `agcbench sample` (whole samples), one handle per thread (`run_n.sh`; raw timings `logs/N*/agc_*_t1.log`, `agc_*_t16.log`;
  refrel3 rows `logs/N*/rr3_*`); counters `logs/N*/counters.txt`, `cnt_*.tsv`; profile `logs/N*/perf_*`; silence gate
  `logs/N*/silence_*.txt`; driver log `logs/run_all.log`.
- Local paths: in 111 text files the home directory `/home/aeterna` is replaced by `~`; nothing else is changed.
  `PROVENANCE.tsv` lists every file with the SHA-256 of the original (as in the aceapex tag) and of this copy.
