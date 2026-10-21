# refrel v2 prototype at N = 558: window law, round trip, negative controls (2026-10-21)

Logs behind [`research/refrel_v2/RESULTS_V2.md`](../../research/refrel_v2/RESULTS_V2.md). The v2 container is a
research prototype and is **not frozen**: it is unstable and may stop being readable. refrel3 v1 is untouched. The
**encoder is not published**, so the compression cannot be repeated from this repository; decoding, the window law and
the controls below can.

Host of every run: ace-core, AMD EPYC 4344P (8 cores / 16 threads), Ubuntu 22.04.5, kernel 5.15.0-163, g++ 11.4.0.
Peak RSS is the `Maximum resident set size` of `/usr/bin/time -v`, in kbytes.

## `law-2026-10-21-rv2dec/` - the headline run

The window law (model B, hw-apex protocol v1.1) on the 558 v2 archives with the decoder of `research/refrel_v2`
(binary `6e1f91b7...`, `build.sh`): D_Q 3.399679 Gbases/s, c0 3.266 us, errors -7.52 / -5.02 / +19.77 % at
1 / 8 / 64 KiB, **PASS**, peak RSS 36.6 GB, wall 1 h 48 min. Requests and truth are the files of the published v1
verdict (`logs/verdict-n558-v1`), same SHA-256. One thread on CPU 3, silence gate before the run
(`silence_*.txt`), 10 000 requests per W, every answer's SHA-256 == truth (bad 0), every full decode's XXH3 == truth.
`raw_*.tsv.gz`: W, request id, nanoseconds, 40 000 rows. Inputs and the binary: `*_inputs.sha256`.

Command: `bash research/refrel_v2/run_law.sh <rv2dec> <t2t.fa> <v1 dir> <anchors.txt> <rv2 dir> <names.txt>
<requests.tsv> <truth.tsv> <out dir>`.

## `law-2026-10-09-v2/` and `law-2026-10-10-v3/` - earlier prototype runs

The same archives measured with the prototype binary (the one that still contained the encoder): version 2 of the
block table (`67c0ccb7...`: D_Q 3.423414, c0 3.159 us, -5.44 / -2.11 / +18.54 %, PASS) and version 3
(`7b4152f1...`: D_Q 3.432540, c0 3.542 us, -4.27 / -2.82 / +19.32 %, PASS). The difference between these and the
2026-10-21 run is tabulated in RESULTS_V2.md; its cause is not measured.

## `roundtrip-n558/` - every assembly decoded in full

`n558_k8_m3_judge.tsv`: 558 / 558 PASS. Each assembly decoded from its archive (the 8 anchors from their v1 q16k
archives, the other 550 from v2) to FASTA; FASTA bytes and FASTA XXH3 == the source (`MANIFEST.tsv`) and SHA-256 of
the upper-case bases == the truth of the AGC comparison. Wall 5 768 s, peak RSS 40.2 GB. `n558_k8_m3_sizes.txt`:
2 972 864 737 B over 558 = 5.3277 MB per assembly, anchors counted as members at their v1 size without block hashes.

## `archives/` - what to fetch and check

- `anchors.txt` - the 8 anchor assemblies, in the order the decoder expects.
- `SHA256SUMS.anchors` - their v1 q16k archives, which are files of the published dataset (DATASET.md, Zenodo record
  10.5281/zenodo.23232317; "Rebuilding a tar from its parts"). The files themselves are not in this repository.
- `SHA256SUMS.rv2` - the 550 v2 archives of the run (2.86 GB in total). Not in this repository either; five of them
  are a release asset, see below.

## `subset5/` - the five assemblies the controls run on

`names5.txt` (five year-1 assemblies), `SHA256SUMS` (their five v2 archives, 22 752 605 B in total),
`req5.tsv.gz` and `truth5.tsv.gz`. The requests come from the frozen generator of the hw-apex job
(`tools/window_cost/window_requests.py`, seed 20261003, 10 000 per W); the truth comes from
`tools/window_cost window_cost truth`, that is from the **v1** archives with every v1 check - it does not use any v2
code, so it is an independent oracle for the v2 decoder.

## `negative_controls/` - five deliberate faults, each caught

Each fault is one line of decoder code behind a compile-time flag of `research/refrel_v2/build.sh`; without the flag
no fault code is compiled (the clean binary is `6e1f91b7...` in every run, `clean_build.txt`). Driver:
`research/refrel_v2/run_controls.sh`. For each variant the five assemblies are decoded in full (`rv2dec judge`,
FASTA bytes + XXH3 against `MANIFEST.tsv`) and the window law is measured (`rv2dec law`, every answer's SHA-256 and
every full decode's XXH3 against the subset truth).

| variant | what it breaks | judge | window law | binary |
|---|---|---|---|---|
| clean | nothing | **pass 5 / fail 0** | **PASS**, bad 0 | `6e1f91b7...` |
| `break_lift` | the lifted diagonal is one base off | fail 5 | refused / exit 4 | `a1229d6b...` |
| `break_abs` | an absolute copy into an anchor is one base off | fail 5 | refused / exit 4 | `4a84d9a5...` |
| `break_bt3` | the per-target block-table prediction is one base off | fail 5 | refused at open: `block table target` | `fb28b20e...` |
| `break_ref` | an absolute copy into the reference is one base off | fail 5 | refused / exit 4 | `73a9aada...` |
| `break_self` | a self copy resolves one base earlier | fail 5 | refused / exit 4 | `3680a895...` |

Totals per variant are `SUMMARY.tsv`; the per-assembly lines and the timings are the `*_judge.txt` and `*_law.txt`
files (the per-request and per-decode lines are dropped, the result lines kept).

**What this proves and what it does not.** It proves that the judge notices a decoder that is wrong by one base on any
of the five copy paths. It does **not** prove that the judge would notice a corrupted archive: producing a corrupted
v2 archive needs the encoder, which is not published. A sixth fault (`break_lower`, lower-case runs applied one base
late) was written and then dropped: the five assemblies of the subset carry no lower-case runs at all (`rv2dec ledger`
reports 8 bytes of case-run data each, the count only), so the fault had nothing to shift and the judge had nothing to
catch. Keeping it would have meant publishing a control that passes for a reason unrelated to the judge.
