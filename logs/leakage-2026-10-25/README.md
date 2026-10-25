# Leakage by locus: what a train / test split of a pangenome actually separates (2026-10-25)

A training window and a test window that answer to the same stretch of the reference are not independent samples. If
a model is trained on windows of some assemblies and tested on windows of others, a test window can still cover the
very locus the model has already read in another haplotype. This directory measures how often that happens, and what
a different way of splitting changes.

The measurement is possible at all because a refrel3 v1 archive already carries the alignment: every copy records
where in the reference it came from, so the reference interval of a window is a read of the archive, not a new
alignment. The tool is `scripts/locate.cpp`, built against the v1 reader; `tilemap` tiles every contig into windows
and writes each window's interval, `leak` does the whole measurement in the tool for a cohort too large to dump.

A window is used only if **both** of its ends have a reference coordinate and the hull of the two is within 10 % of
the window length. A window that crosses a rearrangement or a long insertion is not one locus, and it is counted
separately rather than pretended into one. Overlap is the share of a test window's reference interval that the union
of the training windows' intervals covers.

Run by the author on ace-core (AMD EPYC 4344P, Ubuntu 22.04.5, g++ 11.4.0). HPRC assemblies: BioProject PRJNA730823;
their use here implies no endorsement by the HPRC.

## Four assemblies, W = 8 192 (`n4/`)

`y1_HG00438.1`, `y1_HG00438.2`, `y1_HG00621.1`, `y1_HG00621.2` - the cohort of the N = 4 window-law verdict.
1 463 063 windows mapped in 0.8 s on 8 threads, **1 876 206 windows/s** (two `locate` calls per window, so 3.75 M
coordinate lookups per second).

| | windows | share |
|---|---:|---:|
| tiled | 1 463 063 | 100 % |
| both ends have a reference coordinate | 1 457 730 | 99.635 % |
| of those, hull within 10 % of W - used below | 1 354 262 | 92.90 % |
| an end falls in a literal, no coordinate exists | 5 333 | 0.365 % |
| hull far from W: a rearrangement or a long insertion | 103 468 | 7.10 % |

| split | train / test windows | train covers T2T | overlap >= 50 % | overlap >= 90 % | median overlap | zero overlap |
|---|---|---:|---:|---:|---:|---:|
| random, by sample (train HG00438.1 + HG00621.1) | 669 173 / 685 089 | 92.4 % | **99.305 %** | **98.954 %** | 100.00 % | 0.429 % |
| by individual (both haplotypes of HG00438 in train) | 685 157 / 669 105 | 91.8 % | **98.966 %** | **98.651 %** | 100.00 % | 0.804 % |
| by locus, T2T in 10 Mb chunks | 674 871 / 679 391 | 46.4 % | **0.024 %** | **0.002 %** | 0.00 % | 99.937 % |
| by locus, T2T in 1 Mb chunks | 681 023 / 673 239 | 47.0 % | **0.219 %** | **0.014 %** | 0.00 % | 99.401 % |

Splitting by sample separates almost nothing: **99.3 % of test windows are covered more than half by some training
window, 99.0 % more than nine tenths, and the median overlap is 100 %**. Holding out a whole individual, both
haplotypes, changes it by a third of a percent. Splitting by locus removes it: 0.024 % at 10 Mb chunks. The residue
is the rearranged windows - cut the reference finer, into 1 Mb, and it rises to 0.219 %, because there are more
boundaries to cross.

What this does **not** say: nothing here is a statement about any model's accuracy. It is a property of the sampling,
measured on one cohort against one reference, and the overlap criterion is the interval of a window, not its sequence
identity - two windows over the same locus in different haplotypes are not the same string.

## All 558 assemblies (`n558/leak558.log`)

234 individuals behind the 558 haplotype assemblies; the split by individual puts 117 of them, 272 assemblies, in
train. The locus split cuts T2T into 312 chunks of 10 Mb and sends 149 of them to train. Windows tiled: 205 093 303 at
W = 8 192 and 51 243 747 at W = 32 768 per pass; **3.3 to 4.8 million windows per second** on 16 threads, so the
whole cohort is mapped in 15 to 85 s per run.

| W | split | train covers T2T | test windows | overlap >= 50 % | overlap >= 90 % | mean overlap | zero overlap |
|---:|---|---:|---:|---:|---:|---:|---:|
| 8 192 | by individual | 98.8 % | 96 774 185 | **99.9955 %** | **99.9901 %** | 99.9946 % | 0.0019 % |
| 8 192 | by locus, 10 Mb | 47.4 % | 100 029 548 | **0.0414 %** | **0.0085 %** | 0.0414 % | 99.9177 % |
| 32 768 | by individual | 98.3 % | 24 004 831 | **99.9932 %** | **99.9878 %** | 99.9926 % | 0.0034 % |
| 32 768 | by locus, 10 Mb | 47.3 % | 24 816 467 | **0.1634 %** | **0.0330 %** | 0.1638 % | 99.6721 % |

Windows dropped as not being one locus: at W = 8 192, 15 199 486 of 205 093 303 for a hull far from W (7.4 %) and
1 067 986 for an end in a literal (0.5 %); at W = 32 768, 4 078 603 (8.0 %) and 327 336 (0.6 %).

At cohort scale the picture is sharper than on four assemblies, not softer. Holding out **117 of 234 people** leaves
**99.9955 %** of the test windows covered by a training window more than half, and 99.9901 % more than nine tenths:
with 272 assemblies the training side covers 98.8 % of the reference, so there is almost nowhere left for a test
window to be new. The locus split keeps leakage at 0.04 % for 8 KiB windows. It rises to 0.16 % at 32 KiB, four times
higher, for the reason the four-assembly table already showed: a longer window is likelier to reach across a chunk
boundary or to sit on a rearrangement that maps into a training chunk.

So the choice that matters for a pangenome training set is not how many individuals are held out. It is whether the
split is made in reference coordinates at all.

## Files

- `scripts/locate.cpp` - the tool. `locate` maps single positions, `tilemap` writes the window table of a small
  cohort, `leak` computes the whole measurement in process for a large one. Build: see the header of the file.
- `scripts/leak_n4.py` - the split and overlap arithmetic for the four-assembly table above.
- `n4/leak_n4.log` - its output, the table above verbatim.
- `n4/tiles8k.head.tsv` - the first rows of the window table, for the format.
- `n4/names4.txt` - the four assemblies.
- `n558/leak558.log` - the four runs on the whole cohort, verbatim, and the accessor self-test that ran before them.
- `n558/names558.txt` - the cohort order.
- `scripts/run_leak558.sh` - the driver. It waits for the window-law measurement of `research/refrel_v2` to finish
  before it starts: that measurement times single requests on a pinned core, and load from these runs would corrupt
  it.
- `inputs.sha256` - the tool, the reference and both name lists.
