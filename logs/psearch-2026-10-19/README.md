# Pattern search over the refrel3 v1 cohort: reference transfer and residual windows (2026-10-19)

A refrel3 v1 archive stores an assembly as copies from the decoded reference plus its own literal bytes, so a window
of 23 bases either lies inside one reference-contiguous copy or it does not. The first kind can be answered by scanning
the reference once and carrying the answer through the copy coordinates ("transfer"); only the second kind has to be
looked at in the assembly itself ("residual"). This directory records that decomposition for the 558-assembly cohort
and the sequences it turns up. It is a measurement of the representation, not a biological study.

Run by the author on ace-core (AMD EPYC 4344P, Ubuntu 22.04.5, g++ 11.4.0), 16 threads, behind a silence gate.
Search rule: canonical domain = upper-case bases of the contig records (lower case folded up); a window of 23 bases is
considered only if all 23 are A/C/G/T; windows never cross contigs; both strands; for guide matching the PAM is NGG
and up to 6 substitutions are allowed in the 20-base spacer, the PAM itself exact.

## Correctness of the decomposition

The decomposition was checked against searches that do not use it:

| check | scale | result |
|---|---|---|
| transfer + residual vs every window of the fully decoded archive | 558 assemblies, 2 993 arbitrary 23-mers | identical, 3 445 753 630 hits, digest 686af4ef… |
| transfer + residual vs every window of the source FASTA (`.fa.gz`) | 50 assemblies | identical, 309 628 026 hits |
| transfer + residual vs Cas-OFFinder 3.0.0 | T2T reference, 4 guides, ≤6 substitutions | identical sets, 121 079 sites |
| transfer + residual vs Cas-OFFinder 3.0.0 | 50 assemblies, 4 guides, ≤6 substitutions | identical sets, 5 916 825 sites |

Cas-OFFinder (github.com/snugel/cas-offinder @0a9ac005, built against the Khronos ICD loader and pocl on CPU) folds
case the same way. Three negative controls were also run and each was detected: a broken reverse-complement transfer,
a shifted context at a block boundary, and a shifted self-copy source.

Wall time at 50 assemblies, same guides and threshold: 23.7 s for the decomposition against 929 s of search plus 581 s
of decompression for Cas-OFFinder, which reads uncompressed FASTA. The two tools do different amounts of I/O and this
is not a benchmark of either.

## Guides

12 SpCas9 guides with an NGG PAM, taken unchanged from two public repositories (`guides.tsv`):

| source | guides |
|---|---|
| aryeelab/guideseq @997b89242dec2ac4a1a63273025aad371e448b8b | EMX1 (`test/test_manifest.yaml`), VEGFA site 1–3 (`test/large_test/reference_output`) |
| tsailabSJ/genetic_variation_guideseq @31f8ffd04ce6b04a52eb8d33466735dc59ea8299 | CTLA4.s9, AAVS1.s14, CXCR4.s8, LAG3.s9, TRAC.s1, CCR5.s8, AAVS1.s2, TRAC_site_2 (`All_supplemental_tables_260728_AF.xlsx`, sheet `Supp_Table1|gRNA`) |

The CHANGE-seq panel of Lazzarotto et al. 2020 is not in `tsailabSJ/changeseq`: its test configuration carries two
guides that duplicate EMX1 and VEGFA site 1. The full spacer list is in that paper's supplementary material.

## Files

- `guides.tsv` — guide, spacer, PAM class, source repository and commit.
- `per_guide.tsv` — per guide: all sites, sites carried from the reference scan, residual sites and their share,
  unique residual 23-mers, how many of those are absent from T2T or present elsewhere in it, the absent ones split by
  substitution count, and the count of absent sequences at ≤2 substitutions that are not simple repeats.
- `absent_from_T2T_le3_mismatches.tsv` — every unique residual 23-mer absent from T2T at ≤3 substitutions: sequence,
  substitutions, occurrences, how many assemblies carry it, repeat class with its two measures, and one example
  location (assembly, contig, position, strand).
- `summary.json` — the cohort totals quoted below.

Positions are zero-based offsets inside the named contig of the named assembly. A sequence is reported in canonical
form, the lesser of the window and its reverse complement, so the two strands are one entry.

## Repeat classification rule

A 23-mer is called a **simple repeat** if either of the following holds, otherwise a **unique sequence**:

- the Shannon entropy of its single-base composition is below 1.5 bits (the maximum for four bases is 2.0), or
- for some period p in {1, 2, 3}, at least 90 % of the positions equal the base p positions ahead.

Both measures are recorded per sequence in `absent_from_T2T_le3_mismatches.tsv`, so a different threshold can be
applied to the same data. The rule is a coarse filter chosen for this table; it is not a standard definition.

## Cohort totals (558 assemblies, 12 guides, NGG, ≤6 substitutions)

| quantity | value |
|---|---:|
| sites found | 145 014 299 |
| of them carried from the single reference scan | 134 144 968 |
| of them residual | 10 869 331 |
| unique residual 23-mers | 118 132 |
| absent from T2T | 96 616 |
| present in T2T at another location | 21 516 |
| absent and classed simple repeat | 38 583 |
| absent and classed unique sequence | 58 033 |

Residual share per guide runs from 2.6 % (TRAC.s1) to 17.1 % (VEGFA site 3). The highest share belongs to the guide
whose spacer is itself a CA/CT repeat.

## Observations

These are observations on one cohort of 558 HPRC haplotype assemblies against one reference (T2T-CHM13), under the
search rule stated above. They are not a conclusion about the safety of any editing experiment, not a prediction of
cleavage, and not an evaluation of any search tool. Whether a sequence is cut in a cell depends on chromatin, delivery,
expression and assay sensitivity, none of which is measured here.

- At 0 substitutions no site is absent from the reference: for all 12 guides every exact spacer-plus-NGG match seen in
  an assembly is also present in T2T.
- At ≤2 substitutions, 23 unique sequences are absent from T2T; 9 of them are classed as unique sequence rather than
  simple repeat, and they belong to 7 of the 12 guides. Three of those nine occur in 115 to 189 of the 558 assemblies,
  the rest in 2 to 36.
- About 40 % of all sequences absent from T2T (38 583 of 96 616) are classed as simple repeats, so a repeat filter
  changes the size of such a list substantially.
- 21 516 unique residual sequences do occur in T2T, but at a different location, so a reference-only scan reports them
  with coordinates that do not apply to the assembly.

What this does not establish: the cohort is 558 assemblies and not a population sample; a sequence absent from T2T may
be present in references other than T2T; the counts depend on the 23-base window, the NGG PAM and the substitution
threshold; and no experimental off-target data was used anywhere in this directory.
