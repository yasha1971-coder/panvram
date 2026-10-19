# psearch and crispr: pattern search over a refrel3 v1 cohort

Two tools that search a cohort of refrel3 v1 archives three ways and compare the answers. They include only
`csrc/pv_v1.h` from this repository, the public v1 reader; nothing from an encoder or a research prototype.

    bash tools/psearch/build.sh [out dir]     # needs g++ >= 11, libzstd-dev, zlib

Common rule for both tools: the canonical domain is the upper-case bases of the contig records, lower case folded up;
a window of 23 bases is used only if all 23 are A/C/G/T; windows never cross contigs; both strands are reported.

Modes, in order of what they trust:

- `fa` — every window of the source FASTA (`.fa.gz`), independent of the archives.
- `v1` — every window of the fully decoded archive.
- `parse` — the reference is scanned once; per assembly only the residual windows are examined (literals, self copies,
  windows straddling the boundary of a reference-contiguous run), and the rest of the answers are carried over through
  the copy coordinates.
- `parse2` — as `parse` but residual bytes are resolved on demand instead of executing the copies; byte-exact and
  slower, kept as a measured negative result.

Passing several modes runs each and compares them per assembly by the SHA-256 of its sorted hit list, printing one
`COMPARE` line per pair and up to 20 differing hits.

## psearch: exact 23-mers

    psearch gen <reference.fa> <archive dir> <names.txt> <seed> <patterns.tsv>
    psearch run <reference.fa> <archive dir> <names.txt> <patterns.tsv> <modes> <threads> [fasta dir]

`gen` samples patterns from the reference and from the archives in seven classes: reference, literal, run boundary,
short insert, self copy, contig end, block edge. `run` searches them; with mode `origin` it also splits each hit into
carried-from-the-reference and residual, and counts unique +/- 20 base contexts of residual hits.

    # verify the decomposition against the source FASTA on 50 assemblies
    psearch run t2t.fa cohort/v1 names50.txt patterns.tsv fa,v1,parse,parse2 16 fasta_dir
    # origin split, writing the per-pattern table
    PSEARCH_ORIGIN_OUT=origin.tsv psearch run t2t.fa cohort/v1 names558.txt patterns.tsv origin 16

## crispr: spacer plus PAM with substitutions

A site at window start p is `spacer` with at most M substitutions followed by any base and `GG` on the plus strand,
or `CC`, any base, and the reverse complement of the spacer on the minus strand. The PAM is matched exactly.

    crispr guides   <reference.fa> <archive dir> <names.txt> <patterns.tsv> <seed> <n_ref> <n_asm> <guides.tsv>
    crispr run      <reference.fa> <archive dir> <names.txt> <guides.tsv> <M> <modes> <threads> [fasta dir]
    crispr refscan  <reference.fa> <guides.tsv> <M> <threads> <out.tsv>
    crispr sites    <reference.fa> <archive dir> <names.txt> <guides.tsv> <M> <threads> <out.tsv>
    crispr gorigin  <reference.fa> <archive dir> <names.txt> <guides.tsv> <M> <threads>
    crispr residual <reference.fa> <archive dir> <names.txt> <guides.tsv> <M> <threads> <out.tsv>
    crispr kmerscan <reference.fa> <sequences.txt> <threads> <out.tsv>

`refscan` dumps the sites of the reference alone, `sites` dumps every site of every assembly with contig names (the
form used to compare with another tool), `gorigin` prints the per-guide split into carried and residual with the
substitution histogram, `residual` collects the unique canonical 23-mers of residual sites with one example location
each, and `kmerscan` counts how often given 23-mers occur in the reference, either strand.

    # compare with an external tool
    crispr sites t2t.fa cohort/v1 names50.txt guides.tsv 6 16 ours.tsv
    # per-guide origin split and the not-in-reference filter
    crispr gorigin  t2t.fa cohort/v1 names558.txt guides.tsv 6 16
    crispr residual t2t.fa cohort/v1 names558.txt guides.tsv 6 16 residual.tsv
    cut -f2 residual.tsv | tail -n +2 | sort -u > seqs.txt
    crispr kmerscan t2t.fa seqs.txt 16 refcount.tsv

## Deliberate-error switches

Set only to check that a comparison would notice. Each makes the corresponding comparison fail; logs of all four are
in `logs/psearch-2026-10-19/negative_controls/`.

| variable | what it breaks |
|---|---|
| `PSEARCH_BREAK_RC` | the transfer of a reverse-complement run keeps the strand instead of flipping it |
| `PSEARCH_BREAK_CTX` | the context carried across a block boundary is shifted by one base (`parse2`) |
| `PSEARCH_BREAK_SELF` | a self copy resolves its source one base off (`parse2`) |
| `CRISPR_BREAK_MM` | a carried site with three substitutions is reported as four |

Other variables: `PSEARCH_PER_ASM` and `CRISPR_PER_ASM` write the per-assembly digests, `PSEARCH_ORIGIN_OUT` the
per-pattern origin table, `PSEARCH_CHAIN_STATS` is not needed (chain counters are always reported for `parse2`).
