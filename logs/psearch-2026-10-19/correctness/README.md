# Per-assembly digests

Every file lists one assembly per line: its name, then for each search mode the number of hits and the SHA-256 of its
sorted hit list. A hit is `(contig index, position in contig, pattern or guide id, strand)`, and for CRISPR also the
substitution count; the list is sorted and hashed as fixed-width records. Two modes agree for an assembly exactly when
the two digests are equal, which is how the `COMPARE` verdicts in the run logs are decided.

| file | search | assemblies | modes |
|---|---|---:|---|
| `exact_n50_fa_v1_parse_parse2.tsv` | 2 993 arbitrary 23-mers | 50 | source FASTA, full decode, transfer+residual, on-demand resolution |
| `exact_n558_v1_parse_parse2.tsv` | 2 993 arbitrary 23-mers | 558 | full decode, transfer+residual, on-demand resolution |
| `crispr_n50_fa_v1_parse.tsv` | 120 guides, NGG, <= 6 substitutions | 50 | source FASTA, full decode, transfer+residual |
| `crispr_n558_v1_parse.tsv` | 120 guides, NGG, <= 6 substitutions | 558 | full decode, transfer+residual |
| `recheck_exact_n50_current_binary.tsv` | 2 993 arbitrary 23-mers | 50 | source FASTA, transfer+residual, rerun with the committed build |
| `recheck_crispr_n5_current_binary.tsv` | 120 guides | 5 | full decode, transfer+residual, rerun with the committed build |
| `patterns_2993.tsv` | the patterns themselves | — | id, class, sequence, where it was sampled from |

In every file the digests of all modes of a line are equal: no assembly differs between a search that uses the
decomposition and one that does not. The source FASTA columns are the ones that depend on nothing in this repository.
The two recheck files exist because the archived passes used earlier builds; see `../REPRODUCIBILITY.md`.
