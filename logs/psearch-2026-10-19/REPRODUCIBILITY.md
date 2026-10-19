# Which binary produced which log

`bash tools/psearch/build.sh` is deterministic on this host: three consecutive builds give the same bytes, and a build
from the sources as committed here reproduces the binaries in use at the end of the work:

| tool | SHA-256 of a build from this branch |
|---|---|
| psearch | b7176b0ff76a8c589ae41fd7f059f968db434e5a34e9783a355c15df4802f0c8 |
| crispr | 89f5534977b8afa1458be734bcd6f709e321b4c9f11087d8309f11f5e99ed14e |

The tools grew while the work went on, so earlier passes ran earlier builds. Each run receipt records the binary it
used; collected here:

| run | tool | SHA-256 recorded in the run | equals the build above |
|---|---|---|---|
| exact search, 558 assemblies, origin split | psearch | b7176b0f… | yes |
| exact search, 558 assemblies, `v1,parse,parse2` | psearch | cacb31ff… | no |
| exact search, 50 assemblies, `fa,v1,parse,parse2` | psearch | 0d91e213… | no |
| CRISPR, 50 and 558 assemblies, 120 guides | crispr | 4a98afa6… | no |

Because three of those differ, the passes were repeated with the current binaries and the per-assembly SHA-256 of the
sorted hit lists compared with the archived ones:

| recheck | result |
|---|---|
| psearch, 50 assemblies, modes `fa` and `parse`, 2 993 patterns | 50 assemblies, 0 digests differ from `correctness/exact_n50_fa_v1_parse_parse2.tsv` |
| crispr, 5 assemblies, 120 guides, modes `v1` and `parse` | 5 assemblies, 0 digests differ from `correctness/crispr_n50_fa_v1_parse.tsv` |

Receipts of the rechecks are `correctness/recheck_exact_n50_current_binary.tsv` and
`correctness/recheck_crispr_n5_current_binary.tsv`. So the earlier builds and the committed one return the same hits;
what changed between them was added commands and counters, not the search.

The intermediate source states are not kept, so a build of this branch will match the last column of the first table
and not the three earlier hashes. The digest rechecks are what ties the archived numbers to the committed source.
