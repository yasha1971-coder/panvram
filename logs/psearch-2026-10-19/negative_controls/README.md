# Deliberate errors

Four single-line faults, each in a different part of the decomposition, run on 5 assemblies so that the comparison has
something to catch. All four are caught: the per-assembly digests differ for every assembly and the pass reports FAIL.
`00_clean_PASS.log` is the same comparison with nothing broken.

| log | what is broken | verdict |
|---|---|---|
| `00_clean_PASS.log` | nothing | PASS, 0 of 5 assemblies differ |
| `01_broken_reverse_complement_transfer.log` | carrying a hit through a reverse-complement run keeps the strand instead of flipping it | FAIL, 5 of 5 |
| `02_shifted_block_boundary_context.log` | the context carried across a 16 KiB block boundary is shifted by one base | FAIL, 5 of 5 |
| `03_shifted_self_copy_source.log` | a self copy resolves its source one base off | FAIL, 5 of 5 |
| `04_shifted_transferred_mismatch_count.log` | a carried CRISPR site with three substitutions is reported as four | FAIL, 5 of 5 |

Each log keeps its first six differing hits, which name the pattern class or guide involved. The switches are
environment variables listed in `tools/psearch/README.md`; they are off unless set.
