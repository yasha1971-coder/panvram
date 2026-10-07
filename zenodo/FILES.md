# Zenodo dataset record - files (draft; not uploaded)

Metadata: `zenodo/dataset.zenodo.json`. Nine files, 20 429 332 480 B of tars + 0.43 MB (limit 50 GB per record). Built 2026-10-07 on ace-core in
~/outgoing/zenodo (not uploaded); every member read back: 558/558 SHA-256 == MANIFEST.tsv in both tars; determinism: a
second build of cohort558_q16k.tar gave the same SHA-256.

| # | file | size | made by |
|---:|---|---:|---|
| 1 | cohort558_q4k.tar | 12 330 864 640 B; sha256 d0a4a0baab5bba5e573c8f2ae5dd18f2a66f6ac222c29e9809e38ad99b26ecdc | `scripts/make_dataset_tars.sh` |
| 2 | cohort558_q4k.tar.sha256 | 84 B | same |
| 3 | cohort558_q16k.tar | 8 098 467 840 B; sha256 3b8c497e9193984724128e0fd18df02d131a05c341df41df44fe3ed89a8114d5 | same |
| 4 | cohort558_q16k.tar.sha256 | 85 B | same |
| 5 | MANIFEST.tsv | 196 949 B | repository |
| 6 | SHA256SUMS.q4k | 56 583 B | repository |
| 7 | SHA256SUMS.q16k | 57 141 B | repository |
| 8 | DATASET.md | 118 030 B | repository |
| 9 | FORMAT.md | (format of the archives) | repository |

Not included: the reference (T2T-CHM13v2.0, NCBI GCA_009914755.4; md5 in DATASET.md), the source FASTA (HPRC indexes,
URL and hash per row of MANIFEST.tsv), the software (its own record, `.zenodo.json`).
Tars: archives at the tar root, deterministic (`--sort=name --owner=0 --group=0 --numeric-owner --mtime='2026-10-04 00:00Z'`);
each checked by reading every member back and comparing its SHA-256 with MANIFEST.tsv.
