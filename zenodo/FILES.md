# Zenodo dataset record - files (draft; not uploaded)

Metadata: `zenodo/dataset.zenodo.json`. Nine files, about 20.43 GB in total (limit 50 GB per record):

| # | file | size | made by |
|---:|---|---:|---|
| 1 | cohort558_q4k.tar | ~12.33 GB (12 330 370 123 B of archives + SHA256SUMS + tar headers) | `scripts/make_dataset_tars.sh` |
| 2 | cohort558_q4k.tar.sha256 | 84 B | same |
| 3 | cohort558_q16k.tar | ~8.10 GB (8 097 986 851 B of archives + SHA256SUMS + tar headers) | same |
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
