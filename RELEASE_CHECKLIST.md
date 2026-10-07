# Release checklist - panvram 1.0.0 (prepared locally, not pushed, not public)

Checked 2026-10-07 on commit-to-be (git ls-files of the working tree):

- [x] No secrets: no API keys / tokens / private keys / passwords (grep: ghp_, github_pat_, AKIA..., PRIVATE KEY, api key, password, token=) - 0 hits ("token" elsewhere = ACGT tokens of the API).
- [x] No personal gmail: 0 hits for gmail / yasha1971@. Contact e-mail: aceapex.contour@outlook.com (CITATION.cff, .zenodo.json).
- [x] No local paths: 0 hits for /home/aeterna and "aeterna" in tracked files (Colab logs show /content paths; the clean-room round-2 report shows the implementer's own sandbox /home/claude/cr2).
- [x] License MIT (full text, Yakiv Shavidze); xxHash BSD-2-Clause notice kept in csrc/xxhash.h.
- [x] CITATION.cff: Yakiv Shavidze, ORCID 0009-0008-3622-3448.
- [x] HPRC data use in README (no sale, no re-identification, acknowledgement PRJNA730823 + NHGRI); full quotes in DATASET.md.
- [x] README numbers only from files in logs/; Limitations: size 13.04 vs AGC 5.40 MB per assembly (558), 16 GB by memory measurement only, byte determinism with the gcc 11.4 recipe, haplotype coordinates only.
- [x] Public notebook notebooks/quickstart_558.ipynb without Google Drive (archives by MANIFEST.tsv, SHA-256 checked).
- [ ] Before publishing (user's decision): dataset record and package URL placeholders in quickstart_558.ipynb and gate558_colab.ipynb, DOI in .zenodo.json / CITATION.cff; the aceapex links in README / EVIDENCE.md point to branch refrel (local commits not pushed).
