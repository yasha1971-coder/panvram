# Security

panvram reads untrusted archive files. Every archive is checked before any decode (FORMAT.md section 5: header XXH3,
reference SHA-256 and size, section sizes, model tables, block table); decoding is bounded (every copy inside the block
and the reference, every rANS stream consumed exactly; the GPU kernel's waits are bounded and report a status instead of
hanging). Block XXH3 is checked on the CPU path (`fasta`, `windows(..., verify)`), not in the GPU kernel.

Reporting: please report a vulnerability (a crash, hang, out-of-bounds access or silently wrong output on some archive)
privately through the repository's GitHub security advisory ("Report a vulnerability"), with the archive or a script
that produces it. Do not open a public issue for it.

Supported versions: the latest release.
