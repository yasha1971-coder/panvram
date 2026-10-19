# Independent baseline: Cas-OFFinder

Cas-OFFinder 3.0.0, github.com/snugel/cas-offinder at commit `0a9ac005859641ca215add6c6dbed147ae80f7c0`, used as an
outside check that the reference-transfer plus residual decomposition returns the same site set as a plain scan.

No OpenCL runtime is installed on ace-core, so the stack was assembled from source and from pip:

- OpenCL headers: KhronosGroup/OpenCL-Headers (shallow clone, header only).
- ICD loader: KhronosGroup/OpenCL-ICD-Loader, built with CMake against those headers, giving `libOpenCL.so`.
- Implementation: `pocl-binary-distribution` from pip, which ships `libpocl` and an ICD entry; the ICD file was
  rewritten with the absolute path of that library and selected through `OCL_ICD_VENDORS`.
- Device: CPU (`cas-offinder <input> C <output>`); no GPU was used.

Build of Cas-OFFinder:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DOpenCL_INCLUDE_DIR=<OpenCL-Headers> -DOpenCL_LIBRARY=<OpenCL-ICD-Loader/build/libOpenCL.so>
    cmake --build build -j8

Input file for the reference run is `input_reference.txt`: the genome directory, the pattern
`NNNNNNNNNNNNNNNNNNNNNGG`, then one line per guide as `<spacer>NNN 6`. Cohort runs used the same pattern and guides
with the assembly FASTA in place of the reference; `run_casoff_asm.sh` is the driver. It decompresses one assembly,
searches it and deletes it again, because 50 assemblies uncompressed need about 150 GB while 22 GB were free; the
decompression time is reported separately in `times_n50_per_assembly.txt`.

Both tools fold lower case to upper (Cas-OFFinder in `read_fasta.cpp`), and the T2T FASTA does contain soft-masked
regions, so that rule is exercised by the reference comparison.

## Result

Sets compared as tuples: reference run `(guide, contig, position, strand, substitutions)`, cohort runs
`(assembly, contig, position, strand, substitutions, guide)`. Positions are zero-based inside the named contig; the
tuples are sorted and hashed, digests in `set_digests.json`.

| data | sites | sets identical | SHA-256 of the sorted set |
|---|---:|---|---|
| T2T reference, 4 guides | 121 079 | yes | a8bc09fdef6484a2… |
| 5 assemblies, 4 guides | 587 044 | yes | 69bfcf4501ff599b… |
| 50 assemblies, 4 guides | 5 916 825 | yes | db6352e200ca6172… |

Guides are the four of `aryeelab/guideseq`; PAM NGG; up to 6 substitutions in the 20-base spacer.

## Time

| run | Cas-OFFinder | this decomposition |
|---|---:|---:|
| T2T reference | 19.3 s | 10.7 s total, of which 0.8 s is the scan |
| 50 assemblies | 929 s search plus 581 s decompression | 23.7 s |

The two do different amounts of input work: Cas-OFFinder reads uncompressed FASTA, the decomposition reads the
archives and decodes them. This is a check that the answers agree, not a benchmark of either tool.
