# Linux CPU wheels (TestPyPI only)

The `pypi` branch builds x86_64 manylinux_2_28 wheels for CPython 3.10 to
3.14. The native extension is compiled against torch 2.14.1+cpu; wheel
metadata pins the release, torch==2.14.1, because its C++ ABI is not stable
between releases. The pin carries no local label, so torch 2.14.1 from PyPI
(the default CUDA build) or 2.14.1+cpu from download.pytorch.org satisfies it. No GPU or CUDA toolkit is needed. These wheels do not contain the
optional CUDA kernel. Building from source with `PANVRAM_CUDA=1` remains
available; CUDA wheels are outside this workflow's scope.

Install runtime dependencies first (TestPyPI is not a dependency mirror):

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install torch==2.14.1+cpu --index-url https://download.pytorch.org/whl/cpu
pip install -i https://test.pypi.org/simple panvram
```

Or in one step, with dependencies from PyPI (this installs torch 2.14.1 from PyPI if another version is present):

```sh
pip install -i https://test.pypi.org/simple --extra-index-url https://pypi.org/simple panvram
```

The CPU extension loads torch's libraries from the pinned runtime package.
auditwheel bundles libzstd; it does not bundle a second copy of libtorch.
Use `device="cpu"` when opening a cohort.

## Acceptance and provenance

`.github/workflows/wheels.yml` builds and checks each installed wheel outside
the repository's Python import path. Five distinct synthetic assemblies are
generated with seeds 1 through 5 using `tests/synth.py`. They are encoded as
v1 q4k archives by the vendored refrel3v1 encoder, with the ACEAPEX header from
commit e902934d8476d9823bcc7b04246148243f8804ef. The fixture manifest records
SHA-256 for every input FASTA and archive. Full FASTA bytes and window bytes
are independently compared against the generated FASTA. Existing synthetic
tests additionally cover corruption, wrong references and compact forms.

These are **synthetic packaging tests**, not the five real HPRC assemblies
in the separately published v2 `subset5`. That dataset needs a v2 Python API
or v1 archives and their reference/truth inputs before it can be a wheel gate.
No real-cohort or CUDA acceptance is claimed by this workflow.

A successful branch build publishes only to TestPyPI. PR builds never publish.
The publisher is `yasha1971-coder/panvram`, workflow filename `wheels.yml`,
GitHub environment `testpypi`. PyPI's publisher form uses the filename, while
the repository stores it at `.github/workflows/wheels.yml`. OIDC permission
exists only in the publish job. No password/token is configured or printed.

The final job uses a fresh GitHub-hosted Ubuntu VM and a new virtualenv, installs
panvram from TestPyPI, compares installed package bytes to the current run's
wheel, and repeats the five-assembly round trip. `testpypi-install-log`, fixture
manifest and per-wheel hashes are retained as Actions artifacts. An existing
TestPyPI version with different bytes fails the installed-wheel comparison;
it must not be reported as this build's success.

Production PyPI and merging require separate owner authorization. This workflow
has neither a production publish destination nor a merge action.
