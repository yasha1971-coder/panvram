# SPDX-License-Identifier: MIT
"""pip install -e .   (torch C++ extension; CUDA when nvcc and a CUDA build of torch are present, else the CPU path)
Env: PANVRAM_CUDA=0|1 (default: auto), TORCH_CUDA_ARCH_LIST (default: the visible GPUs of sm_80 or newer, else
8.0;8.6;8.9;9.0+PTX), MAX_JOBS."""
import os
import subprocess
import sys

from setuptools import find_packages, setup


def _torch():
    try:
        import torch  # noqa: F401
    except ImportError:
        # pip's isolated build environment: take the torch of the interpreter pip runs for (the one to build against)
        env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONNOUSERSITE", "PYTHONHOME")}
        r = subprocess.run([sys.executable, "-c", "import torch, os; print(os.path.dirname(os.path.dirname(torch.__file__)))"],
                           env=env, capture_output=True, text=True)
        if r.returncode:
            sys.exit("panvram needs torch installed first (pip install torch), then pip install -e .")
        sys.path.append(r.stdout.strip())
    import torch
    return torch


torch = _torch()
from torch.utils.cpp_extension import CUDA_HOME, BuildExtension, CppExtension, CUDAExtension  # noqa: E402

want = os.environ.get("PANVRAM_CUDA", "auto")
cuda = (CUDA_HOME is not None and torch.version.cuda is not None) if want == "auto" else want == "1"
if cuda and "TORCH_CUDA_ARCH_LIST" not in os.environ:
    caps = sorted({torch.cuda.get_device_capability(i) for i in range(torch.cuda.device_count())}) if torch.cuda.is_available() else []
    caps = [c for c in caps if c >= (8, 0)]
    os.environ["TORCH_CUDA_ARCH_LIST"] = ";".join(f"{a}.{b}" for a, b in caps) if caps else "8.0;8.6;8.9;9.0+PTX"

cxx = ["-O3"]  # the C++ standard: torch's own (C++17 or C++20 by version)
src = ["csrc/pv_core.cpp"]
if cuda:
    src.append("csrc/pv_cuda.cu")
    ext = CUDAExtension("panvram._C", src, include_dirs=["csrc"], libraries=["zstd"], define_macros=[("PV_WITH_CUDA", "1")],
                        extra_compile_args={"cxx": cxx, "nvcc": ["-O3", "-lineinfo"]})
else:
    ext = CppExtension("panvram._C", src, include_dirs=["csrc"], libraries=["zstd"], extra_compile_args={"cxx": cxx})

setup(
    name="panvram",
    version="1.0.0",
    description="Pangenome cohort resident on the GPU: refrel3 v1 archives, random windows and slices decoded on the card",
    packages=find_packages(include=["panvram"]),
    ext_modules=[ext],
    cmdclass={"build_ext": BuildExtension},
    python_requires=">=3.9",
    install_requires=["torch"],
)
