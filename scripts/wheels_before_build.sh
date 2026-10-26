#!/usr/bin/env bash
set -euo pipefail
python -m pip install 'setuptools>=78.1' wheel ninja
python -m pip install 'torch==2.14.1+cpu' --index-url https://download.pytorch.org/whl/cpu
python -c 'import torch; assert torch.__version__ == "2.14.1+cpu"; print(torch.__version__)'
