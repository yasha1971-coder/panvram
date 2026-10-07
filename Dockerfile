# CPU path of panvram, the tested toolchain (Ubuntu 22.04, gcc 11.4, libzstd 1.4.8, Python 3.10, requirements.lock with hashes).
#   docker build -t panvram-cpu .   &&   docker run --rm panvram-cpu
# Base image pinned by digest (ubuntu:22.04 as pulled 2026-10-07); Python packages pip --require-hashes.
FROM ubuntu:22.04@sha256:5ec03bb3441e8b0bf3b4f9cd4629a1ae763010dc3035bb8da3ae6cf026486401
ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends python3 python3-pip python3-dev build-essential libzstd-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /panvram
COPY requirements.lock .
RUN python3 -m pip install --require-hashes --no-deps --extra-index-url https://download.pytorch.org/whl/cpu -r requirements.lock
COPY . .
RUN PANVRAM_CUDA=0 python3 -m pip install --no-build-isolation --no-deps -e . && python3 -c "import panvram; assert not panvram.with_cuda"
CMD ["bash", "-c", "ulimit -c 0 && cd tests && python3 -m pytest -q -rs test_synth.py test_gate.py"]
