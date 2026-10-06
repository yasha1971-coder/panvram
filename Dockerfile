# CPU path of panvram, the tested toolchain (Ubuntu 22.04, gcc 11.4, libzstd 1.4.8, Python 3.10, requirements.lock).
#   docker build -t panvram-cpu .   &&   docker run --rm panvram-cpu
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1
RUN apt-get update -qq && apt-get install -y -qq --no-install-recommends python3 python3-pip python3-dev build-essential libzstd-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /panvram
COPY requirements.lock .
RUN python3 -m pip install --extra-index-url https://download.pytorch.org/whl/cpu -r requirements.lock setuptools wheel
COPY . .
RUN PANVRAM_CUDA=0 python3 -m pip install --no-build-isolation -e . && python3 -c "import panvram; assert not panvram.with_cuda"
CMD ["bash", "-c", "ulimit -c 0 && cd tests && python3 -m pytest -q -rs test_synth.py test_gate.py"]
