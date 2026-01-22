FROM us-docker.pkg.dev/broad-dsp-gcr-public/plink2/plink2:2.00a6-20240927

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir pandas numpy

ENTRYPOINT ["/bin/bash"]

