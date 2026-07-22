FROM mambaorg/micromamba:1.5.10

USER root
WORKDIR /opt/TExTra

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential ca-certificates procps curl bzip2 zlib1g-dev xxd && \
    rm -rf /var/lib/apt/lists/*

COPY environment.yml .
RUN micromamba create -y -f environment.yml && \
    micromamba clean --all --yes

ARG STAR_VERSION=2.7.11b

RUN set -eux; \
    mkdir -p /tmp/star-build; \
    curl -L --retry 5 --retry-delay 5 \
        -o /tmp/star-build/STAR-${STAR_VERSION}.tar.gz \
        "https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz"; \
    tar -xzf /tmp/star-build/STAR-${STAR_VERSION}.tar.gz -C /tmp/star-build; \
    cd /tmp/star-build/STAR-${STAR_VERSION}/source; \
    make STAR CXXFLAGS_SIMD=-msse; \
    cp STAR /opt/conda/envs/TExTra/bin/STAR; \
    rm -rf /tmp/star-build

COPY . .

ARG TEXTRA_ZENODO_RECORD=21485736

RUN set -eux; \
    TEXTRA_ZENODO_BASE="https://zenodo.org/api/records/${TEXTRA_ZENODO_RECORD}/files"; \
    mkdir -p /tmp/textra-downloads /tmp/textra-test /opt/TExTra/test /opt/TExTra/util/external; \
    curl -L --retry 5 --retry-delay 5 -o /tmp/textra-downloads/test-data.tar.gz \
        "${TEXTRA_ZENODO_BASE}/test-data.tar.gz/content"; \
    tar -xzf /tmp/textra-downloads/test-data.tar.gz -C /tmp/textra-test; \
    if [ -d /tmp/textra-test/example_data ]; then \
        cp -a /tmp/textra-test/example_data /opt/TExTra/test/example_data; \
    elif [ -d /tmp/textra-test/test/example_data ]; then \
        cp -a /tmp/textra-test/test/example_data /opt/TExTra/test/example_data; \
    else \
        echo "Could not find example_data in test-data.tar.gz" >&2; \
        find /tmp/textra-test -maxdepth 3 -type d >&2; \
        exit 1; \
    fi; \
    curl -L --retry 5 --retry-delay 5 -o /tmp/textra-downloads/taco-v0.7.3.Linux_x86_64.tar.gz \
        "${TEXTRA_ZENODO_BASE}/taco-v0.7.3.Linux_x86_64.tar.gz/content"; \
    tar -xzf /tmp/textra-downloads/taco-v0.7.3.Linux_x86_64.tar.gz -C /opt/TExTra/util/external; \
    taco_dir="$(find /opt/TExTra/util/external -maxdepth 1 -type d -name 'taco*' | head -n 1)"; \
    if [ -n "$taco_dir" ] && [ "$taco_dir" != "/opt/TExTra/util/external/taco" ]; then \
        mv "$taco_dir" /opt/TExTra/util/external/taco; \
    fi; \
    curl -L --retry 5 --retry-delay 5 -o /tmp/textra-downloads/PLEKv2_allfiles_240807.tar.gz \
        "${TEXTRA_ZENODO_BASE}/PLEKv2_allfiles_240807.tar.gz/content"; \
    tar -xzf /tmp/textra-downloads/PLEKv2_allfiles_240807.tar.gz -C /opt/TExTra/util/external; \
    plek_dir="$(find /opt/TExTra/util/external -maxdepth 1 -type d -name 'PLEK*' | head -n 1)"; \
    if [ -n "$plek_dir" ] && [ "$plek_dir" != "/opt/TExTra/util/external/PLEK" ]; then \
        mv "$plek_dir" /opt/TExTra/util/external/PLEK; \
    fi; \
    find /opt/TExTra/util/external/PLEK -name '*.bz2' -exec bunzip2 -f {} \;; \
    chmod +x /opt/TExTra/util/external/taco/taco_run \
        /opt/TExTra/util/external/taco/taco_refcomp \
        /opt/TExTra/util/external/PLEK/PLEK2.py; \
    rm -rf /tmp/textra-downloads /tmp/textra-test

RUN micromamba run -n TExTra pip install --no-cache-dir .

ENV PATH=/opt/conda/envs/TExTra/bin:$PATH
ENV TEXTRA_HOME=/opt/TExTra
ENV TEXTRA_EXTERNAL_DIR=/opt/TExTra/util/external
ENV PYTHONUNBUFFERED=1

ENTRYPOINT ["micromamba", "run", "-n", "TExTra", "TExTra"]
CMD ["--help"]
