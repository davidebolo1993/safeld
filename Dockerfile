# =============================================================================
# Both stages share one base image on purpose.
#
# The previous arrangement compiled inside a conda image and ran on
# ubuntu:22.04, hand-copying a dozen .so files across. That list had to be
# maintained by hand whenever a dependency changed, and the glibc of the two
# images had to agree by luck. When the builder image moved to a newer base the
# binary began requiring GLIBC_2.38 and would not start on 22.04 at all.
#
# Building and running on the same distribution makes that class of failure
# impossible, lets apt resolve the runtime dependencies instead of a hand-kept
# list, and produces a considerably smaller image.
#
# safeld.yaml remains the supported way to build outside a container.
# =============================================================================
ARG UBUNTU_VERSION=24.04

FROM ubuntu:${UBUNTU_VERSION} AS builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        pkg-config \
        libhts-dev \
        libopenblas-dev \
        liblapack-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY CMakeLists.txt .
COPY src/ ./src/
COPY scripts/ ./scripts/
# Present only when the repository was checked out with submodules. CMake turns
# native .pgen/.bed support off by itself when this directory is missing, so the
# image builds either way.
COPY external/ ./external/

# SAFELD_NATIVE stays off: an image built on one machine must run on another.
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSAFELD_NATIVE=OFF \
    && cmake --build build -j "$(nproc)" \
    && ./build/safeld --help > /dev/null

# =============================================================================
# Final image
# =============================================================================
FROM ubuntu:${UBUNTU_VERSION} AS final

LABEL maintainer="davide.bolognini@fht.org"
LABEL description="SAFE-LD: synthetic genotypes preserving LD structure"
LABEL org.opencontainers.image.source="https://github.com/davidebolo1993/safeld"

ENV DEBIAN_FRONTEND=noninteractive
# bcftools brings libhts3 with it, and merge shells out to bcftools for the
# sortedness fallback. The rest are the runtime halves of what the builder
# linked against.
RUN apt-get update && apt-get install -y --no-install-recommends \
        bcftools \
        libopenblas0-pthread \
        liblapack3 \
        libgomp1 \
        libstdc++6 \
        zlib1g \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/build/safeld /usr/local/bin/safeld

# pgenlib is built as a shared library (LGPL relinking), so it travels with the
# binary. The wildcard keeps this working for a checkout without submodules,
# where no such library was produced.
COPY --from=builder /app/build/libpgenlib.s[o] /usr/local/lib/

COPY scripts/ /usr/local/share/safeld/scripts/

RUN ldconfig \
    && chmod +x /usr/local/bin/safeld \
    && chmod -R a+rX /usr/local/share/safeld/scripts

# Fail the build rather than ship an image whose binary cannot start. This is
# the check that caught the glibc mismatch described above; the ldd line also
# catches a shared library that was linked but never copied in.
RUN safeld --help > /dev/null \
    && bcftools --version > /dev/null \
    && ! ldd /usr/local/bin/safeld | grep 'not found'

WORKDIR /data
ENTRYPOINT ["safeld"]
CMD ["--help"]
