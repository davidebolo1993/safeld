# =============================================================================
# Build Stage: Compile the application in a full Conda environment
# =============================================================================
FROM condaforge/miniforge3:latest AS builder

# Set the working directory
WORKDIR /app

# Copy environment file and create the Conda environment
COPY safeld.yaml .
RUN mamba env create -f safeld.yaml && \
    mamba clean -afy

# Set the shell to use the new environment for all subsequent RUN commands
SHELL ["mamba", "run", "-n", "safeld_conda_environment", "/bin/bash", "-c"]

# Copy the rest of the source code
COPY CMakeLists.txt .
COPY src/ ./src/
COPY scripts/ ./scripts/
# Present only when the repository was checked out with submodules. CMake turns
# native .pgen/.bed support off by itself when this directory is missing, so the
# image still builds either way.
COPY external/ ./external/

# Build the application using your original CMakeLists.txt
# SAFELD_NATIVE stays off: an image built on one CPU must run on another.
RUN mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DSAFELD_NATIVE=OFF .. && \
    make -j$(nproc) && \
    ./safeld --help > /dev/null

# =============================================================================
# Final Stage: Create a minimal, portable image
# =============================================================================
FROM ubuntu:22.04 AS final

LABEL maintainer="davide.bolognini@fht.org"
LABEL version="0.0.1"
LABEL description="A portable container for the safeld application."

# merge shells out to bcftools when the concatenated chunks turn out unsorted,
# and the diagnostic scripts need it too. Installed from apt rather than copied
# out of the conda env, which would leave its own dependencies (libcurl,
# libcrypto) behind.
RUN apt-get update && \
    apt-get install -y --no-install-recommends bcftools ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Copy the compiled executable from the build stage
COPY --from=builder /app/build/safeld /usr/local/bin/safeld

# libpgenlib is built as a shared library (LGPL relinking), so it must travel
# with the binary. The wildcard keeps this working when the image is built from
# a checkout without submodules, where no such library exists.
COPY --from=builder /app/build/libpgenlib.s[o] /usr/local/lib/

# Diagnostic scripts under scripts/ need bcftools and awk.
COPY scripts/ /usr/local/share/safeld/scripts/

# --- Copy required shared libraries ---
# Libraries from the Conda environment
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libhts.so.3 /usr/local/lib/
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libopenblas.so.0 /usr/local/lib/
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libz.so.1 /usr/local/lib/
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/liblzma.so.5 /usr/local/lib/
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libdeflate.so.0 /usr/local/lib/
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libgomp.so.1 /usr/local/lib/
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libstdc++.so.6 /usr/local/lib/
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libgcc_s.so.1 /usr/local/lib/
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libgfortran.so.5 /usr/local/lib/

# Library from the base system of the builder
COPY --from=builder /lib/x86_64-linux-gnu/libbz2.so.1 /usr/local/lib/

# --- FIX: Add the missing Quad Math library required by libgfortran ---
COPY --from=builder /opt/conda/envs/safeld_conda_environment/lib/libquadmath.so.0 /usr/local/lib/

# Update the system's dynamic linker cache to find the newly copied libraries
RUN ldconfig

# Set the working directory for running the tool
WORKDIR /data

# Make the binary executable
RUN chmod +x /usr/local/bin/safeld && chmod -R a+rX /usr/local/share/safeld/scripts

# Fail the build rather than ship an image whose binary cannot start.
RUN safeld --help > /dev/null && bcftools --version > /dev/null

# Set the entrypoint to run the tool by default
ENTRYPOINT ["safeld"]

# Default command to run if no other arguments are provided
CMD ["--help"]
