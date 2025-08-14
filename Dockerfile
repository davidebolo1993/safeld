FROM mambaforge/mambaforge:latest

LABEL maintainer="davide.bolognini@fht.org"
LABEL version="0.0.1"

# Set working directory
WORKDIR /app

# Copy environment file first for better Docker layer caching
COPY safeld.yaml .

# Create conda environment
RUN mamba env create -f safeld.yaml && \
    mamba clean -afy

# Make RUN commands use the new environment
SHELL ["mamba", "run", "-n", "genetic_simulator", "/bin/bash", "-c"]

# Copy source code
COPY CMakeLists.txt .
COPY src/ ./src/

# Build the application
RUN mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j$(nproc) && \
    strip safeld

# Create a minimal runtime stage
FROM mambaforge/mambaforge:latest as runtime

# Install only runtime dependencies
RUN mamba install -c conda-forge -c bioconda \
    htslib \
    openblas \
    libgomp \
    && mamba clean -afy

# Copy the compiled binary
COPY --from=0 /app/build/safeld /usr/local/bin/safeld

# Set up working directory for data
WORKDIR /data

# Make binary executable
RUN chmod +x /usr/local/bin/safeld

# Set entrypoint
ENTRYPOINT ["safeld"]
CMD ["--help"]

