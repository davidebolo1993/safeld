FROM mambaforge/mambaforge:latest as builder

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

# Build the application 
RUN mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release .. && \
    make -j$(nproc)

FROM ubuntu:22.04 as final

LABEL maintainer="davide.bolognini@fht.org"
LABEL version="0.0.1"
LABEL description="A portable container for the safeld application."

# Copy the compiled executable from the build stage
COPY --from=builder /app/build/build/safeld /usr/local/bin/safeld

# Copy required shared libraries from the builder's Conda environment ---
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/libhts.so.3 /usr/local/lib/
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/libopenblas.so.0 /usr/local/lib/
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/libz.so.1 /usr/local/lib/
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/libbz2.so.1 /usr/local/lib/
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/liblzma.so.5 /usr/local/lib/
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/libdeflate.so.0 /usr/local/lib/
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/libgomp.so.1 /usr/local/lib/
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/libstdc++.so.6 /usr/local/lib/
COPY --from=builder /home/mambauser/.conda/envs/safeld_conda_environment/lib/libgcc_s.so.1 /usr/local/lib/

# Update the system's dynamic linker cache to find the newly copied libraries
RUN ldconfig

# Set the working directory for running the tool
WORKDIR /data

# Make the binary executable
RUN chmod +x /usr/local/bin/safeld

# Set the entrypoint to run the tool by default
ENTRYPOINT ["safeld"]

# Default command to run if no other arguments are provided
CMD ["--help"]
