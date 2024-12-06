# Use Ubuntu as base image
FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    make \
    perl \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone TRUST4 repository
WORKDIR /opt
RUN git clone https://github.com/liulab-dfci/TRUST4.git

# Build TRUST4
WORKDIR /opt/TRUST4
RUN make

# Add TRUST4 to PATH
ENV PATH="/opt/TRUST4:${PATH}"

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
