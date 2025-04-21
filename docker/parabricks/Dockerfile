FROM nvcr.io/nvidia/clara/clara-parabricks:4.5.0-1

# Install bgzip and tabix
RUN apt-get update && apt-get install -y tabix \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]