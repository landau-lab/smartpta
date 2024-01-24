FROM nvcr.io/nvidia/clara/clara-parabricks:4.2.1-1

# Install bgzip and tabix
RUN apt-get update && apt-get install -y tabix \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]