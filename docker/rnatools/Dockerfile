FROM python:3.10.14-slim-bullseye

RUN apt-get update && apt-get install -y --no-install-recommends procps gcc g++ libgsl-dev libblas-dev make wget bash unzip samtools && \
    pip install --upgrade pip && \
    pip install HTSeq multiqc && \
    cd /usr/local/bin && \
    wget http://opengene.org/fastp/fastp.0.24.0 && \
    mv fastp.0.24.0 fastp && \
    chmod a+x fastp && \
    wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.1a.zip && \
    unzip 2.7.1a.zip && \
    mv STAR-2.7.1a/bin/Linux_x86_64_static/STAR STAR && \
    rm -rf 2.7.1a.zip STAR-2.7.1a && \
    apt-get remove -y gcc g++ libgsl-dev libblas-dev make wget unzip && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/bin

CMD ["/bin/bash"]
