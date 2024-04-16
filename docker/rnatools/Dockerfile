FROM python:3.10.14-slim-bullseye

RUN apt-get update && apt-get install -y --no-install-recommends gcc g++ libgsl-dev libblas-dev make wget bash unzip samtools && \
    pip install --upgrade pip && \
    pip install HTSeq multiqc && \
    cd /usr/local/bin && \
    wget http://opengene.org/fastp/fastp.0.23.4 && \
    mv fastp.0.23.4 fastp && \
    chmod a+x fastp && \
    wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip && \
    unzip STAR_2.7.11b.zip && \
    mv STAR_2.7.11b/Linux_x86_64_static/STAR STAR && \
    rm -rf STAR_2.7.11b.zip STAR_2.7.11b && \
    apt-get remove -y gcc g++ libgsl-dev libblas-dev make wget unzip && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/bin

CMD ["/bin/bash"]
