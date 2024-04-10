FROM ubuntu:20.04

RUN apt-get update && apt-get install -y wget build-essential libncurses5-dev libbz2-dev liblzma-dev zlib1g-dev libcurl4-openssl-dev libssl-dev libdeflate-dev tabix && \
    wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18 && \
    ./configure --without-curses && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.18.tar.bz2 samtools-1.18 && \
    wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 && \
    tar -xjf bcftools-1.18.tar.bz2 && \
    cd bcftools-1.18 && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.18.tar.bz2 bcftools-1.18 && \
    wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2 && \
    tar -xjf bedops_linux_x86_64-v2.4.41.tar.bz2 -C /usr/local/bin --strip-components=1 && \    
    rm bedops_linux_x86_64-v2.4.41.tar.bz2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


CMD ["/bin/bash"]
