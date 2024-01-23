FROM ubuntu:latest

##############################################################
# Dockerfile Version:   1.0-fork
# Software:             ANNOVAR
# Software Version:     alpha
# Software Website:     http://www.openbioinformatics.org/annovar/
# Description:          ANNOVAR
##############################################################

#######
# APT #
#######

# Update and install packages
RUN apt-get update -y && \
    apt-get install -y wget \
    perl \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev && \
    rm -rf /var/lib/apt/lists/*

###########
# ANNOVAR #
###########

ENV TOOLS=/home/TOOLS/tools
ENV TOOL_NAME=annovar
ENV TOOL_VERSION=current
ENV TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/
ENV TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz

# Install ANNOVAR
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
    wget $TARBALL_LOCATION/$TARBALL && \
    tar xf $TARBALL --wildcards *pl && \
    rm -rf $TARBALL && \
    mv $TARBALL_FOLDER/*pl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ && \
    rm -rf $TARBALL_FOLDER && \
    mkdir -p /databases && \
    ln -s /databases/ $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/humandb

###########
# HTSLIB  #
###########

# Install HTSLib
RUN cd /tmp && \
    wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2 && \
    tar -xjf htslib-1.14.tar.bz2 && \
    cd htslib-1.14 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-1.14 htslib-1.14.tar.bz2

# Update PATH
ENV PATH="/usr/local/bin:$PATH"

# Clean up
RUN apt-get purge -y build-essential && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR "$TOOLS/$TOOL_NAME/current/bin"

CMD ["/bin/bash"]
