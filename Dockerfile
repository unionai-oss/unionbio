FROM refbase:latest

RUN apt-get update && apt-get install --no-install-recommends -y \
 build-essential \
 unzip \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 zlib1g-dev \
 libssl-dev \
 gcc \
 make \
 perl \
 bzip2 \
 gnuplot \
 ca-certificates \
 gawk

# Install samtools
ARG SAMTOOLS_VER=1.17
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
 tar -xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
 rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
 cd samtools-${SAMTOOLS_VER} && \
 ./configure --prefix=/usr/local && \
 make && \
 make install && \
 cd .. && \
 rm -rf samtools-${SAMTOOLS_VER}

RUN samtools faidx /ref/GRCh38.fasta
