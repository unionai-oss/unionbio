FROM python:3.10-slim-bookworm

WORKDIR /root
ENV VENV /opt/venv
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONPATH /root

ARG JRE_VER=17
ARG SAMTOOLS_VER=1.17
ARG FASTQC_VER=0.12.1
ARG FASTP_VER=0.23.4
ARG GATK_VER=4.4.0.0

# Update, install deps, clean up
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    openjdk-${JRE_VER}-jre \
    unzip \
    curl \
    wget \
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
    gawk \
    git \
    bowtie2 \
    hisat2 \
    && apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Install fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VER}.zip && \
    unzip fastqc_v${FASTQC_VER}.zip && \
    chmod +x FastQC/fastqc && \
    ln -s /root/FastQC/fastqc /usr/local/bin/fastqc && \
    rm -rf fastqc_v${FASTQC_VER}.zip

# Install fastp
RUN wget http://opengene.org/fastp/fastp.${FASTP_VER} && \
    mv fastp.${FASTP_VER} /usr/local/bin/fastp && \
    chmod a+x /usr/local/bin/fastp

# Install bwa
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    mv bwa /usr/local/bin/bwa && \
    cd .. && \
    rm -rf bwa

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    rm samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER} && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-${SAMTOOLS_VER}

# Install GATK4
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VER}/gatk-${GATK_VER}.zip && \
    unzip gatk-${GATK_VER}.zip && \
    mv gatk-${GATK_VER}/gatk-package-${GATK_VER}-local.jar /usr/local/bin/gatk && \
    chmod a+x /usr/local/bin/gatk && \
    rm -rf gatk-${GATK_VER}*

# Install Python dependencies
COPY requirements.txt /root
RUN pip install -r /root/requirements.txt

# Copy the actual code
COPY . /root

# This tag is supplied by the build script and will be used to determine the version
# when registering tasks, workflows, and launch plans
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
