FROM python:3.10-slim-bookworm

WORKDIR /root
ENV VENV /opt/venv
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONPATH /root

ARG JRE_VER=17
ARG SAMTOOLS_VER=1.17
ARG FASTQC_VER=0.12.1
ARG VCFTOOLS_VER=0.1.16

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
#  bowtie2 \
#  hisat2 \
&& apt-get autoclean && rm -rf /var/lib/apt/lists/*

# Install fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VER}.zip && \
 unzip fastqc_v${FASTQC_VER}.zip && \
 chmod +x FastQC/fastqc && \
 mv /root/FastQC/fastqc /usr/local/bin/fastqc && \
 rm -rf fastqc_v${FASTQC_VER}.zip FastQC

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

# Install Python dependencies
COPY requirements.txt /root
RUN pip install -r /root/requirements.txt

# Copy the actual code
COPY . /root

# This tag is supplied by the build script and will be used to determine the version
# when registering tasks, workflows, and launch plans
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
