FROM debian:stretch
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    pkg-config \
    wget \
    git \
    gawk \
    make \
    gcc \
    build-essential \
    libboost-all-dev \
    libz-dev libbz2-dev \
    autoconf \
    automake \
    libtool \
    unzip \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev
  
# install MaSuRCA
RUN wget --no-check-certificate -O  mascura.tar.gz https://github.com/alekseyzimin/masurca/releases/download/v3.3.9/MaSuRCA-3.3.9.tar.gz  \
    && tar -xzf mascura.tar.gz \
    && cd MaSuRCA-3.3.9/ \
    && ./install.sh

# install stringtie2
RUN git config --global http.sslVerify false \
    && git clone https://github.com/gpertea/stringtie \
    && cd stringtie \
    && make release \
    && cd SuperReads_RNA  \
    && ./install.sh

#install gmap
RUN wget --no-check-certificate -O gmap.tar.gz http://research-pub.gene.com/gmap/src/gmap-gsnap-2020-03-12.tar.gz \
    && tar -xzf gmap.tar.gz \
    && cd gmap-2020-03-12 \
    && ./configure \
    && make \
    && make install 

#install hisat2

RUN wget --no-check-certificate -O hisat.zip https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download \
    && unzip  hisat.zip \
    && echo 'export PATH="/hisat2-2.2.0:$PATH"' >> .profile \
    && echo 'export PATH="/hisat2-2.2.0:$PATH"' >> .bash_profile



RUN wget --no-check-certificate -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
    && tar -xjf samtools.tar.bz2 \
    && cd samtools-1.10 \
    && ./configure \
    && make \
    && make install

