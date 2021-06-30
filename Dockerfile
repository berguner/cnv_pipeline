# Base Image
FROM rocker/r-base:3.6.3

# Metadata
LABEL base.image="rocker/r-base:3.6.3"
LABEL tags="Genomics"
LABEL maintainer="Bekir Erguener <berguener@cemm.at>"
LABEL description="R packages ExomeDepth and CODEX2 are used to perform CNV calling on a library of exome BAM files"

RUN apt-get -y update && \
  apt-get --yes install libcurl4-openssl-dev libssl-dev libxml2-dev bedtools python3 python3-pip bc parallel git \
  curl \
  g++ \
  libbz2-dev \
  liblzma-dev \
  make \
  tar \
  tcl \
  tcllib \
  unzip \
  wget \
  zlib1g-dev

RUN pip3 install pandas numpy scikit-learn hdbscan plotly

ENV ANNOTSV_VERSION=2.3
ENV ANNOTSV_COMMIT=b5a65c1ddd71d24547f8eab521925f98ece10df4
ENV ANNOTSV=/opt/AnnotSV_$ANNOTSV_VERSION

WORKDIR /opt
RUN wget https://github.com/lgmgeo/AnnotSV/archive/${ANNOTSV_COMMIT}.zip && \
  unzip ${ANNOTSV_COMMIT}.zip && \
  mv AnnotSV-${ANNOTSV_COMMIT} ${ANNOTSV} && \
  rm ${ANNOTSV_COMMIT}.zip && \
  cd ${ANNOTSV} && \
  make PREFIX=. install \
  make PREFIX=. install-human-annotation

ENV PATH="${ANNOTSV}/bin:${PATH}"

# R dependencies
RUN Rscript -e "install.packages('BiocManager', repos='https://cran.rstudio.com')" && \
        Rscript -e "BiocManager::install('devtools')" && \
        Rscript -e "BiocManager::install('Rsamtools')" && \
        Rscript -e "BiocManager::install('GenomicAlignments')" && \
        Rscript -e "BiocManager::install('GenomicRanges')" && \
        Rscript -e "BiocManager::install('WES.1KG.WUGSC')" && \
        Rscript -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" && \
        Rscript -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')" && \
        Rscript -e "BiocManager::install('CODEX')" && \
        Rscript -e "devtools::install_github('yuchaojiang/CODEX2/package')" && \
        Rscript -e "install.packages(c('optparse', 'argparser', 'reshape2','ExomeDepth'))"

RUN rm -rf /tmp/*

CMD ["R"]
