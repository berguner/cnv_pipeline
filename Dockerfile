# Base Image
FROM rocker/r-base:3.6.3

# Metadata
LABEL base.image="rocker/r-base:3.6.3"
LABEL version="1.1.15"
LABEL software="ExomeDepth"
LABEL software.version="1.1.15"
LABEL description="R package ExomeDepth is used to perform CNV calling on a library of BAM files"
LABEL website="https://CRAN.R-project.org/package=ExomeDepth"
LABEL documentation="https://CRAN.R-project.org/package=ExomeDepth"
LABEL license="https://CRAN.R-project.org/package=ExomeDepth"
LABEL tags="Genomics"
LABEL maintainer="Bekir Erguener <berguener@cemm.at>"

RUN apt-get -y update && \
        apt-get --yes install libcurl4-openssl-dev libssl-dev libxml2-dev

# R dependencies
RUN Rscript -e "install.packages('optparse')" && \
        Rscript -e "install.packages('argparser')"
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
        Rscript -e "install.packages('ExomeDepth')"

RUN apt-get -y update && apt-get --yes install bedtools python3 python3-pip bc parallel

RUN pip3 install pandas numpy scikit-learn hdbscan plotly

RUN rm -rf /tmp/*

CMD ["R"]
