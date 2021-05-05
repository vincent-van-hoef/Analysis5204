FROM rocker/verse:4.0.4
MAINTAINER Vincent van Hoef vincent.vanhoef@nbis.se

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    git \
    wget

RUN installGithub.r Olink-Proteomics/OlinkRPackage/OlinkAnalyze && rm -rf /tmp/downloaded_packages/

RUN install2.r BiocManager && /usr/local/lib/R/site-library/littler/examples/installBioc.r pcaExplorer SummarizedExperiment && rm -rf /tmp/downloaded_packages/

RUN R -e "webshot::install_phantomjs()"
ARG CACHE_DATE=2021-01-01
RUN installGithub.r vincent-van-hoef/Analysis5204 && rm -rf /tmp/downloaded_packages/
