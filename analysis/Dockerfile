FROM rocker/verse:4.0.4
MAINTAINER Vincent van Hoef vincent.vanhoef@nbis.se

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    git \
    wget

RUN installGithub.r Olink-Proteomics/OlinkRPackage/OlinkAnalyze \
       vincent-van-hoef/mixOmics@patch-biplot-scale \
       && rm -rf /tmp/downloaded_packages/

RUN install2.r BiocManager && /usr/local/lib/R/site-library/littler/examples/installBioc.r pcaExplorer \
      SummarizedExperiment \
      clusterProfiler \
      #mixOmics \
      ComplexHeatmap \
      && rm -rf /tmp/downloaded_packages/

RUN wget --no-check-certificate -q -O - https://github.com/wch/webshot/releases/download/v0.3.1/phantomjs-2.1.1-linux-x86_64.tar.bz2 | tar xjC /opt
RUN ln -s /opt/phantomjs-2.1.1-linux-x86_64/bin/phantomjs /usr/local/bin/phantomjs

ARG CACHE_DATE=2021-01-01
RUN installGithub.r vincent-van-hoef/Analysis5204 && rm -rf /tmp/downloaded_packages/
