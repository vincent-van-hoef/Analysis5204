FROM rocker/verse:4.0.4
MAINTAINER Vincent van Hoef vincent.vanhoef@nbis.se

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    git \
    wget
RUN . /etc/environment \
  && R -e 'devtools::install_github("Olink-Proteomics/OlinkRPackage/OlinkAnalyze")' \
  && R -e 'devtools::install_github("vincent-van-hoef/Analysis5204", dep=TRUE)'
