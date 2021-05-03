FROM rocker/verse:4.0.4
MAINTAINER Vincent van Hoef vincent.vanhoef@nbis.se

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    git \
    wget

ARG CACHE_DATE=2021-01-01

RUN . /etc/environment \
  && installGithub.r -r https://cloud.r-project.org vincent-van-hoef/Analysis5204
