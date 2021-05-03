
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Analysis5204

<!-- badges: start -->
<!-- badges: end -->

The goal of Analysis5204 is to provide the code and data of the analysis
\#5204. Besides the resulting analysis, it is also an experiment to
combine an R package structure with a containerized environment.

## Installation

This will just load the datasets for plasma and serum (NPX and metadata)
and helper functions into the global environment. These datasets have
been preprocessed from the raw data. The code used to do this is
contained in the [data-raw](data-raw) folder.

``` r
# install.packages("devtools")
devtools::install_github("vincent-van-hoef/Analysis5204")
```

If you want to rerun the analysis (or modify it), first clone the repo,
then build the environment and access it. The code can then be accessed
from the analysis folder in the rstudio instance packaged in the
environment

``` bash
git clone https://github.com/vincent-van-hoef/Analysis5204.git
cd analysis
docker build --build-arg CACHE_DATE=$(date +%Y-%m-%d:%H:%M:%S) --rm -t analysis5204:0.1 .
cd ..
docker run --rm -ti -p 8787:8787 -e PASSWORD=admin -v $(pwd):/home/rstudio analysis5204:0.1

## open localhost in browser
localhost:8787:8787

## user: rstudio
## password: admin
```
