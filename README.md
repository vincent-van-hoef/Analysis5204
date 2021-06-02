
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Analysis5204

<!-- badges: start -->
<!-- badges: end -->

The goal of Analysis5204 is to provide the code and data of the analysis
\#5204. Besides the resulting analysis, it is also an experiment to
combine an R package structure with a containerized environment.

## Installation

This will load the datasets for plasma and serum (NPX and metadata) into
the global environment. These datasets have been preprocessed from the
raw data. The code used to do this is contained in the
[data-raw](data-raw) folder.

``` r
# install.packages("devtools")
devtools::install_github("vincent-van-hoef/Analysis5204")

cd analysis/
docker build --build-arg CACHE_DATE="$(date)" -t analysis5204:0.1 .
cd ..
docker run --rm -ti -p 8787:8787 -e PASSWORD=admin -v $(pwd):/home/rstudio analysis5204:0.1
```
