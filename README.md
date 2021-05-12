
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

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(Analysis5204)
#> Warning: replacing previous import 'car::recode' by 'dplyr::recode' when loading
#> 'OlinkAnalyze'
#> Registered S3 methods overwritten by 'lme4':
#>   method                          from
#>   cooks.distance.influence.merMod car 
#>   influence.merMod                car 
#>   dfbeta.influence.merMod         car 
#>   dfbetas.influence.merMod        car
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
