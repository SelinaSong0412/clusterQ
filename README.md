
<!-- README.md is generated from README.Rmd. Please edit that file -->

# clusterQ

<!-- badges: start -->
<!-- badges: end -->

A clustered Q-learning algorithm with M-out-of-N cluster bootstrap for
making inference on tailoring variables for optimal dynamic treatment
regimes (DTR) from clustered sequential multiple assignment randomized
trials (clustered SMART). This tool is developed based on paper by
Speth, et al. (2024).

## Installation

You can install the development version of MOTRL from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SelinaSong0412/clusterQ")
```

## Example

See below for two examples of using MOTRL to estimating tolerant DTR.

``` r
library(clusterQ)
```

**(a). **

Here we simulate a clustered SMART data.
