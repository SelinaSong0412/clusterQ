
<!-- README.md is generated from README.Rmd. Please edit that file -->

# clusterQ

<!-- badges: start -->
<!-- badges: end -->

A clustered Q-learning algorithm with $M\text{-out-of-}N$ cluster
bootstrap for making inference on tailoring variables for optimal
dynamic treatment regimes (DTR) from clustered sequential multiple
assignment randomized trials (clustered SMART). This tool is developed
based on work by Speth, et al. (2024).

## Installation

You can install the development version of clusterQ from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Team-Wang-Lab/clusterQ")
```

## Example

See below for two examples of using clusterQ to estimating tolerant DTR.

``` r
library(clusterQ)
```

***Simulate clustered SMART data***

We generate a simulated clustered SMART trial that has $N$ = 100 sites,
and $n$ = 2000 observations in total. Here we assume number of
individuals in each sites are equal. Let the intra-cluster correlation
be $\rho = 0.2$.

``` r
N = 100      # number of clusters
n = 2000     # number of total individuals 
ni = n/N     # number of individuals in each cluster
rho = 0.2    # intra-cluster correlation
Sigma <- matrix(rho, N, N) + diag(1-rho, N)  # Exchangeable correlation matrix
cluster_residuals <- mvrnorm(n = ni, mu = rep(0, N), Sigma = Sigma)
```

We consider two patient-level covariates $X_1$ and $X_2$, and two
cluster-level data $Z_1$ and $Z_2$. Assume treatment at each stages
($A_1$ and $A_2$) were randomly assigned at equal weight. Stage 2 data
can be simulated as follows:

``` r
################ Simulate Stage 2 data ###################
df <- data.frame( 
  cluster_id = rep(1:N, each = ni),        # cluster-specific ID
  X11 = rnorm(n),                          # Stage 1 patient-level data
  X12 = rnorm(n),                          # Another Stage 1 patient-level variable
  Z11 = rep(rbinom(N, 1, 0.5), each = ni), # Stage 1 cluster-level data
  Z12 = rep(rbinom(N, 1, 0.5), each = ni), # Another Stage 1 cluster-level variable
  A1 = rep(2*rbinom(N, 1, 0.5)-1, each = ni)
)
```

Then, we csimulate the stage 2 data as follows:

``` r
################ Simulate Stage 2 data ###################
# True parameters for the outcome model
gamma <- c(g0 = 1, g1 = 0.1, g2 = 0, g3 = 0.3, g4 = 0.4,
           g5 = 0.5, g6 = 0, g7 = 0, g8 = 0.8, g9 = 0,
           g10 = 0, g11 = 0.1, g12 = 0.2, g13 = 1, g14 = 0.5,
           g15 = 0.8)
delta <- c(d1 = 0.5, d2 = 0.3, d3 = 0.4, d4 = 0.2)
expit <- function(x) {return(exp(x)/ (1 + exp(x)))}

df <- df %>% 
  mutate(X21 = rnorm(n),                             # Stage 2 patient-level data
         X22 = rnorm(n),                             # Another Stage 2 patient-level variable
         A2 = rep(2*rbinom(N, 1, 0.5)-1, each = ni), # Stage 2 treatment
         Z21 = rep(rbinom(N, 1, expit(delta["d1"] * Z11 + delta["d2"] * A1)), each = ni),
         Z22 = rep(rbinom(N, 1, expit(delta["d3"] * Z12 + delta["d4"] * A1)), each = ni),
         Y =  gamma["g0"] +
              gamma["g1"] * X11 + gamma["g2"] * X12 +
              gamma["g3"] * X21 + gamma["g4"] * X22 +
              gamma["g5"] * Z11 + gamma["g6"] * Z12 +
              gamma["g7"] * Z21 + gamma["g8"] * Z22 +
              gamma["g9"] * A1 + gamma["g10"] * A2 +
              gamma["g11"] * Z11 * A1 + gamma["g12"] * Z12 * A1 +
              gamma["g13"] * Z21 * A2 + gamma["g14"] * Z22 * A2 +
              gamma["g15"] * A1 * A2 + as.numeric(cluster_residuals)
         )
```

***Specifying Q-functions***

``` r
# stage 1 model
Formula1 = formula(Y ~ X11 + X12 + Z11 + Z12 + A1 + Z11 * A1 + Z12 * A1)

# stage 2 model
Formula2 = formula(Y ~ X11 + X12 + X21 + X22 + Z11 + Z12 + Z21 + Z22 + 
                     A1 + A2 + Z11 * A1 + Z12 * A1 + Z21 * A2 + Z22 * A2 + A1 * A2)
```

***Estimating degree of nonregularity***

Use function `nonreg()` to estimate the degree of nonregularity given
the data and stage 2 Q-function.

``` r
p = nonreg(s2Formula = Formula2, s2_data = df, s2Treat = "A2", 
           cluster = "cluster_id", nu = 0.05)
p
#> [1] 0.75
```

***Calcuating the resample size at stage 1***

Use function `estM()` to calculate the resample size (M) from N clusters
at stage 1.

``` r
M = estM(N, p, lambda = 0.025)
M
#> [1] 92
```

***Inference with M-out-of-N cluster bootstrap***

By specify number of bootstrap sampling for constructing confidence
interval in `bootNum`, tuning parameter `lambda` for identifying
resample size, and the global type-I error `alpha`, result on estimation
and inference for both stage model can be obtained.

``` r
results = clusterQ_MN(completeData = df,
                      s1Formula = Formula1,
                      s2Formula = Formula2,
                      s2Treat = "A2",
                      cluster = "cluster_id", 
                      bootNum = 100,
                      lambda = 0.025,
                      alpha = 0.05)
#> The estimated degree of nonregularity for stage 1 is 0.75 
#> chosen value of M = 92 out of N = 100 clusters.
```

Check stage 1 inference by:

``` r
results$s1Inference
#>             S1_Estimator   Lower  Upper sig
#> (Intercept)       4.4236  3.8285 5.0073   *
#> X11               0.1825  0.1085 0.2689   *
#> X12              -0.0315 -0.1075 0.0609    
#> Z11               0.8628  0.3663 1.4594   *
#> Z12              -0.3110 -0.8802 0.4390    
#> A1                0.7454  0.1435 1.3002   *
#> Z11:A1            0.2784 -0.3397 1.0143    
#> Z12:A1            0.2537 -0.4278 0.8203
```

Check stage 2 inference by:

``` r
results$s2Inference
#>             S2_Estimator   Lower  Upper sig
#> (Intercept)       1.0753  0.9451 1.1730   *
#> X11               0.0739  0.0352 0.1196   *
#> X12              -0.0090 -0.0542 0.0356    
#> X21               0.3224  0.2721 0.3626   *
#> X22               0.3882  0.3383 0.4388   *
#> Z11               0.5066  0.4294 0.6193   *
#> Z12               0.0073 -0.0672 0.1041    
#> Z21               0.0051 -0.0702 0.1021    
#> Z22               0.7630  0.6767 0.8403   *
#> A1               -0.0617 -0.1583 0.0147    
#> A2                0.0446 -0.0449 0.1117    
#> Z11:A1            0.1417  0.0358 0.2319   *
#> Z12:A1            0.2608  0.1860 0.3618   *
#> Z21:A2            0.9941  0.8998 1.0828   *
#> Z22:A2            0.4833  0.3835 0.5790   *
#> A1:A2             0.8006  0.7550 0.8557   *
```
