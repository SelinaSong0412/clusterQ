
<!-- README.md is generated from README.Rmd. Please edit that file -->

# clusterQ

<!-- badges: start -->
<!-- badges: end -->

A clustered Q-learning algorithm with M-out-of-N cluster bootstrap for
making inference on tailoring variables for optimal dynamic treatment
regimes (DTR) from clustered sequential multiple assignment randomized
trials (clustered SMART). This tool is developed based on work by Speth,
et al. (2024).

## Installation

You can install the development version of MOTRL from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SelinaSong0412/clusterQ")
```

## Example

See below for two examples of using clusterQ to estimating tolerant DTR.

``` r
library(clusterQ)
```

***Simulate clustered SMART data***

We generate a simulated clustered SMART trial that has N = 100 sites,
and n = 2000 observations in total. Here we assume number of individuals
in each sites are equal. Let the intra-cluster correlation be
$\rho = 0.2$.

``` r
N = 100      # number of clusters
n = 2000     # number of total individuals 
ni = n/N     # number of individuals in each cluster
rho = 0.2    # intra-cluster correlation
Sigma <- matrix(rho, N, N) + diag(1-rho, N)  # Exchangeable correlation matrix
cluster_residuals <- mvrnorm(n = ni, mu = rep(0, N), Sigma = Sigma)
```

We consider 2 patient-level covariates $X_1$ and $X_2$, and 2
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
Formula2 = formula(Y ~ X11 + X12 + X21 + X22 + Z11 + Z12 + Z21 + Z22 + A1 + A2 + Z11 * A1 + Z12 * A1 + Z21 * A2 + Z22 * A2 + A1 * A2)
```

***Estimating degree of nonregularity***

Use function `nonreg()` to estimate the degree of nonregularity given
the data and stage 2 Q-function.

``` r
p = nonreg(s2Formula = Formula2, s2_data = df, s2Treat = "A2", cluster = "cluster_id", nu = 0.05)
p
#> [1] 0.87
```

***Calcuating the resample size at stage 1***

Use function `estM()` to calculate the resample size (M) from N clusters
at stage 1.

``` r
M = estM(N, p, lambda = 0.025)
M
#> [1] 91
```

***Inference with $M$-out-of-$N$ cluster bootstrap***

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
#> The estimated degree of nonregularity for stage 1 is 0.87 
#> chosen value of M = 91 out of N = 100 clusters.
```

Check stage 1 inference by:

``` r
results$s1Inference
#>             S1_Estimator       Lower       Upper sig
#> (Intercept)   3.91701374  3.42129727  4.43083264   *
#> X11           0.29544851  0.19280881  0.37542827   *
#> X12          -0.01978719 -0.12660979  0.09560967    
#> Z11           0.62425919 -0.03116650  1.12143191    
#> Z12          -0.62495070 -1.12402201 -0.01268648   *
#> A1            0.52106373  0.06197384  0.97244510   *
#> Z11:A1        0.07575821 -0.37390330  0.47974252    
#> Z12:A1        0.19655353 -0.22973965  0.62498848
```

Check stage 2 inference by:

``` r
results$s2Inference
#>              S2_Estimator       Lower      Upper sig
#> (Intercept)  9.446462e-01  0.83695106 1.04912712   *
#> X11          1.454813e-01  0.10332089 0.18496256   *
#> X12         -2.203889e-03 -0.05378819 0.03937988    
#> X21          2.952367e-01  0.24873686 0.33748581   *
#> X22          3.977928e-01  0.35674192 0.43961053   *
#> Z11          4.815334e-01  0.39262997 0.56049903   *
#> Z12          8.860307e-05 -0.07884948 0.08519691    
#> Z21          2.186006e-02 -0.06120048 0.11225471    
#> Z22          7.478997e-01  0.66070564 0.84609196   *
#> A1          -2.847638e-02 -0.09220293 0.03182853    
#> A2          -2.734502e-02 -0.09837143 0.07710229    
#> Z11:A1       1.254810e-01  0.04402056 0.20808459   *
#> Z12:A1       2.334790e-01  0.14464795 0.32167582   *
#> Z21:A2       9.672511e-01  0.89485346 1.05950271   *
#> Z22:A2       5.856729e-01  0.49868921 0.65272278   *
#> A1:A2        8.014366e-01  0.76122619 0.84148385   *
```
