
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
gamma <- c(g0 = 1, g1 = 0.1, g2 = 0.2, g3 = 0.3, g4 = 0.4,
           g5 = 0.5, g6 = 0.6, g7 = 0.7, g8 = 0.8, g9 = 0.9,
           g10 = 1, g11 = 0.1, g12 = 0.2, g13 = 0.3, g14 = 0.5,
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
#> [1] 0.71
```

***Calcuating the resample size at stage 1***

Use function `nonreg()` to calculate the resample size (M) from N
clusters at stage 1.

``` r
M = estM(N, p, lambda = 0.025)
M
#> [1] 93
```

***Inference with $M$-out-of-$N$ cluster bootstrap***

By specify number of boostrap sampling for constructing confidence
interval in `bootNum`, tuning parameter `lambda` for identifying
resample size, and the global type-I error `alpha`, result on estimation
and inference for both stage model can be obtained.

``` r
results = clusterQ::clusterQ_MN(completeData = df,
                      s1Formula = Formula1,
                      s2Formula = Formula2,
                      s2Treat = "A2",
                      cluster = "cluster_id", 
                      bootNum = 100,
                      lambda = 0.025,
                      alpha = 0.05)
#> The estimated degree of nonregularity for stage 1 is 0.71 
#> chosen value of M = 93 out of N = 100 clusters.
```

Check stage 1 inference by:

``` r
results$s1Inference
#>             S1_Estimator     Lower     Upper sig
#> (Intercept)    4.7089434 4.0255665 5.3066673   *
#> X11            0.2916834 0.1991565 0.3998238   *
#> X12            0.4516492 0.3133602 0.5491326   *
#> Z11            1.1899747 0.6013415 1.8372131   *
#> Z12            2.1365684 1.4754089 2.7422353   *
#> A1             2.3756589 1.7869391 3.1379072   *
#> Z11:A1         0.8090647 0.1413179 1.3766566   *
#> Z12:A1         0.8162722 0.1534394 1.4272422   *
```

Check stage 2 inference by:

``` r
results$s2Inference
#>             S2_Estimator      Lower     Upper sig
#> (Intercept)    0.9616716 0.86413356 1.0800984   *
#> X11            0.1328785 0.09653853 0.1763828   *
#> X12            0.2226868 0.17137026 0.2562611   *
#> X21            0.3090488 0.26197014 0.3564879   *
#> X22            0.4019560 0.35750736 0.4450930   *
#> Z11            0.4736197 0.39563307 0.5626672   *
#> Z12            0.6197995 0.52298225 0.7084738   *
#> Z21            0.6633120 0.58433910 0.7543839   *
#> Z22            0.8424813 0.75091983 0.9208518   *
#> A1             0.9184800 0.84147585 1.0085285   *
#> A2             0.9595998 0.88618361 1.0369087   *
#> Z11:A1         0.1152546 0.02769461 0.2150141   *
#> Z12:A1         0.1891494 0.08295130 0.2965337   *
#> Z21:A2         0.3637038 0.27759153 0.4603982   *
#> Z22:A2         0.4827637 0.39764602 0.5685890   *
#> A1:A2          0.7725062 0.72906496 0.8206921   *
```
