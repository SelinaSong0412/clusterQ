#' Estimate degree of nonregularity
#' 
#' @description A function that estimate the degree of nonregularity using the next stage information.
#'
#' @param s2Formula A formula that specifying the regression model for Q-function of the next stage.
#' @param s2_data A matrix contains the next stage data including covariates, outcome and treatment
#' @param s2Treat A string that records variable name of next stage treatment in `s2_data`
#' @param cluster A string that records variable name of the cluster identifier
#' @param nu A number that indicating the desired global type I error for t-test; default value is 0.05
#' @param ... Additional parameters in lm
#' 
#' @return \item{p:}{ A scalar of estimated degree of nonregularity at stage 1}
#' 
#' @importFrom MASS mvrnorm
#' @importFrom stats model.matrix lm coef model.frame model.response vcov qt
#' @importFrom nlme lme fixef corCompSymm
#' @importFrom dplyr %>% distinct select
#' @importFrom tidyselect all_of
#' 
#' @examples 
#' library(dplyr)
#' expit <- function(x) {return(exp(x)/ (1 + exp(x)))}
#' N = 100
#' n = 2000
#' ni = n/N
#' rho = 0.4
#' 
#' # Coefficients for the outcome model
#' gamma <- c(g0 = 1, g1 = 0.1, g2 = 0.2, g3 = 0.3, g4 = 0.4,
#'            g5 = 0.5, g6 = 0.6, g7 = 0.7, g8 = 0.8, g9 = 0.9,
#'            g10 = 1, g11 = 0, g12 = 0, g13 = 0, g14 = 0,
#'            g15 = 0)
#' 
#' # Parameters for Z21 and Z22 simulation
#' delta <- c(d1 = 0.5, d2 = 0.3, d3 = 0.4, d4 = 0.2)
#' 
#' set.seed(412) # For reproducibility
#' 
#' # Simulate Stage 1 data
#' df <- data.frame(
#'   cluster_id = rep(1:N, each = ni),
#'   X11 = rnorm(n), # Stage 1 patient-level data
#'   X12 = rnorm(n), # Another Stage 1 patient-level variable
#'   X21 = rnorm(n), # Stage 2 patient-level data
#'   X22 =  rnorm(n), # Another Stage 2 patient-level variable
#'   Z11 = rep(rbinom(N, 1, 0.5), each = ni),
#'   Z12 = rep(rbinom(N, 1, 0.5), each = ni),
#'   A1 = rep(2*rbinom(N, 1, 0.5)-1, each = ni),
#'   A2 = rep(2*rbinom(N, 1, 0.5)-1, each = ni)
#' )
#' # Generating residuals with intra-cluster correlation (ICC)
#' Sigma <- matrix(rho, N, N) + diag(1-rho, N)  # Exchangeable correlation matrix
#' cluster_residuals <- MASS::mvrnorm(n = ni, mu = rep(0, N), Sigma = Sigma)
#' 
#' df <- df %>% 
#'   mutate(Z21 = rep(rbinom(N, 1, expit(delta["d1"] * Z11 + delta["d2"] * A1)), each = ni),
#'          Z22 = rep(rbinom(N, 1, expit(delta["d3"] * Z12 + delta["d4"] * A1)), each = ni),
#'          Y =  gamma["g0"] +
#'            gamma["g1"] * X11 + gamma["g2"] * X12 +
#'            gamma["g3"] * X21 + gamma["g4"] * X22 +
#'            gamma["g5"] * Z11 + gamma["g6"] * Z12 +
#'            gamma["g7"] * Z21 + gamma["g8"] * Z22 +
#'            gamma["g9"] * A1 + gamma["g10"] * A2 +
#'            gamma["g11"] * Z11 * A1 + gamma["g12"] * Z12 * A1 +
#'            gamma["g13"] * Z21 * A2 + gamma["g14"] * Z22 * A2 +
#'            gamma["g15"] * A1 * A2 + as.numeric(cluster_residuals)
#'   )
#' # stage 1 model
#' Formula1 = formula(Y ~ X11 + X12 + Z11 + Z12 + A1 + Z11 * A1 + Z12 * A1)
#' 
#' # stage 2 model
#' Formula2 = formula(Y ~ X11 + X12 + X21 + X22 + Z11 + Z12 + Z21 + Z22 + A1 + A2 + Z11 * A1 + Z12 * A1 + Z21 * A2 + Z22 * A2 + A1 * A2)
#' p = nonreg(s2Formula = Formula2, s2_data = df, s2Treat = "A2", cluster = "cluster_id", nu = 0.05)
#' p
#' 
#' @export

nonreg <- function(s2Formula,s2_data,
                   s2Treat, cluster, nu = 0.05,...) {
  
  interact <- interact.terms(s2Formula, s2Treat)
  n2 = nrow(s2_data)
  # Yao modify
  colnames(s2_data)[colnames(s2_data) == cluster] <- "cluster"
  n_cluster = length(unique(s2_data$cluster))
  
  stage2.model <- lme(s2Formula, data = s2_data,
                      random = ~ 1 | cluster, 
                      correlation = corCompSymm(form = ~ 1 | cluster), 
                      method = "REML",...) # Use REML for estimation
  s2_cf <- fixef(stage2.model)
  Sigma2 <- vcov(stage2.model)
  
  X2 <- model.matrix(s2Formula,data=s2_data)
  s2_var<-names(as.data.frame(X2))[-1]
  interaction<-paste(interact,s2Treat,sep=":")
  interaction2<-paste(s2Treat,interact,sep=":")
  reverse_match<-match(interaction2,s2_var)
  if(sum(!is.na(reverse_match))>0) {
    match_num<-which(!is.na(reverse_match))
    interaction[match_num]<-interaction2[match_num]
  }
  
  # identifying nonreg terms
  nonreg_term <- c(s2Treat,interaction)
  nonreg_col <- match(nonreg_term,names(as.data.frame(X2)))
  Sigma_nonreg <- Sigma2[nonreg_col,nonreg_col] 
  
  H_interact_main <- s2_data %>%
    dplyr::select(cluster, all_of(interact)) %>% 
    distinct(cluster, .keep_all = TRUE) %>%
    dplyr::select(-cluster)

  Yprime <- as.numeric(s2_cf[s2Treat] + as.matrix(H_interact_main) %*%s2_cf[interaction])
  h <- cbind(1, as.matrix(H_interact_main)) # h matrix for each patients
  sigma2 <- n_cluster*diag(h%*%Sigma_nonreg%*%t(h))
  TS <- abs(Yprime)/sqrt(sigma2)
  cutoff <- qt(p = 1-nu/(2*n_cluster), df = n2/n_cluster - 1)
  
  # sigma2 <- diag(h%*%Sigma_nonreg%*%t(h))
  # TS <- abs(Yprime)/sqrt(sigma2)
  # cutoff <- qt(p = 1-nu/(2*n_cluster), df = n2/n_cluster - 1)
  
  nonregularity <- ifelse(TS <= cutoff, 1, 0)
  p <- mean(nonregularity)   # estimated value of the degree of non-regularity
  return(p)
}
