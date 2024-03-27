#' Inference under nonregrularity with M-out-of-N cluster bootstrap 
#' 
#' @description A function that use clustered Q-learning to produce estimation and inference on regression parameters for two stages. It uses M-out-of-N cluster bootstrap (M resample out of N samples) to adjust for potential nonregrularity for the first stage inference.
#'
#' @param completeData A matrix contains all data including covariates, outcome and treatment
#' @param s1Formula A formula that specifying the regression model for stage 1 Q-function
#' @param s2Formula A formula that specifying the regression model for stage 2 Q-function
#' @param s2Treat A string that records variable name of stage 2 treatment in s2_data
#' @param cluster A string that records variable name of the cluster identifier
#' @param bootNum A number indicates number of bootstrap sampling in constructing CIs
#' @param alpha A number indicates the type I error of inference, the default value is 0.05
#' @param nu A number indicates the desired global type I error rate for t-test, the default value is 0.05
#' @param lambda A number between 0 and 1 that governs the smallest acceptable resample size. Typical choice are between [0.025, 0.1], default value is 0.025 
#' @param ... Additional parameters in lm or lme
#' 
#' @return A list containing 
#' \itemize{
#'   \item{s1Coefficients:}{ Stage 1 regression coefficients}
#'   \item{s2Coefficients:}{ Stage 2 regression coefficients}
#'   \item{s1Inference:}{ Stage 1 coefficients and confidence interval}
#'   \item{s2Inference:}{ Stage 2 coefficients and confidence interval}
#'   \item{estM:}{ Estimated stage 1 bootstrap resample size}
#' }
#' 
#' @importFrom MASS mvrnorm
#' @importFrom stats model.matrix as.formula lm coef quantile model.frame 
#' @importFrom dplyr %>% mutate case_when
#' @importFrom nlme lme fixef corCompSymm
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
#' df = mutate(df, Z21 = rep(rbinom(N, 1, expit(delta["d1"] * Z11 + delta["d2"] * A1)), each = ni))
#' df = mutate(df, Z22 = rep(rbinom(N, 1, expit(delta["d3"] * Z12 + delta["d4"] * A1)), each = ni))
#' df = mutate(df, Y =  gamma["g0"] +
#'            gamma["g1"] * X11 + gamma["g2"] * X12 +
#'            gamma["g3"] * X21 + gamma["g4"] * X22 +
#'            gamma["g5"] * Z11 + gamma["g6"] * Z12 +
#'            gamma["g7"] * Z21 + gamma["g8"] * Z22 +
#'            gamma["g9"] * A1 + gamma["g10"] * A2 +
#'            gamma["g11"] * Z11 * A1 + gamma["g12"] * Z12 * A1 +
#'            gamma["g13"] * Z21 * A2 + gamma["g14"] * Z22 * A2 +
#'            gamma["g15"] * A1 * A2 + as.numeric(cluster_residuals))
#'   
#' # stage 1 model
#' Formula1 = formula(Y ~ X11 + X12 + Z11 + Z12 + A1 + Z11 * A1 + Z12 * A1)
#' 
#' # stage 2 model
#' Formula2 = formula(Y ~ X11 + X12 + X21 + X22 + Z11 + Z12 + Z21 + Z22 + A1 + A2 + 
#' Z11 * A1 + Z12 * A1 + Z21 * A2 + Z22 * A2 + A1 * A2)
#' result = clusterQ_MN(completeData = df,s1Formula = Formula1,s2Formula = Formula2, s2Treat = "A2",cluster = "cluster_id", bootNum =100)
#' result$s1Inference # check stage 1 inference
#' result$s2Inference # check stage 2 inference
#' @export
clusterQ_MN <- function(completeData,
                s1Formula,# s1Contrast,# not necessary delete
                s2Formula,# s2Contrast,
                s2Treat,# interact,# s2Indicator,# (delete)
                cluster, 
                bootNum =1000,
                alpha=0.05,
                nu = 0.05,
                lambda = 0.025,...) {
  
  interact <- interact.terms(s2Formula, s2Treat)
  colnames(completeData)[colnames(completeData) == cluster] <- "cluster" # change colname of cluster indicator for using lme function
  n_cluster = length(unique(completeData$cluster))
  
  s2_data<-subset(completeData, !is.na(get(s2Treat))) # YS edited
  n1 <- dim(completeData)[1]
  n2 <- dim(s2_data)[1]
  p <- nonreg(s2Formula,s2_data,s2Treat,cluster,nu) 
  M <- estM(n_cluster, p, lambda)
  
  s2_var<-names(as.data.frame(model.matrix(s2Formula,data=s2_data,...)))[-1]
  interaction<-paste(interact,s2Treat,sep=":")
  interaction2<-paste(s2Treat,interact,sep=":")
  reverse_match<-match(interaction2,s2_var)
  if (sum(!is.na(reverse_match))>0) {
    match_num<-which(!is.na(reverse_match))
    interaction[match_num]<-interaction2[match_num]
  }
  s2_var_noint<-s2_var[-match(s2Treat,s2_var)]
  s2_var_noint<-s2_var_noint[-match(interaction,s2_var_noint)]
  peudo_formula<- as.formula(paste(" ~ ", paste(s2_var_noint, collapse= "+")))
  MM <- model.matrix(peudo_formula,data=completeData)
  a<-match(names(as.data.frame(MM)),names(completeData))
  
  if (sum(is.na(a))>0) {
    MM <- MM[match(rownames(completeData),rownames(MM)),] 
    rownames(MM) <- rownames(completeData)
    extra<-MM[,is.na(a)]
    completeData<-cbind(completeData,extra)
  }
  
  ############## Stage 2 analysis  ###############
  X2<-model.matrix(s2Formula,data=s2_data,...)
  s2Contrast<-diag(ncol(X2))
  stage2.model <- lme(s2Formula, data = s2_data,
               random = ~ 1 | cluster, 
               correlation = corCompSymm(form = ~ 1 | cluster), 
               method = "REML",...) # Use REML for estimation
  s2_cf <- fixef(stage2.model)

  # Initialization
  stage2 <- matrix(0,nrow=dim(t(s2Contrast))[1],ncol=3)
  bootest2 <- matrix(0,nrow=dim(t(s2Contrast))[1],ncol=bootNum)
  for (i in 1:bootNum) {
    index <- sample(1:n2,n2,replace=TRUE)
    bootsamp <- s2_data[index,]
    stage2.model.boot <- lme(s2Formula, data = bootsamp,
                        random = ~ 1 | cluster, 
                        correlation = corCompSymm(form = ~ 1 | cluster), 
                        method = "REML",...) # Use REML for estimation
    s2_cf.boot <- fixef(stage2.model.boot)
    bootest2[,i] <- s2Contrast%*%(2*s2_cf - s2_cf.boot)
  }
  for (i in 1:nrow(t(s2Contrast))) {
    stage2[i,2:3] <- quantile(bootest2[i,], probs=c(alpha/2,1-alpha/2), na.rm=TRUE)
  }
  stage2[,1]<-s2Contrast%*%s2_cf
  colnames(stage2) <- c("S2_Estimator", "Lower", "Upper")
  rownames(stage2) <- names(s2_cf) 
  stage2 <- as.data.frame(stage2)
  stage2 <- round(stage2,4)
  stage2 <- stage2 %>%
    mutate(sig= case_when(Lower*Upper > 0 ~ "*",
                          TRUE ~ ""))
  
  ############## Stage 1 analysis  ###############
  s1_var<-names(model.frame(s1Formula,data=completeData,...))
  X1<-model.matrix(s1Formula,data=completeData,...)
  s1Contrast<-diag(ncol(X1))

  # Construct pseudo outcome
  peudo_data<-completeData
  peudo_data[,s1_var[1]]<-peudo_data[,s1_var[1]]+s2_cf[1]+
    as.matrix(peudo_data[,s2_var_noint])%*%s2_cf[s2_var_noint]+
    abs(s2_cf[s2Treat]+
          as.matrix(peudo_data[,interact])%*%s2_cf[interaction])
  stage1.model <- lme(s1Formula, data = peudo_data,
                      random = ~ 1 | cluster, 
                      correlation = corCompSymm(form = ~ 1 | cluster), 
                      method = "REML",...) # Use REML for estimation
  s1_cf <- fixef(stage1.model)

  # Initialization
  stage1 <- matrix(0,nrow=dim(t(s1Contrast))[1],ncol=3)
  bootest1 <- matrix(0,nrow=dim(t(s1Contrast))[1],ncol=bootNum)
  k <- 0
  while ( k < bootNum ) {
    index <- sample(1:n_cluster,M,replace=TRUE)
    bootsamp <- completeData[0,]
    for (i in index) {
      bootsamp <- rbind(bootsamp, completeData[completeData$cluster == i,])
    }
    bootsamp_s2 <- subset(bootsamp, !is.na(get(s2Treat))) # YS edited
    stage2.model.boot <- lme(s2Formula, data = bootsamp_s2,
                        random = ~ 1 | cluster, 
                        correlation = corCompSymm(form = ~ 1 | cluster), 
                        method = "REML",...) # Use REML for estimation
    stage2cf <- fixef(stage2.model.boot)
    bootsamp[,s1_var[1]]<-bootsamp[,s1_var[1]]+stage2cf[1]+
      as.matrix(bootsamp[,s2_var_noint])%*%stage2cf[s2_var_noint]+
      abs(stage2cf[s2Treat]+as.matrix(bootsamp[,interact])%*%stage2cf[interaction])
    
    if (sum(is.na(bootsamp[,s1_var[1]])) < (M * n1 / n_cluster)) {
      k <- k+1
      stage1.model.boot <- lme(s1Formula, data = bootsamp,
                          random = ~ 1 | cluster, 
                          correlation = corCompSymm(form = ~ 1 | cluster), 
                          method = "REML",...) # Use REML for estimation
      stage1cf <- fixef(stage1.model.boot)
      bootest1[,k] <- s1Contrast%*%(2*s1_cf-stage1cf)
    }
  }
  
  for (i in 1:nrow(t(s1Contrast))) {
    stage1[i,2:3] <- quantile(bootest1[i,], probs=c(alpha/2,1-alpha/2), na.rm=TRUE)
  }
  stage1[,1]<-s1Contrast%*%s1_cf
  colnames(stage1)<-c("S1_Estimator", "Lower", "Upper")
  rownames(stage1) <- names(s1_cf) 
  stage1 <- as.data.frame(stage1)
  stage1 <- round(stage1,4)
  stage1 <- stage1 %>%
    mutate(sig= case_when(Lower*Upper > 0 ~ "*",
                          TRUE ~ ""))
  
  cat("The estimated degree of nonregularity for stage 1 is", p, "\n")
  cat("chosen value of M =", M, "out of N =", n_cluster , "clusters.", "\n")
  object <- list(s1Coefficients = s1_cf, 
                 s1Inference = stage1,
                 s2Coefficients = s2_cf, 
                 s2Inference = stage2,
                 estM = M) 
  object$call <- match.call()
  object
}

