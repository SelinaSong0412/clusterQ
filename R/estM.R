#' Estimate the resample size M from N sample
#' 
#' @description A function that estimate the resample size M from N sample with given degree of nonregularity in the stage and a prespecified parameter \code{lambda}
#'
#' @param N A number of total sample size
#' @param p A number of degree of nonregularity
#' @param lambda A number between 0 and 1 that governs the smallest acceptable resample size. Typical choice are between [0.025, 0.1], default value is 0.025 
#' 
#' @return \item{M:}{ An estimated resample size}
#' 
#' @examples 
#' M_hat <- estM(N = 100, p = 0.45, lambda = 0.02) 
#' M_hat
#' 
#' @export

estM <- function(N, p, lambda = 0.025) {
  M <- ceiling(N^(1-p*(lambda/(1+lambda))))
  return(M)
}
