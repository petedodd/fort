## various utilities that are not exported
##' Hi Bounds
##'
##' 97.5 percent quantile
##' @title hi
##' @param x a vector of inputs
##' @return the 97.5% quantile
##' @author Pete Dodd
hi <- function(x)quantile(x,0.975)

##' Lo Bounds
##'
##' 2.5 percent quantile
##' @title lo
##' @param x a vector of inputs
##' @return the 2.5% quantile
##' @author Pete Dodd
lo <- function(x)quantile(x,0.025)

##' Logit Transformation
##'
##' Taken from `logitnorm` package
##' @title logit
##' @param x a vector of inputs
##' @return logit(x)
##' @import logitnorm
##' @author Pete Dodd
logit <- function(x) logitnorm::logit(x)

##' Inverse Logit Transformation
##'
##' Taken from `logitnorm` package
##' @title expit
##' @param x a vector of inputs
##' @return 1/(1+exp(-x))
##' @import logitnorm
##' @author Pete Dodd
expit <- function(x) logitnorm::invlogit(x)


##' Parametrize Lognormal Distribution
##'
##' Lognormal distribution parameters from moment-matching
##' @title lnparz
##' @param mn mean of target 
##' @param S SD of target
##' @return list with mu (mu) and sigma (sdlog) for lognormal distribution
##' @author Pete Dodd
lnparz <- function(mn,S){ #get log normal parms
  list(mu=log(mn^2/sqrt(mn^2+S^2)),sdlog=sqrt(log(1+S^2/mn^2)))
}

##' Parametrize Lognormal Distribution
##'
##' Lognormal distribution parameters from moment-matching
##' @title getLNmu
##' @param m mean of target 
##' @param S SD of target
##' @return  mu for lognormal distribution
##' @author Pete Dodd
getLNmu <- function(m,S){
  log(m^2/sqrt(m^2+S^2))
}

##' Parametrize Lognormal Distribution
##'
##' Lognormal distribution parameters from moment-matching
##' @title getLNsig
##' @param m mean of target 
##' @param S SD of target
##' @return  sigma (sdlog) for lognormal distribution
##' @author Pete Dodd
getLNsig <- function(m,S){
  sqrt(log(1+S^2/m^2))
}

##' Lognormal Distribution Mean and SD
##'
##' Compute the mean and SD for a lognormal distribution with given parameters
##' @title exparz
##' @param mu lognormal distribution mu parameter
##' @param sdlog lognormal distribution sigma parameter
##' @return list with mn (mean) and S (SD) of associated lognormal distribution
##' @author Pete Dodd
exparz <- function(mu,sdlog){
  list(mn=exp(mu+sdlog^2/2),S=sqrt(exp(sdlog^2)-1)*exp(mu+sdlog^2/2))
}

##' Parametrize Logitnormal Distribution
##'
##' Parametrize logitnormal distribution based on mean/SD. Wrapper for `logitnorm::twCoefLogitnormCi`
##' @title lgtparz
##' @param mn mean of target
##' @param S SD of target
##' @import logitnorm
##' @return list with mu (mu) and sigma (sdlog) parameters for logitnormal distribution
##' @author Pete Dodd
lgtparz <- function(mn,S){ #get logit normal parms
  A <- cbind(mn,S)
  notNA <- !is.na(mn+S)
  U <- pmin(mn[notNA]+1.0*S[notNA],1)
  L <- pmax(mn[notNA]-1.0*S[notNA],0)
  A[notNA,] <- logitnorm::twCoefLogitnormCi(lower=L,upper=U,perc = 0.6827) #w/i +/- SD
  A[!notNA,] <- rep(NA_real_,2)
  list(mu=A[,1],sdlog=A[,2])
}
