## various utilities that are not exported
hi <- function(x)quantile(x,0.975)
lo <- function(x)quantile(x,0.025)
## logit <- function(x)log(x/(1-x))
## expit <- function(x) 1/(1+exp(-x))
logit <- function(x) logitnorm::logit(x)
expit <- function(x) logitnorm::invlogit(x)
gh <- function(x) glue(here(x))
lnparz <- function(mn,S){ #get log normal parms
  list(mu=log(mn^2/sqrt(mn^2+S^2)),sdlog=sqrt(log(1+S^2/mn^2)))
}
getLNmu <- function(m,S){
  log(m^2/sqrt(m^2+S^2))
}
getLNsig <- function(m,S){
  sqrt(log(1+S^2/m^2))
}
exparz <- function(mu,sdlog){
  list(mn=exp(mu+sdlog^2/2),S=sqrt(exp(sdlog^2)-1)*exp(mu+sdlog^2/2))
}
lgtparz <- function(mn,S){ #get logit normal parms
  A <- cbind(mn,S)
  notNA <- !is.na(mn+S)
  U <- pmin(mn[notNA]+1.0*S[notNA],1)
  L <- pmax(mn[notNA]-1.0*S[notNA],0)
  A[notNA,] <- logitnorm::twCoefLogitnormCi(lower=L,upper=U,perc = 0.6827) #w/i +/- SD
  A[!notNA,] <- rep(NA_real_,2)
  list(mu=A[,1],sdlog=A[,2])
}
