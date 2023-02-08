## various utilities that are not exported
hi <- function(x)quantile(x,0.975)
lo <- function(x)quantile(x,0.025)
logit <- function(x)log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))
gh <- function(x) glue(here(x))
lnparz <- function(mn,S){
  list(mu=log(mn^2/sqrt(mn^2+S^2)),sdlog=sqrt(log(1+S^2/mn^2)))
}
exparz <- function(mu,sdlog){
  list(mn=exp(mu+sdlog^2/2),S=sqrt(exp(sdlog^2)-1)*exp(mu+sdlog^2/2))
}
