## currently slow simulation-based projector

## not exported as temporary: a utility function used to create noisy data replicates
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Simulating and fitting replicates with structural model
##' @param yrz 
##' @param mnz 
##' @param sdz 
##' @param nrep 
##' @param runs 
##' @param trnsfm 
##' @return data.table of results
##' @author Pete Dodd
##' @import data.table
##' @import imputeTS
noisyex <- function(yrz,mnz,sdz,nrep,runs=TRUE,trnsfm=FALSE){
  if(trnsfm){
    tf <- lnparz(mnz,sdz)
    mnz <- tf$mu
    sdz <- tf$sdlog
  }
  ## add noise
  Xreps <- matrix(mnz,nrow=length(mnz),ncol=nrep)
  upto <- which(!is.na(Xreps[,1]))
  for(i in 1:nrep) Xreps[upto,i] <- Xreps[upto,i] + rnorm(length(upto),sd=sdz[upto])
  ## imput
  imp <- imputeTS::na_kalman(Xreps)
  if(trnsfm) imp <- exp(imp)
  dimp <- data.table::as.data.table(imp)
  dimp$year <- yrz
  dimp <- data.table::melt(dimp,id='year')
  if(runs!=TRUE){
    dimp <- dimp[,.(value=mean(value),s=sd(value),
                    lo=lo(value),hi=hi(value)),
                 by=year]
  }
  dimp
}


##' Main projection function to be exported
##'
##' .. content for \details{projections} ..
##' @title projections
##' @param year vector of years for which there is data or an estimate is required
##' @param Ihat Incidence midpoints (NA if projection/imputation needed)
##' @param sEI Incidence uncertainty as SD (NA if projection/imputation needed)
##' @param Nhat Notifications midpoints (NA if projection/imputation needed)
##' @param sEN Notifications uncertainty as SD (NA if projection/imputation needed)
##' @param Mhat Deaths midpoints (NA if projection/imputation needed)
##' @param sEM Deaths uncertainty as SD (NA if projection/imputation needed)
##' @param Phat Prevalence midpoints (NA if projection/imputation needed)
##' @param sEP Prevalence uncertainty as SD (NA if projection/imputation needed)
##' @param nrep number of replicates used - NOTE likely to be temporary
##' @return A data.frame/data.table with the projections and uncertainty
##' @author Pete Dodd
##' @import data.table
##' @export
projections <- function(year,
                        Ihat,sEI,
                        Nhat,sEN,
                        Mhat,sEM,
                        Phat,sEP,
                        nrep=500
                        ){
  ## incidence
  suppressWarnings({RI <- noisyex(year,Ihat,sEI,nrep,runs=FALSE,TRUE)})
  names(RI) <- c('year','I.mid','I.sd','I.lo','I.hi')
  ## notifications
  suppressWarnings({RN <- noisyex(year,Nhat,sEN,nrep,runs=FALSE,TRUE)})
  names(RN) <- c('year','N.mid','N.sd','N.lo','N.hi')
  ## mortality
  suppressWarnings({RM <- noisyex(year,Mhat,sEM,nrep,runs=FALSE,TRUE)})
  names(RM) <- c('year','M.mid','M.sd','M.lo','M.hi')
  ## prevalence NOTE not currently implemented
  RP <- data.table::data.table(year=year,
                               P.mid=rep(NA_real_,length(year)),
                               P.sd=rep(NA_real_,length(year)),
                               P.lo=rep(NA_real_,length(year)),
                               P.hi=rep(NA_real_,length(year)))
  ## merge and output
  ANS1 <- merge(RI,RN,by='year')
  ANS2 <- merge(RM,RP,by='year')
  merge(ANS1,ANS2,by='year')
}
