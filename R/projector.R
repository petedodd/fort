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
##' @examples
##' tmp <- structure(list(year = c(2000, 2001, 2002, 2003, 2004, 2005, 2006,
##' 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017,
##' 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026),
##' Ihat = c(288.533950004829,287.93322935538, 286.858763024248, 285.088943283027, 282.422826602471,
##' 278.659650849076, 273.695148906889, 267.782862966495, 261.238262319895,
##' 254.351459991824, 247.410089922074, 240.673040367854, 234.284291413532,
##' 228.355885062355, 222.940587325497, 216.730763815714, 210.520940305932,
##' 204.489042196973, 198.629971526196, 192.938776399049, NA, NA,NA, NA, NA, NA, NA),
##' sEI = c(82.4505044782574, 82.2788459607041,81.971812555474, 81.4660788612642,
##' 80.7042243167865, 79.6288801248133, 78.2102512754221, 76.5207883326334,
##' 74.6506378711457, 72.6827008981588, 70.6991701325661, 68.7740241983581,
##' 66.9484062245379, 65.2543329198881,50.0532444478606, 45.4237612044787,
##' 39.3249190909355, 36.161588804207,35.1254802829674, 34.1190584700803,
##' NA, NA, NA, NA, NA, NA, NA),
##' Nhat = c(105.597560166138, 100.937200476189, 97.0396342996481,
##' 96.5595726784802, 100.580594702341, 100.752701403095, 105.434702191619,
##' 109.52777439355, 110.960319162264, 111.019454219404, 108.554358614967,
##' 105.891527759511, 101.900468674008, 97.1161848558792, 124.231710859869,
##' 127.247487405267, 133.171236539199, 123.233182638182, 141.107743265778,
##' 158.247577690289, NA, NA, NA, NA, NA, NA, NA),
##' sEN = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, NA,
##' NA, NA, NA, NA, NA),
##' Mhat = c(67.3923314966283, 64.9639449981795,
##' 62.7591216508678, 59.9862798209979, 56.8774350966658, 55.2906791480839,
##' 53.3253778160102, 52.5869006070828, 49.256969022149, 47.0520258424162,
##' 44.4951773213113, 42.955216531147, 41.6499649964408, 40.4263075158909,
##' 37.4464515389009, 35.8778532613629, 34.6111975994581, 34.0979910751733,
##' 33.4345311257622, 32.5925643081621, NA, NA, NA, NA, NA, NA, NA),
##' sEM = c(5.50426791938695, 4.99926659859985, 4.62031120755379,
##' 4.27306764998549, 3.9284689776149, 3.88167033820785, 3.73219080309218,
##' 2.33566032352505, 3.18594599695718, 2.89690790680736, 2.67494608005328,
##' 2.52339290625221, 2.49578415017652, 1.43500737514435, 1.10097123857461,
##' 1.14876519324219, 1.22099474332274, 1.21995493069537, 1.2210104574169,
##' 1.21500236771799, NA, NA, NA, NA, NA, NA, NA)), row.names = c(NA,-27L), class = "data.frame")
##'
##' test <- projections(tmp$year,tmp$Ihat,tmp$sEI,tmp$Nhat,tmp$sEN,tmp$Mhat,tmp$sEM)
##' 
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
