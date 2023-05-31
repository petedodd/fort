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
##' @param trnsfm 0 = no transformation; 1 = log transformation; 2 = logit transformation
##' @return data.table of results
##' @author Pete Dodd
##' @import data.table
##' @import imputeTS
noisyex <- function(yrz,mnz,sdz,nrep,runs=TRUE,trnsfm=0){
  if(trnsfm==1){
    tf <- lnparz(mnz,sdz)
    mnz <- tf$mu
    sdz <- tf$sdlog
  }
  if(trnsfm==2){
    tf <- lgtparz(mnz,sdz)
    mnz <- tf$mu
    sdz <- tf$sdlog
  }
  ## add noise
  Xreps <- matrix(mnz,nrow=length(mnz),ncol=nrep)
  upto <- which(!is.na(Xreps[,1]))
  for(i in 1:nrep) Xreps[upto,i] <- Xreps[upto,i] + rnorm(length(upto),sd=sdz[upto])
  ## imput
  imp <- imputeTS::na_kalman(Xreps)
  if(trnsfm==1) imp <- exp(imp)
  if(trnsfm==2) imp <- expit(imp)
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
##' @param TXf Treatment fatality midpoint (NA if projection/imputation needed, assumed 0 if not given)
##' @param HRd Hazard ratios to apply to detection hazard (assumed 1 if not given)
##' @param HRi Hazard ratios to apply to incidence (assumed 1 if not given)
##' @param ORt Hazard ratios to apply to treatment fatality (assumed 1 if not given)
##' @param nrep number of replicates used - NOTE likely to be temporary
##' @param output Return data + projection ('projection') or fit + projection ('fit')
##' @param modeltype Which version to use:
##'   * `failsafe': the failsafe model using simulation independent timeseries models
##' @return A data.frame/data.table with the projections and uncertainty
##' @examples
##'
##' # Some fake data for basic examples:
##' 
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
##' # Default projections using failsafe method:
##'
##' ans0 <- projections(tmp$year,tmp$Ihat,tmp$sEI,tmp$Nhat,tmp$sEN,tmp$Mhat,tmp$sEM)
##'
##'
##' # Create some treatment CFR data:
##' tmp$TXf <- tmp$Nhat/1e3
##'
##' # Failsafe run using treatment CFR data:
##' ans1 <- projections(year=tmp$year,
##'          Ihat=tmp$Ihat,sEI=tmp$sEI,
##'          Nhat=tmp$Nhat,sEN=tmp$sEN,
##'          Mhat=tmp$Mhat,sEM=tmp$sEM,
##'          Phat=tmp$Mhat,sEP=tmp$sEM,
##'          TXf = tmp$TXf)
##'
##' # Create intervention data:
##' HRd <- HRi <- ORt <- rep(1,length(tmp$Nhat))
##' HRd[is.na(tmp$Nhat)] <- HRi[is.na(tmp$Nhat)] <- 1.1
##' ORt[is.na(tmp$Nhat)] <- 0.7
##'
##' # Failsafe run with interventions:
##' ans2 <- projections(year=tmp$year,
##'          Ihat=tmp$Ihat,sEI=tmp$sEI,
##'          Nhat=tmp$Nhat,sEN=tmp$sEN,
##'          Mhat=tmp$Mhat,sEM=tmp$sEM,
##'          Phat=tmp$Mhat,sEP=tmp$sEM,
##'          TXf = tmp$TXf,ORt=ORt,
##'          HRd=HRd,HRi=HRi)
##' ans2
##'
##' @author Pete Dodd
##' @import data.table
##' @export
projections <- function(year,
                        Ihat,sEI,
                        Nhat,sEN,
                        Mhat,sEM,
                        Phat,sEP,
                        ...,
                        nrep=500,
                        output='projection', #TODO
                        modeltype='failsafe'
                        ){
  ## argument management
  arguments <- list(...)
  list2env(arguments,envir = environment())                   #boost ... to this scope
  ## completions
  if(!'TXf' %in% names(arguments)) TXf <- rep(0,length(year))
  if(!'HRd' %in% names(arguments)) HRd <- rep(1,length(year))
  if(!'HRi' %in% names(arguments)) HRi <- rep(1,length(year))
  if(!'ORt' %in% names(arguments)) ORt <- rep(1,length(year))
  if(modeltype=='failsafe'){ #============== FAILSAFE MODEL ==============
    ## incidence
    suppressWarnings({RI <- noisyex(year,Ihat,sEI,nrep,runs=FALSE,1)})
    names(RI) <- c('year','I.mid','I.sd','I.lo','I.hi')
    ## notifications
    suppressWarnings({RN <- noisyex(year,Nhat,sEN,nrep,runs=FALSE,1)})
    names(RN) <- c('year','N.mid','N.sd','N.lo','N.hi')
    ## mortality
    suppressWarnings({RM <- noisyex(year,
                                    Mhat-TXf*Nhat, #take off mortality on treatment (see below addon)
                                    sEM,nrep,runs=FALSE,1)})
    names(RM) <- c('year','M.mid','M.sd','M.lo','M.hi')
    ## fraction HIV+
    if('Hhat' %in% names(arguments) & 'sEH' %in% names(arguments)){ #HIV
      cat('Running HIV component...\n')
      suppressWarnings({RH <- noisyex(year,arguments$Hhat,arguments$sEH,nrep,runs=FALSE,2)})
      names(RH) <- c('year','H.mid','H.sd','H.lo','H.hi')
    } else { #HIV not given
      RH <- data.table::data.table(year=year,
                                   H.mid=rep(0.0,length(year)),
                                   H.sd=rep(0.0,length(year)),
                                   H.lo=rep(0.0,length(year)),
                                   H.hi=rep(0.0,length(year)))
    }
    ## prevalence NOTE not currently implemented
    RP <- data.table::data.table(year=year,
                                 P.mid=rep(NA_real_,length(year)),
                                 P.sd=rep(NA_real_,length(year)),
                                 P.lo=rep(NA_real_,length(year)),
                                 P.hi=rep(NA_real_,length(year)))
    ## merge and output
    ANS1 <- merge(RI,RN,by='year')
    ANS2 <- merge(RM,RP,by='year')
    ANS <- merge(ANS1,ANS2,by='year')
    if('Hhat' %in% names(arguments) & 'sEH' %in% names(arguments)) ANS <- merge(ANS,RH,by='year')
    ## HR interventions NOTE the notif one will have limited validity
    ANS[,c('I.mid','I.lo','I.hi'):=.(I.mid*HRi,I.lo*HRi,I.hi*HRi)] #OK
    ANS[,c('N.mid','N.lo','N.hi'):=.(N.mid*HRi*HRd,N.lo*HRi*HRd,N.hi*HRi*HRd)] #approx
    ANS[,c('M.mid','M.lo','M.hi'):=.(M.mid/HRd,M.lo/HRd,M.hi/HRd)]             #approx
    ## add treated mortality back
    ## project TXf
    if(any(is.na(TXf))){
      TXf <- expit( imputeTS::na_kalman(logit(TXf)) + log(ORt) )
    } else { TXf <- expit( logit(TXf) + log(ORt) );} #apply effect
    ANS[,c('M.mid','M.lo','M.hi'):=.(M.mid + TXf*N.mid,M.lo + TXf*N.mid,M.hi + TXf*N.mid)] #NOTE no extra uncertainty
    ## computing this using duration assumption -
    ## WHO methods appendix: tx ~ U[0.2,2]; ut ~ U[1,4]
    tx.mid <- (2+0.2)/2; ut.mid <- (4+1)/2 #midpoints
    tx.sd <- (2-0.2)/3.92; ut.sd <- (4-1)/3.92 #SD
    ANS[,P.mid:=N.mid*tx.mid + (I.mid-N.mid)*ut.mid]
    ## P.sd^2 = (N.mid*tx.m)^2 * ((N.sd/N.mid)^2+(tx.sd/tx.mid)^2) +
    ##     ((I.mid-N.mid)*ut.m)^2 * ( (ut.sd/ut.m)^2 + (I.sd^2+N.sd^2)/(I.mid-N.mid)^2 )
    ANS[,P.sd:=sqrt(
    (N.mid*tx.mid)^2 * ((N.sd/N.mid)^2+(tx.sd/tx.mid)^2)+
    ((I.mid-N.mid)*ut.mid)^2 * ( (ut.sd/ut.mid)^2 + (I.sd^2+N.sd^2)/(I.mid-N.mid)^2 )
    )]
    ANS[,c('P.lo','P.hi'):=.(pmax(0,P.mid-1.96*P.sd),P.mid+1.96*P.sd)]
  } else {
    stop(paste0("modeltype = ",modeltype," is not defined!"))
  }
  ## output
  ANS
}




##' Main projection function to be exported
##'
##' .. content for \details{Cprojections} ..
##' @title Cprojections
##' @param year vector of years for which there is data or an estimate is required
##' @param Ihat Incidence midpoints (NA if projection/imputation needed)
##' @param sEI Incidence uncertainty as SD (NA if projection/imputation needed)
##' @param Nhat Notifications midpoints (NA if projection/imputation needed)
##' @param sEN Notifications uncertainty as SD (NA if projection/imputation needed)
##' @param Mhat Deaths midpoints (NA if projection/imputation needed)
##' @param sEM Deaths uncertainty as SD (NA if projection/imputation needed)
##' @param Phat Prevalence midpoints (NA if projection/imputation needed)
##' @param sEP Prevalence uncertainty as SD (NA if projection/imputation needed)
##' @return A data.frame/data.table with the projections and uncertainty
##' @examples
##' 
##' @author Pete Dodd
##' @import data.table
##' @import bssm
##' @export
Cprojections <- function(year,
                        Ihat,sEI,
                        Nhat,sEN,
                        Mhat,sEM,
                        Phat,sEP,
                        ...
                        ){

  pntrs <- create_xptrs() #create pointers
  initial_theta <- c(logSI = 1)

  known_params <- c(mI0 = Ihat[1],mP0 = Phat[1],
                    mN0 = Nhat[1],mD0 = Mhat[1],
                    momega = -1.1, mdelta = 0, mpsi = logit(0.5),
                    sI0 = Ihat[1]/1e2,sP0 = Phat[1]/1e2,
                    sN0 = Nhat[1]/1e2,sD0 = Mhat[1]/1e2,
                    somega = 0.7,sdelta = sdelta0,spsi = 0.3)

  known_tv_params <- cbind(sEI,sEP,
                           sEN,sEM)

  ## tests:
  state <- c(100, 100, 100, 10,
             log(3),log(1),logit(0.5))

  ## create model
  model <- bssm::ssm_nlg(y = Yhat, a1=pntrs$a1, P1 = pntrs$P1, 
                   Z = pntrs$Z_fn, H = pntrs$H_fn, T = pntrs$T_fn, R = pntrs$R_fn, 
                   Z_gn = pntrs$Z_gn, T_gn = pntrs$T_gn,
                   theta = initial_theta, log_prior_pdf = pntrs$log_prior_pdf,
                   known_params = known_params, known_tv_params = known_tv_params,
                   n_states = 7, n_etas = 4,
                   state_names = SNMZ)

  ## inference
  mcmc.type <- 'ekf'              #change inference type
  ITER <- 6000 ; BURN <- 1000 
  mcmc_fit <- bssm::run_mcmc(model, iter = ITER, burnin = BURN,mcmc_type = "ekf")
  ## summary(mcmc_fit, return_se = TRUE)
  outs <- data.table::as.data.table(mcmc_fit,variable='states')
  outs <- outs[,.(mid=mean(value),lo=quantile(value,0.0275),hi=quantile(value,0.975)),
               by=.(variable,time)]
  outs

}
