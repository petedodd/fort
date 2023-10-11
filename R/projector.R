## currently slow simulation-based projector

## not exported as temporary: a utility function used to create noisy data replicates
##' Utility For Simulation-based Structural Time Series Fitting & Projection
##'
##' Samples timeseries within uncertainty bounds and uses a Kalman filter from imputeTS to project.
##' @title Simulating and fitting replicates with structural model
##' @param yrz vector of years
##' @param mnz vector of means
##' @param sdz vector of SDs
##' @param nrep number of replicates to use
##' @param runs logical: whether to return individual run data
##' @param trnsfm 0 = no transformation; 1 = log transformation; 2 = logit transformation
##' @return data.table of results
##' @author Pete Dodd
##' @import data.table
##' @import imputeTS
noisyex <- function(yrz,mnz,sdz,nrep,runs=TRUE,trnsfm=0){
  ## R CMD check safety
  value <- NULL
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
    dimp <- dimp[,list(value=mean(value),s=sd(value),
                    lo=lo(value),hi=hi(value)),
                 by=year]
  }
  dimp
}
##' MCMC Summary Computer
##'
##' Computes mid (mean) and hi/lo (95% quantiles) by variable and time and also adds in exponentiated variables if their names begin with log.
##' @title MCMC to summary
##' @param fit an MCMC fit object from fitting with bssm
##' @return a data.table with mid,lo,hi by time and variable
##' @import data.table
##' @import bssm
##' @author Pete Dodd
mcmcsmry <- function(fit){
  out <- data.table::as.data.table(fit,variable='states')
  tmpo <- out[grepl('log',variable)] #the logged
  tmpo[,value:=exp(value)]
  tmpo[,variable:=gsub('log','',variable)]
  out <- rbind(out,tmpo)
  out[,list(mid=mean(value),lo=lo(value),hi=hi(value)),by=list(variable,time)]
}


##' Main Projection Function
##'
##' Top level main function for use in generating projections of TB notifications, incidence, prevalence and mortality.
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
##' @param ... Additional parameters to help tune away from defaults
##' @param output Return data + projection ('projection') or fit + projection ('fit')  or futureonly
##' @param modeltype Which version to use:
##'   * `failsafe': the failsafe model using simulation independent timeseries models
##'   * `rwI': SSM with random walk for incidence. NOTE no indirect impact
##'   * `IP': SSM with AR(1) on a mixture of log(Incidence) and log(Prevalence).
##' @param returninternalfit (default FALSE) return the internal estimated states (only when returntype=='fit')
##' @param verbose (Default=FALSE) Give more output for use in debugging.
##' @return A data.frame/data.table with the projections and uncertainty
##' @importFrom stats deltat frequency quantile rnorm sd ts tsp
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
##' \dontrun{
##' ans0 <- projections(tmp$year,tmp$Ihat,tmp$sEI,tmp$Nhat,tmp$sEN,tmp$Mhat,tmp$sEM)
##' }
##'
##'
##' # Create some treatment CFR data:
##' tmp$TXf <- tmp$Nhat/1e3
##'
##' # Failsafe run using treatment CFR data:
##' \dontrun{
##' ans1 <- projections(year=tmp$year,
##'          Ihat=tmp$Ihat,sEI=tmp$sEI,
##'          Nhat=tmp$Nhat,sEN=tmp$sEN,
##'          Mhat=tmp$Mhat,sEM=tmp$sEM,
##'          Phat=tmp$Mhat,sEP=tmp$sEM,
##'          TXf = tmp$TXf)
##' }
##' 
##' # Create intervention data:
##' HRd <- HRi <- ORt <- rep(1,length(tmp$Nhat))
##' HRd[is.na(tmp$Nhat)] <- HRi[is.na(tmp$Nhat)] <- 1.1
##' ORt[is.na(tmp$Nhat)] <- 0.7
##'
##' # Failsafe run with interventions:
##' \dontrun{
##' ans2 <- projections(year=tmp$year,
##'          Ihat=tmp$Ihat,sEI=tmp$sEI,
##'          Nhat=tmp$Nhat,sEN=tmp$sEN,
##'          Mhat=tmp$Mhat,sEM=tmp$sEM,
##'          Phat=tmp$Mhat,sEP=tmp$sEM,
##'          TXf = tmp$TXf,ORt=ORt,
##'          HRd=HRd,HRi=HRi)
##' ans2
##' }
##'
##' ## example for rwI SSM model
##' \dontrun{
##' ans3 <- projections(year=tmp$year,
##'           Ihat=tmp$Ihat,sEI=tmp$sEI,
##'           Nhat=tmp$Nhat,sEN=tmp$sEN,
##'           Mhat=tmp$Mhat,sEM=tmp$sEM,
##'           Phat=tmp$Mhat,sEP=tmp$sEM,
##'           TXf = tmp$TXf,ORt=ORt,
##'           HRd=HRd,HRi=HRi,
##'           modeltype = 'rwI',
##'           output='projection')
##'
##' ans3
##' }
##'
##' ## example for IP SSM model
##' \dontrun{
##' ans4 <- projections(year=tmp$year,
##'           Ihat=tmp$Ihat,sEI=tmp$sEI,
##'           Nhat=tmp$Nhat,sEN=tmp$sEN,
##'           Mhat=tmp$Mhat,sEM=tmp$sEM,
##'           Phat=tmp$Mhat,sEP=tmp$sEM,
##'           TXf = tmp$TXf,
##'           modeltype = 'IP',
##'           output='projection')
##'
##' ans4
##' }
##'
##' @author Pete Dodd
##' @import data.table
##' @export
projections <- function(year,
                        Ihat,sEI,
                        Nhat,sEN,
                        Mhat,sEM,
                        Phat,sEP,
                        TXf,HRd,HRi,ORt,
                        ...,
                        nrep=500,
                        output='projection', #TODO
                        modeltype='failsafe',
                        returninternalfit=FALSE,
                        verbose=FALSE
                        ){
  ## argument management
  arguments <- list(...)
  list2env(arguments,envir = environment())                   #boost ... to this scope
  ## R CMD check safety
  I.sd <- I.mid <- N.mid <- I.lo <- I.hi <- N.lo <- N.hi <- M.mid <- M.sdx3.92 <- NULL
  M.lo <-  M.hi <- P.mid <- P.sd <- N.sd <- variable <- NULL
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
    ## calculate a version of the CFR to carry fwd
    maxidxnotna <- sum(!is.na(Ihat + Nhat + Mhat)) #last index with all necessary data provided
    cfrz <- (Mhat-TXf*Nhat) / (Ihat - Nhat) #(untreated mort'y)/(untreated inc)
    ## CFR <- cfrz[maxidxnotna] #last
    CFR <- mean(cfrz,na.rm=TRUE) #mean
    if(is.na(CFR) | CFR<0 | CFR >1){
      warning('Untreated CFR implied by data provided is pathological!\nUsing CFR = 0.5')
      CFR <- 0.5                       #safety
    }
    ## calculate last CDR to carry fwd
    CDR <- Nhat[maxidxnotna] / Ihat[maxidxnotna]
    if(is.na(CDR)  | CDR<0 | CDR >1 ){
      warning('Untreated CDR implied by data provided is pathological!\nUsing CDR = 0.7')
      CDR <- 0.7                       #safety
    }
    ## replace incidence forecasts with version using N/CDR
    RI[,I.sd:=I.sd/I.mid]        #make proportion
    RI[(maxidxnotna+1):nrow(RI),
       c('I.mid'):=RN[(maxidxnotna+1):nrow(RI),list(N.mid/CDR)]] #CDR-based incidence
    RI[,I.sd:=I.sd * I.mid]        #make not a proportion
    RI[(maxidxnotna+1):nrow(RI),
       c('I.lo','I.hi'):=list(pmax(0,I.mid - 1.96 * I.sd),I.mid + 1.96 * I.sd)]
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
    ANS[,c('I.mid','I.lo','I.hi'):=list(I.mid*HRi,I.lo*HRi,I.hi*HRi)] #OK
    ANS[,c('N.mid','N.lo','N.hi'):=list(N.mid*HRi*HRd,N.lo*HRi*HRd,N.hi*HRi*HRd)] #approx
    ## ANS[,c('M.mid','M.lo','M.hi'):=list(M.mid/HRd,M.lo/HRd,M.hi/HRd)]             #approx
    ## mortality approximation: Untreated' = (Incidence' - Notifications') x CFR
    pb <- maxidxnotna+1; pe <- nrow(ANS)                           #begin/end of projection
    ANS[pb:pe,M.mid:=(I.mid-N.mid)*CFR]                            #mean
    ANS[,M.sdx3.92:=sqrt((I.lo-I.hi)^2+(N.lo-N.hi)^2)*CFR]         #uncertainty measure
    ANS[pb:pe,c('M.lo','M.hi'):=list(pmax(0,M.mid-M.sdx3.92/2),M.mid+M.sdx3.92/2)]
    ANS[,M.sdx3.92:=NULL]
    ## add treated mortality back
    ## project TXf
    if(any(is.na(TXf))){
      TXf <- expit( imputeTS::na_kalman(logit(TXf)) + log(ORt) )
    } else { TXf <- expit( logit(TXf) + log(ORt) );} #apply effect
    ANS[,c('M.mid','M.lo','M.hi'):=list(M.mid + TXf*N.mid,M.lo + TXf*N.mid,M.hi + TXf*N.mid)] # addon mortality on treatment
    ## NOTE no extra uncertainty in line above
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
    ANS[,c('P.lo','P.hi'):=list(pmax(0,P.mid-1.96*P.sd),P.mid+1.96*P.sd)]
  } else { #============== SSM versions ==============
    nahead <- which.max(!is.na(rev(Ihat)))-1 #assume NAs at back
    lastd <- length(Ihat)-nahead
    ## TODO preprocess out NAs by interpolation if needed?
    logIRR <- log(HRi[(lastd+1):length(HRd)])      #IRR on incidence
    logIRRdelta <- log(HRd[(lastd+1):length(HRd)]) #detection
    logORpsi <- log(ORt[(lastd+1):length(HRd)])    #deaths off treatment - not for use
    if(all(is.na(Phat))){
      ## make guess for P
      Phat <- Ihat; sEP <- 2*sEI
      if(verbose) cat('No Phat supplied: making a guess from Ihat!\n')
    }
    didx <- 1:lastd #data range
    if(verbose) cat('...nahead=',nahead,'\n')
    if(verbose) cat('...lastd=',lastd,'\n')
    ANS <- Cprojections(year=year[didx],
                        Ihat=Ihat[didx],sEI=sEI[didx],
                        Nhat=Nhat[didx],sEN=sEN[didx],
                        Mhat=Mhat[didx],sEM=sEM[didx],
                        Phat=Phat[didx],sEP=sEP[didx],
                        nahead=nahead,
                        logIRR=logIRR,
                        logIRRdelta=logIRRdelta,
                        logORpsi=logORpsi,
                        returntype=output,
                        modeltype=modeltype,
                        verbose=verbose,
                        ...
                        )
    if(verbose) cat('...Cprojections returned OK...\n')
  }
  if(modeltype!='failsafe' & !returninternalfit){
    ANS <- data.table::dcast(ANS[variable %in% c('Incidence','Notifications',
                                                 'Prevalence','Deaths','time')],
                             time ~ variable,value.var=c('mid','lo','hi'))
    ## safe renaming
    n <- names(ANS)
    names(ANS)[n=='time'] <- 'year'
    names(ANS)[n=='mid_Incidence'] <- 'I.mid'
    names(ANS)[n=='lo_Incidence'] <- 'I.lo'
    names(ANS)[n=='hi_Incidence'] <- 'I.hi'
    names(ANS)[n=='mid_Notifications'] <- 'N.mid'
    names(ANS)[n=='lo_Notifications'] <- 'N.lo'
    names(ANS)[n=='hi_Notifications'] <- 'N.hi'
    names(ANS)[n=='mid_Deaths'] <- 'M.mid'
    names(ANS)[n=='lo_Deaths'] <- 'M.lo'
    names(ANS)[n=='hi_Deaths'] <- 'M.hi'
    names(ANS)[n=='mid_Prevalence'] <- 'P.mid'
    names(ANS)[n=='lo_Prevalence'] <- 'P.lo'
    names(ANS)[n=='hi_Prevalence'] <- 'P.hi'
    ANS$year <- year[ANS$year]
    ANS <- ANS[!is.na(year)] #removes 1 ahead if 'fit'
    ANS[,c('I.sd','N.sd','M.sd','P.sd'):=
           list((I.hi-I.lo)/3.92,(N.hi-N.lo)/3.92,(N.hi-N.lo)/3.92,(N.hi-N.lo)/3.92)]
    setcolorder(ANS,neworder = c("year",
                                 "I.mid", "I.sd",  "I.lo", "I.hi",
                                 "N.mid", "N.sd",  "N.lo", "N.hi", 
                                 "M.mid", "M.sd",  "M.lo", "M.hi",
                                 "P.mid", "P.sd",  "P.lo", "P.hi" ))
  }
  ## output
  ANS
}




##' SSM-based Projections
##'
##' This wraps the SSM-based models and projections using BSSM.
##' @title Cprojections
##' @param year vector of years for which there is data or an estimate is required
##' @param Ihat Incidence midpoints (NA if projection/imputation needed)
##' @param sEI Incidence uncertainty as SD (NA if projection/imputation needed)
##' @param Nhat Notifications midpoints (NA if projection/imputation needed)
##' @param sEN Notifications uncertainty as SD (NA if projection/imputation needed)
##' @param Mhat Deaths midpoints (NA if projection/imputation needed)
##' @param sEM Deaths uncertainty as SD (NA if projection/imputation needed)
##' @param Phat Prevalence midpoints (NA if projection/imputation needed)
##' @param ... Additional parameters to help tune away from defaults
##' @param sEP Prevalence uncertainty as SD (NA if projection/imputation needed)
##' @param returntype Determines what is returned if projecting
##' @param nahead how many years to project
##' * `futureonly': (default) only returns projection
##' * `fit': returns SSM output for all times
##' * `projection`: returns inputs during data combined with projections after
##' @param modeltype String to specify which SSM variant to use:
##'  * `rwI': (default) random walk for incidence. NOTE no indirect impact
##'  * `IP': AR(1) on a mixture of log(Incidence) and log(Prevalence).
##' @param verbose (Default=FALSE) Give more output for use in debugging.
##' @return A data.frame/data.table with the projections and uncertainty
##' @examples
##' 
##' 1+1 #TODO
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
                        ...,
                        nahead=0,
                        returntype='projection',
                        modeltype='rwI',
                        verbose=FALSE
                        ){

  ## avoiding warnings:
  override <- logIRR <- logIRRdelta <- variable <- value <- time <- NULL
  ## set up
  if(nahead==0 & returntype=='futureonly') stop("Can't have nahead=0 and only return future!")
  if(verbose) cat('Using Cprojections...\n')
  arguments <- list(...)
  list2env(arguments,envir = environment())                   #boost ... to this scope

  cat('...nahead = ',nahead,'...\n')

  ## switch to fit if nahead==0
  if(nahead==0 & returntype=='projection'){
    returntype <- 'fit'
    cat('...NB changing returntype to fit since nahead==0...\n')
  }

  ## model choice safety
  if(! (modeltype=='rwI' | substr(modeltype,1,2)=='IP') ){
    wrn <- paste0('modeltype = ',modeltype,' is unknown! Using rwI\n')
    warning(wrn)
    modeltype <- 'rwI'
  }

  ## processing data
  Yhat <- cbind(Ihat,Phat,Nhat,Mhat)
  NoverI <- Nhat[1]/Ihat[1]

  ## transformations NOTE reconsider
  Yhat <- log(Yhat)
  Vhat <- cbind( rep(1,nrow(Yhat)), #I
                rep(1/3,nrow(Yhat)),  #P
                rep(1/10,nrow(Yhat)),  #N
                rep(1/3,nrow(Yhat)) ) #D
  if('override' %in% names(arguments)){
    if('Vhat' %in% names(override)){
      if(override[['Vhat']]=='literal'){
        Vhat <- cbind(sEI,sEP,sEN,sEM) #NOTE these are in logspace
        if(verbose) cat('...** overriding Vhat **...\n')
      }
    }
  }

  ## Initial states:
  IS <- c(log(c(Ihat[1],Phat[1],Nhat[1],Mhat[1])),
                 c(0.1,0.1,0.1,0.1))
  if('override' %in% names(arguments)){
    if('IS' %in% names(override)){
      IS <- override[['IS']]
      if(verbose) cat('...** overriding IS **...\n')
    }
  }
  names(IS) <- c("mI0","mP0","mN0","mD0","sI0","sP0","sN0","sD0")
  if(verbose) cat('...initial state: IS = \n')
  if(verbose) print(IS)

  ## other prior parameters
  sdelta0 <- 0.5 #NOTE for rwI only prior width for detection rate
  momega0 <- -1.1
  mpsi0 <- logit(0.5)## log(tmp$Mhat[1]/tmp$Ihat[1]) + 1.1, mortality lit
  somega0 <- 0.7
  spsi0 <- 0.3

  ## initial thetas
  initial_theta <- c(logSI = -1,logsdelta=-1,logsomega=-1) #default rwI
  if('override' %in% names(arguments)){
    if('initial_theta' %in% names(override)){
      initial_theta <- override[['initial_theta']]
      if(verbose) cat('...** overriding initial_theta **...\n')
    }
  }
  initial_theta_ip <- c(logSI = -1,logsdelta=-1,logsomega=-1,logphiP=-1,logitpr=-1) #IP
  if('override' %in% names(arguments)){
    if('initial_theta_ip' %in% names(override)){
      initial_theta_ip <- override[['initial_theta_ip']]
      if(verbose) cat('...** overriding initial_theta_ip **...\n')
    }
  }
  known_params <- c(mI0 = (IS['mI0']),
                    mP0 = (IS['mP0']),
                    mN0 = (IS['mN0']),
                    mD0 = (IS['mD0']),
                    momega = momega0,
                    mdelta = log(NoverI), #N/I as proxy for N/P
                    mpsi = mpsi0,## log(tmp$Mhat[1]/tmp$Ihat[1]) + 1.1, mortality lit
                    sI0 = (IS['sI0']), #NOTE
                    sP0 = (IS['sP0']),
                    sN0 = (IS['sN0']),
                    sD0 = (IS['sD0']),
                    somega = somega0,
                    sdelta = sdelta0,
                    spsi = spsi0)
  if('override' %in% names(arguments)){
    if('nathist' %in% names(override)){
      for(nm in names(override[['nathist']])){ #loop over these
        if(verbose) cat('...** overriding ',nm,' **...\n')
        if(!nm %in% names(known_params)) stop('nathist override not found in known_params: probably not the behaviour you were looking for!')
        known_params[nm] <- override[['nathist']][[nm]]
      }
    }
  }
  known_tv_params <- cbind(Vhat,matrix(0,nrow=nrow(Vhat),ncol=2))

  ## for predictions
  future_known_tv_params <- known_tv_params[rep(nrow(Vhat),nahead),]
  ## NOTE interventions
  future_known_tv_params[,5] <- logIRR
  future_known_tv_params[,6] <- logIRRdelta

  tt3SNMZ <- c('logIncidence','logPrevalence','logNotifications','logDeaths',
               'logomega','logdelta','psi')


  ## tests:
  state <- c(log(100), log(100), log(100), log(10),
             log(3),log(1),logit(0.5))

  if(verbose) cat('Creating models...\n')
  ## --- create models
  ## model pointers
  pntrsrw <- create_xptrs() #create pointers for rwI model
  pntrsip <- create_xptrs_ip_all() #create pointers for IP model

  if(verbose){
    cat('Testing pointers:\n')
    cat('rwI...\n')
    (ta1 <- a1_fn(initial_theta,known_params))
    print(c(ta1))
    (tP1 <- P1_fn(initial_theta,known_params))
    print(c(diag(tP1)))
    tmp <- rep(0,7)
    for(i in 1:7) tmp[i] <- exparz(ta1[i],sqrt(tP1[i]))$mn
    print(tmp)
    (tH <- H_fn(1,state,initial_theta,known_params,known_tv_params))
    (tR <- R_fn(1,state,initial_theta,known_params,known_tv_params))
    (tZ <- Z_fn(1,state,initial_theta,known_params,known_tv_params))
    (tdZ <- Z_gn(1,state,initial_theta,known_params,known_tv_params))
    (tT <- T_fn(1,state,initial_theta,known_params,known_tv_params))
    (tdT <- T_gn(1,state,initial_theta,known_params,known_tv_params))
    (log_prior_pdf(initial_theta))
    cat('IP...\n')
    print(names(pntrsip))
    (ta1 <- a1_fn_ip(initial_theta_ip,known_params))
    (tP1 <- P1_fn_ip(initial_theta_ip,known_params))
    (tH <- H_fn_ip(1,state,initial_theta_ip,known_params,known_tv_params))
    (tR <- R_fn_ip(1,state,initial_theta_ip,known_params,known_tv_params))
    (tZ <- Z_fn_ip(1,state,initial_theta_ip,known_params,known_tv_params))
    (tdZ <- Z_gn_ip(1,state,initial_theta_ip,known_params,known_tv_params))
    (tT <- T_fn_ip(1,state,initial_theta_ip,known_params,known_tv_params))
    (tdT <- T_gn_ip(1,state,initial_theta_ip,known_params,known_tv_params))
    cat('--- IP prior variants ---\n')
    print(log_prior_pdf_ip4(initial_theta_ip))
    print(log_prior_pdf_ip3(initial_theta_ip))
    print(log_prior_pdf_ip2(initial_theta_ip))
    print(log_prior_pdf_ip1(initial_theta_ip))
    print(log_prior_pdf_ip0(initial_theta_ip))
    print(log_prior_pdf_ipn4(initial_theta_ip))
    print(log_prior_pdf_ipn3(initial_theta_ip))
    print(log_prior_pdf_ipn2(initial_theta_ip))
    print(log_prior_pdf_ipn1(initial_theta_ip))
    cat('...done.\n')
  }

  ## rwI
  modelrwi <- bssm::ssm_nlg(y = Yhat,
                            a1=pntrsrw$a1, P1 = pntrsrw$P1, 
                            Z = pntrsrw$Z_fn, H = pntrsrw$H_fn, T = pntrsrw$T_fn, R = pntrsrw$R_fn, 
                            Z_gn = pntrsrw$Z_gn, T_gn = pntrsrw$T_gn,
                            theta = initial_theta, log_prior_pdf = pntrsrw$log_prior_pdf,
                            known_params = known_params, known_tv_params = known_tv_params,
                            n_states = 7, n_etas = 4,
                            state_names = tt3SNMZ)


  ## IP
  logIPprior <- pntrsip$log_prior_pdf_ip0 #safety for rwI
  if(modeltype=='IP1'){
    logIPprior <- pntrsip$log_prior_pdf_ip1
  } else if(modeltype=='IP2'){
    logIPprior <- pntrsip$log_prior_pdf_ip2
  } else if(modeltype=='IP3'){
    logIPprior <- pntrsip$log_prior_pdf_ip3
  } else if(modeltype=='IP4'){
    logIPprior <- pntrsip$log_prior_pdf_ip4
  } else if(modeltype=='IPn1'){
    logIPprior <- pntrsip$log_prior_pdf_ipn1
  } else if(modeltype=='IPn2'){
    logIPprior <- pntrsip$log_prior_pdf_ipn2
  } else if(modeltype=='IPn3'){
    logIPprior <- pntrsip$log_prior_pdf_ipn2
  } else if(modeltype=='IPn4'){
    logIPprior <- pntrsip$log_prior_pdf_ipn4
  } else if(modeltype=='IP0' | modeltype=='IP'){
    logIPprior <- pntrsip$log_prior_pdf_ip0
  }

  modelip <- bssm::ssm_nlg(y = Yhat,
                           a1=pntrsip$a1, P1 = pntrsip$P1,
                           Z = pntrsip$Z_fn_ip, H = pntrsip$H_fn_ip,
                           T = pntrsip$T_fn_ip, R = pntrsip$R_fn_ip,
                           Z_gn = pntrsip$Z_gn_ip, T_gn = pntrsip$T_gn_ip,
                           theta = initial_theta_ip, log_prior_pdf = logIPprior,
                           known_params = known_params, known_tv_params = known_tv_params,
                           n_states = 7, n_etas = 4,
                           state_names = tt3SNMZ)

  ## choose model to use

  if(modeltype=='rwI'){
    cat('Using model type: ', modeltype,'\n')
    model <- modelrwi
  } else if(substr(modeltype,1,2)=='IP'){
    cat('Using model type: ', modeltype,'\n')
    model <- modelip
  }

  if(verbose) cat('Using return type: ',returntype,'\n')
  if(verbose) cat('Starting inference...\n')

  ## inference
  mcmc.type <- 'ekf'              #change inference type
  ITER <- 6000 ; BURN <- 1000 
  mcmc_fit <- bssm::run_mcmc(model, iter = ITER, burnin = BURN,mcmc_type = "ekf")

  if(verbose) cat('Postprocessing inference...\n')
  cat('calculating fit summary...\n')
  outsf <- mcmcsmry(mcmc_fit) #summarizer to mid/lo/hi

  ## predict
  if(nahead>1){
    if(verbose) cat('Making predictions...\n')
    if(verbose) cat('nahead: ',nahead,' > 1 ...\n')
    future_model <- model
    future_model$y <- ts(matrix(NA, ncol=4,nrow=nahead),
                         start = tsp(model$y)[2] + deltat(model$y),
                         frequency = frequency(model$y))
    future_model$known_tv_params <- future_known_tv_params
    pred <- predict(mcmc_fit, model = future_model, type = "state", 
                    nsim = 1000)
    mcmc_fit <- pred

    if(verbose) cat('Postprocessing projection results...\n')
    outs <- mcmcsmry(mcmc_fit) #summarizer to mid/lo/hi
  }

  if(returntype=='futureonly'){
    ## No action needed
    cat('future only, no summary...\n')
    return(outs)
  }
  if(returntype=='fit'){
    cat('returning fit summary...\n')
    return(outsf)
    ## outs <- rbind(outsf,outs) #combine with past fit
  }
  if(returntype=='projectionfit'){
    cat('returning fit summary...\n')
    outsb <- rbind(outsf[time<max(time)],outs) #combine with past fit
    return(outsb)
  }
  if(returntype=='projection'){
    cat('projection fit summary...\n')
    ## create same format input data
    inputs.m <- data.table::data.table(Incidence=Ihat,Notifications=Nhat,
                                       Deaths=Mhat,Prevalence=Phat,
                                       time=1:length(year))
    inputs.h <- data.table::data.table(Incidence=Ihat+1.96*sEI,
                                       Notifications=Nhat+1.96*sEN,
                                       Deaths=Mhat+1.96*sEM,
                                       Prevalence=Phat+1.96*sEP,
                                       time=1:length(year))
    inputs.l <- data.table::data.table(Incidence=pmax(0,Ihat-1.96*sEI),
                                       Notifications=pmax(0,Nhat-1.96*sEN),
                                       Deaths=pmax(0,Mhat-1.96*sEM),
                                       Prevalence=pmax(0,Phat-1.96*sEP),
                                       time=1:length(year))
    inputs.m <- melt(inputs.m,id.vars = c('time'))
    names(inputs.m)[3] <- 'mid'
    inputs.l <- melt(inputs.l,id.vars = c('time'))
    names(inputs.l)[3] <- 'lo'
    inputs.h <- melt(inputs.h,id.vars = c('time'))
    names(inputs.h)[3] <- 'hi'
    inputs.a <- merge(inputs.m,inputs.l,id=c('time'))
    inputs.a <- merge(inputs.a,inputs.h,id=c('time'))
    ## ouput
    outs <- rbind(inputs.a,outs)
    return(outs)
  }
}
