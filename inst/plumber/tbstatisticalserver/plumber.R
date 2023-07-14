#* Echo the parameter that was sent in (for testing)
#* @param msg The message to echo back.
#* @get /echo
function(msg=""){
  list(msg = paste0("The message is: '", msg, "'"))
}




#' @filter cors
cors <- function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")

  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods","*")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)

    res$status <- 200

    return(list())
  } else {
    plumber::forward()
  }
}

#* @get /version
#* @serializer unboxedJSON
function() {
  list(title="Fort", version = "0.1.7")
}

#* @post /projection
function(
    year,
    iHat = NA,
    sEi = NA,
    nHat = NA,
    sEn = NA,
    mHat = NA,
    sEm = NA,
    pHat = NA,
    sEp = NA,
    tXf = 0,
    hRd = 1,
    hRi = 1,
    oRt = 1,
    nReplicates = 500,
    output="projection",
    modelType="failsafe") {
  set.seed(173)

  data <- fort::projections(
    year = year,
    Ihat = iHat,
    sEI = sEi,
    Nhat = nHat,
    sEN = sEn,
    Mhat = mHat,
    sEM = sEm,
    Phat = pHat,
    sEP = sEp,
    TXf = tXf,
    HRd = hRd,
    HRi = hRi,
    ORt = oRt,
    nrep = nReplicates,
    output = output,
    modeltype = modelType
  ) |> 
    lapply(unlist)

  list(
    year = data$year,
    iMid =data$I.mid,
    iSd =data$I.sd,
    iLo =data$I.lo,
    iHi =data$I.hi,
    nMid =data$N.mid,
    nSd =data$N.sd,
    nLo =data$N.lo,
    nHi =data$N.hi,
    mMid =data$M.mid,
    mSd =data$M.sd,
    mLo =data$M.lo,
    mHi =data$M.hi,
    pMid =data$P.mid,
    pSd =data$P.sd,
    pLo =data$P.lo,
    pHi =data$P.hi
  )
}
