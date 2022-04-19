#-------------------------------------
# run the model
#-------------------------------------
#' runmod will run the futR.cpp recruitment model
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' 
#' @weblink 
#' @param version     recruitment model to run. Default is futR
#' @param recompile   recompile the .cpp file? default is FALSE; only needed if .cpp file is altered
#' @param dlistIN     list of input data and parameters used for the makeADfun() dependency
#' * dlistIN$parameters  list of parameters for the TMB futR model:
#'    [] datlist$parameters$log_a  starting value for alpha parameter
#'    [] datlist$parameters$log_b  starting value for alpha parameter
#'    [] datlist$parameters$rs_parm  starting value for alpha parameter   
#'    [] datlist$parameters$epsi_s  starting value for alpha parameter
#'    [] datlist$parameters$logsigma  double with starting value ofr sigma
#' * dlistIN$rs_dat
#' * dlistIN$maplist
#' * dlistIN$estparams
#' * dlistIN$phases
#' * dlistIN$inputs
#' @export            returns tmpmod summary of the model including the mle
#' 
#' @email             For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @examples
#' 
#' mod<-estRec(dataINUSE  =  dataIN_admb)
#' 
runmod<-function(
  dlistIN    = dlist,
  version    = "futR",
  silentIN   = TRUE,
  se.fit     = TRUE,
  simulate   = FALSE,
  simnitr    = 1000,
  recompile  = FALSE,
  maxitr     = 10000,
  maxeval    = 10000){
  
  if(recompile){  compile(paste0(version, ".cpp"))  }
  dyn.load( dynlib(version) ) 
  
  
  model<- MakeADFun(
              data                 =  dlistIN$rs_dat,
              parameters           =  dlistIN$parameters,
              DLL                  =  version,
              checkParameterOrder  =  TRUE,
              hessian              =  TRUE,
              map                  =  dlistIN$maplist,
              silent               =  silentIN)
  
  tmpmod<-list()
  tmpmod$model    <-  model
  tmpmod$fit      <-  nlminb(model$env$last.par.best,model$fn,model$gr, control=list(iter.max=maxitr,eval.max=maxeval))
  tmpmod$rep      <-  sdreport(model,getJointPrecision=TRUE)
  tmpmod$out      <-  summary(tmpmod$rep)
  tmpmod$mle      <-  model$env$last.par.best
  tmpmod$Hessian  <- with(model,optimHess(tmpmod$mle,model$fn,model$gr))
  tmpmod$LnDet    <- sum(determinant(tmpmod$Hessian, logarithm=TRUE)$modulus[[1]])
  lp              <-  model$env$last.par
  
  if (!se.fit) {
    pred          <- unlist(model$report(lp))
    tmpmod$pred   <- data.frame(def= names(pred),pred=pred,pred.se=NA)
  } else {
    tmpmod$sdr    <- sdreport(model,tmpmod$mle,hessian.fixed=tmpmod$Hessian,getReportCovariance=F)
    tmpmod$sdrsum <- TMB:::summary.sdreport(tmpmod$sdr, "report")
    pred   <- tmpmod$sdrsum[,"Estimate"]
    se     <- tmpmod$sdrsum[,"Std. Error"]
    tmpmod$pred   <- data.frame(def= names(pred),pred=pred,pred.se=se)
  }
  
  if(simulate){
    sim <- replicate(simnitr, {
      simdata        <-   model$simulate(par=model$par, complete=TRUE)
      simdata$R_obs  <-   simdata$R_hat
      simdata$SSB    <-   simdata$S_hat
      obj2           <-   MakeADFun(
        data=simdata,
        parameters=dlistIN$parameters,
        DLL='futR',
        checkParameterOrder=TRUE,
        hessian = TRUE,
        map=dlistIN$maplist,
        profile = NULL, # TODO: Optionally "beta"
        silent = TRUE)
      nlminb(obj2$par, obj2$fn, obj2$gr)$par
      
    })
    tmpmod$sim <- sim
  }
  
  
  return(tmpmod)
}