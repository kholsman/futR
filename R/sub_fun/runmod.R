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
#'    [] datlist$parameters$sig_sobs  starting value for alpha parameter
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
  dlistIN,
  src_fldr   = "src",
  version    = "futR",
  silentIN   = TRUE,
  se.fit     = TRUE,
  simulate   = FALSE,
  sim_nitr   = 1000,
  recompile  = FALSE,
  maxitr     = 10000,
  maxeval    = 10000){
  
  wd0 <- getwd()
  if(!is.null(src_fldr))
    setwd(src_fldr)
  
  if(recompile){  
    if(file.exists(paste0(version,".o")))
       file.remove(paste0(version,".o"))
    if(file.exists(paste0(version,".so")))
      file.remove(paste0(version,".so"))
    compile(paste0(version, ".cpp")) 
  }
  
  dyn.load( dynlib(version) ) 

  model  <- MakeADFun(
              data                 =  dlistIN$rs_dat,
              parameters           =  dlistIN$parameters,
              DLL                  =  version,
              checkParameterOrder  =  TRUE,
              hessian              =  TRUE,
              map                  =  dlistIN$maplist,
              silent               =  silentIN)
  
  tmpmod  <- list()
  tmpmod$model    <-  model
  tmpmod$fit      <-  nlminb(model$env$last.par.best,
                             model$fn,
                             model$gr, 
                             control=list(iter.max=maxitr,eval.max=maxeval))
  if(tmpmod$fit$objective==Inf) 
    stop("Problem with objective function (Inf)")
  if(is.na(tmpmod$fit$objective)) 
    stop("Problem with objective function (NaN)")
  tmpmod$objFun   <-  model$fn()
  tmpmod$report   <-  model$report()
  tmpmod$sdreport <-  sdreport(model,getJointPrecision=TRUE)
  tmpmod$rep      <-  sdreport(model,getJointPrecision=TRUE)
  tmpmod$input    <-  dlistIN$rs_dat
  tmpmod$out      <-  summary(tmpmod$rep)
  tmpmod$mle      <-  model$env$last.par.best
  tmpmod$Hessian  <-  with(model,optimHess(tmpmod$mle,model$fn,model$gr))
  tmpmod$LnDet    <-  sum(determinant(tmpmod$Hessian, logarithm=TRUE)$modulus[[1]])
  lp              <-  model$env$last.par
  
  if (!se.fit) {
    pred          <- unlist(model$report(lp))
    tmpmod$pred   <- data.frame(def= names(pred),pred=pred,pred.se=NA)
  } else {
    tmpmod$sdr    <- sdreport(model,tmpmod$mle,hessian.fixed=tmpmod$Hessian,getReportCovariance=F)
    #tmpmod$sdrsum <- TMB:::summary.sdreport(tmpmod$sdr, "report")
    tmpmod$sdrsum <- summary(sdreport(model,getJointPrecision=TRUE))
    pred          <- tmpmod$sdrsum[,"Estimate"]
    se            <- tmpmod$sdrsum[,"Std. Error"]
    tmpmod$pred   <- data.frame(def= names(pred),pred=pred,pred.se=se)
  }
  
  if(simulate){
    #tmpmod$simdat <- array(sim_nitr)
    sim <- replicate(sim_nitr, {
      simdata <- model$simulate(par = model$par, complete=TRUE)
      
      obj2    <- MakeADFun(data       = simdata,
                         parameters = dlistIN$parameters,
                         DLL        = version,
                         checkParameterOrder=TRUE,
                         hessian    = TRUE,
                         map        = dlistIN$maplist,
                         profile    = NULL, # TODO: Optionally "beta"
                         silent     = TRUE)
      nlminb(obj2$par, obj2$fn, obj2$gr)$par
    })
    tmpmod$sim_df <- data.frame(estimate=as.vector(sim), parameter=names(model$par)[row(sim)])
    tmpmod$sim    <- sim

  }
  setwd(wd0 )
  
  return(tmpmod)
}