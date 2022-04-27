#' run the futR() recruitment model
#'
#' runmod() will run the futR() recruitment model
#' @import futR
#' @email For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @weblink 
#' @param dlistIN     list of input data and parameters used for the makeADfun() dependency
#' * dlistIN$parameters  list of parameters for the TMB futR model:
#'    [] dlistIN$parameters$log_a    scalar starting value for log_a (intercept)
#'    [] dlistIN$parameters$log_b    scalar starting value for log_b (slope) 
#'    [] dlistIN$parameters$beta     vector of starting values for lambda (env. effects on pre-spawning)   
#'    [] dlistIN$parameters$lambda   vector of starting values for beta (env. effects on post-spawning)   
#'    [] dlistIN$parameters$epsi_s   starting value for epsi_s (error around S estimates) 
#'    [] dlistIN$parameters$logsigma starting value for logsigma (process error)
#'    [] dlistIN$parameters$skipFit  starting value for skipFit (switch used for projecting)
#' * dlistIN$rs_dat list of data to read into the model
#'    --- Observation Error inputs
#'    [] tau        scalar for observation and process error
#'    [] sigMethod  Integer (switch): which method for estimating sigma
#'    --- Covariate inputs
#'    [] ncov       Integer number of environmental covariates
#'    [] cov_type   (defunct) vector of covariate type if 1 then run ration
#'    [] nyrs       Integer number of years in the data
#'    [] lambda_0   matrix (ncov, nyrs), binary (0,1) matrix of covar. inclusion in post-spawning effects
#'    [] beta_0     matrix (ncov, nyrs), binary (0,1) matrix of covar. inclusion in pre-spawning effects
#'    [] sdrs_cov   matrix (ncov, nyrs) of the stdev of environmental covariates
#'    [] rs_cov     matrix (ncov, nyrs) of the of environmental covariates
#'    --- RS data inputs
#'    [] sdR        vector (nyrs) of stdev of recruitment in a given year
#'    [] sdS        vector (nyrs) of stdev of spawnerindex in a given year
#'    [] R_obs      vector (nyrs) of recruitment in a given year
#'    [] S_obs      vector (nyrs) of spawner index in a given year
#'    [] years      vector (nyrs) recruiment year
#'    --- Switches
#'    [] rectype    Integer (switch) which type of recruitment model: 1 = loglinear, 2 = loglinear ~f(S), 3 = Beverton-holt, 4 = Ricker
#'    [] estMode    Integer (switch) should the model run in estimation mode (y=1,no=0)?
#' * dlistIN$maplist   list of the map for the TMB model
#' * dlistIN$estparams vector of T/F of which parameters will be estimated
#' * dlistIN$phases    (defunct) phases for estimating parameters
#' * dlistIN$inputs    archive of intputs for this function
#' @param version     recruitment model to run. Default is 'futR'
#' @param src_fldr    subfolder where model is located. Default is 'src'
#' @param silentIN    T/F , default = 'TRUE'
#' @param se.fit      T/F , default = 'TRUE'
#' @param simulate    T/F , simulate the model? default = 'FALSE'
#' @param sim_nitr    number of iterations for the simulate function
#' @param recompile   recompile the .cpp file? default is 'FALSE'
#' @param maxitr      10000  Max iterations for fitting the objective function
#' @param maxeval     10000  Max evaluations for fitting the objective function
#' @return returns  summary of the model including the mle
#' 
#' @examples
#' datlist <- readMake_futR_data("data/in/futR_Inputs.xlsx" )
#' mm      <- runmod(dlistIN   = datlist, version   = 'futR',recompile = FALSE,simulate  = TRUE,sim_nitr  = 1000)  
#' @export
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