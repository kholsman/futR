#' makeDat set up data list and map for futR
#'
#' This function prepares the ceattle estimated data for futR
#' For more information contact author Kirstin Holsman (kirstin.holsman at noaa.gov)
#'  
#' @param estMode     estimation mode, 0 = don't est, 1 = estimate par
#' @param rec_years   required 'years' that match SSB and Rec vectors
#' @param SSB         required vector of SSB (usually previous yr if Rec = age 1)
#' @param Rec         required vector of Rec (usually previous yr if Rec = age 1)
#' @param sdRec       optional sd of rec vector 
#' @param sdSSB       optional sd of SSB vector a
#' @param sigMethod   default = 1, if set to 2 it will not estimate logsigma
#' @param rectype     default is 3; Recruitment model formulation based Climate enhanced Rec~SSB functions from Holsman et al. 2019; 
#'  1 = Linear with biomass (y-1), 
#'  2 =  linear, 
#'  3 = Beverton holt, 
#'  4 = Ricker
#' @param estparams   vector of TRUE / FALSE for each parm,  if parm will be estimated (default is init vals) c(log_a = TRUE, log_b=TRUE, rs_parm=FALSE, logsigma=TRUE)
#' @param covars      default is NULL; data.frame of covariates (rows) of z-score scaled environmental covariates to evaluate in Rec for each year (columns). Must match  Rec and SSB row index
#' @param covars_sd   default is NULL; data.frame of sd for covariates (rows) of z-score scaled environmental covariates to evaluate in Rec for each year (columns). Must match Rec and SSB row index
#' @param beta_0    = NULL,
#' @param lambda_0  = NULL
#' @param startVal default is NULL; starting values for parameters can be a subset list but names must match
#' @param phases default is 1 for all parms except sigma (phase 2), but can be specified.
#' @param fityrs default is NULL; which years to fit (can be all or subset for out of sample analyses)
#' @param REcompile default is TRUE; if true will recompile futR.cpp in TMB TRUE/FALSE
#' @param sigMethod sigma R method
#' 
#' @export
#' 
#' @returns list with :
#' *  parameters      list of parameters and starting vlaues
#' *  rs_dat          data.frame data.frame of 'years', 'SSB', 'R_obs' and optionally, 'sdSSB','sdR_obs'
#' *  maplist         list with map values for calling tmb
#' *  estparams       vecctor of TRUE/ FALSE for each parameter 
#' *  phases          vector of phases for each parameter
#' *  inputs          List of inputs to makeDat()
#' @examples
#'  makeDat()
#'  

makeDat <- function(
  estMode = 1,  # estimate parameters
  rec_years,
  Rec,
  SSB,
  sdSSB = NULL,
  sdRec = NULL,
  rectype,
  tauIN =  0,
  gammaIN = NULL,
  estparams  = c(
    log_a        = TRUE, 
    log_b        = TRUE, 
    gamma        = FALSE,
    beta         = FALSE,
    lambda       = TRUE,
    epsi_s       = FALSE,
    logsigma     = TRUE), 
  covars    = NULL,
  covars_sd = NULL,
  beta_0    = NULL,
  lambda_0  = NULL,
  startVal  = NULL,
  phases    = NULL,
  fityrs    = NULL ,
  REcompile = TRUE,
  sigMethod = 1) {
  
  if(tauIN == 0) 
    tauIN <- 1e-3  # tau = 0 returns NAN
    datIN <- data.frame(rec_years,SSB,Rec,sdSSB,sdRec)

    inputs<-list(
      dataIN     = datIN, 
      rectype    = rectype,
      tauIN      = tauIN,
      estparams  = estparams,
      covars     = covars,
      covars_sd  = covars_sd,
      beta_0     = beta_0,
      lambda_0   = lambda_0,
      startVal   = startVal,
      phases     = phases,
      fityrs     = fityrs ,
      REcompile  = REcompile,
      sigMethod  = sigMethod)

    allyrs       <-  rec_years
    estparams <- c(estparams,skipFit= FALSE)
    if(estMode==0) estparams["skipFit"] <- FALSE
    # if(is.null(rationIN))
    #   rationIN<-rep(0,length(allyrs))

    ix<-1:length(allyrs)
    if(!is.null(fityrs))  
      ix           <-  which(allyrs%in%fityrs)
  
    if(is.null(rectype)) 
      rectype <-4
    
    rs_dat<-list()
    rs_dat$estMode         <-  
    rs_dat$rectype         <-  
    rs_dat$years           <-  
    rs_dat$S_obs           <-  
    rs_dat$R_obs           <-  
    rs_dat$sdS             <-  
    rs_dat$sdR             <-  
    rs_dat$rs_cov          <-  
    rs_dat$sdrs_cov        <-
    rs_dat$nyrs            <-
    rs_dat$cov_type        <-  
    rs_dat$ncov            <-  
    rs_dat$sigMethod       <- 
    # rs_dat$tMethod         <- 
    rs_dat$tau             <-  NA
    
    rs_dat$sigMethod       <-  sigMethod
    rs_dat$estMode         <-  estMode
    rs_dat$rectype         <-  rectype
    rs_dat$years           <-  allyrs[ix]
    rs_dat$S_obs           <-  SSB[ix]  #Rec[y]~  S_obs[y-1]
    rs_dat$R_obs           <-  Rec[ix]  #Rec[y]~  S_obs[y-1]
    #rs_dat$Ration_scaled   <-  rationIN[ix]
    if(!is.null(covars)){
      rs_dat$rs_cov        <-   covars[,ix]
      ec_use               <-   matrix(covars[,ix],dim(covars)[1],length(ix))
      rs_dat$beta_0        <-   matrix(as.numeric(beta_0)*(ec_use*0+1),dim(covars)[1],length(ix))
      rs_dat$lambda_0      <-   matrix(as.numeric(lambda_0)*(ec_use*0+1),dim(covars)[1],length(ix))
    }

    rs_dat$sdrs_cov        <-  rs_dat$rs_cov*0
    if(!is.null(covars_sd)) 
      rs_dat$sdrs_cov      <-  covars_sd[,ix]
    rs_dat$tau             <-  tauIN
    
    if(is.null(sdSSB))  { rs_dat$sdS <- rep(0,length(ix))}else{rs_dat$sdS <- sdSSB}
    if(is.null(sdRec))  { rs_dat$sdR <- rep(0,length(ix))}else{rs_dat$sdR <- sdRec}
    if(is.null(covars)){
      rs_dat$rs_cov    <- t(data.frame(NOcov=rep(0,length(ix))))
      rs_dat$sdrs_cov  <- t(data.frame(NOcov=rep(0,length(ix))))
      rs_dat$ncov      <- NULL
      rs_dat$beta_0    <- t(data.frame(NOcov=rep(0,length(ix))))
      rs_dat$lambda_0  <- t(data.frame(NOcov=rep(0,length(ix))))
      estparams["beta"]   <- FALSE
      estparams["lambda"] <- FALSE
    }else{
      rs_dat$ncov      <-   dim(covars)[1]
      if(rs_dat$ncov==1){
        rs_dat$rs_cov    <- t(as.matrix(rs_dat$rs_cov))
        #rs_dat$beta_0    <- t(as.matrix(rs_dat$beta_0))
        #rs_dat$lambda_0  <- t(as.matrix(rs_dat$lambda_0))
        rs_dat$sdrs_cov  <- t(as.matrix(rs_dat$sdrs_cov))
      }
    }
    
    gamma <- 3
    if(!is.null(gammaIN)) 
      gamma <- gammaIN
    
    rs_dat$cov_type <- rep(0,length(ix))
    
    if(sigMethod!=1){
      estparams['epsi_s'] <- TRUE
    }else{
      estparams['epsi_s'] <- FALSE
    }
    if(sigMethod!=3){
      estparams['logsigma'] <- TRUE
    }else{
      estparams['logsigma'] <- FALSE
    }
    
    rs_dat$nyrs <-  dim(rs_dat$rs_cov)[2]
    rs_dat$ncov <-  dim(rs_dat$rs_cov)[1]

  #___________________________________________
  # 4.1 set initial conditions:
  
    ln_mn_rec <-  mean(log(Rec[ix]))
    ln_mn_SSB <-  mean(log(SSB[ix]))
    
    if(rectype==1){
      #Linear
        tmp_a       =  (-.09*ln_mn_SSB)
        tmp_b       =  1.0/log(ln_mn_rec) # not used
        aphase      =  1
        bphase      =  -4
        sigma_phase =  2
        estparams["log_b"] <- FALSE
    }
    
    if(rectype==2){
      #Linear with Biomass y-1
      tmp_a       =  (-.09*ln_mn_SSB)
      tmp_b       =  1.0/log(ln_mn_rec)
      aphase      =  1
      bphase      =  1
      sigma_phase =  2
      estparams["log_b"] <- TRUE
    }

    
    if(rectype==3){
      #BH
        tmp_a       =   1.2
        tmp_b       =   -0.8*ln_mn_rec
        bphase      =   1
        aphase      =   1
        sigma_phase =   2
    }  
    
    if(rectype==4){
      #Ricker
      tmp_a       =   1.2
      tmp_b       =  -0.8*ln_mn_rec
      bphase      =   1
      aphase      =   1
      sigma_phase =   2
    }
    
    if(rectype==5){
      tmp_a       =  (.09*ln_mn_SSB)
      tmp_b       =   .4*ln_mn_rec
      bphase      =   1
      aphase      =   1
      sigma_phase =   2
    }

    if(is.null(phases)){
      phases      <-   c(
        log_a        = aphase, 
        log_b        = bphase ,
        beta         = sigma_phase+1,
        lambda       = sigma_phase+1,
        epsi_s       = sigma_phase+2,
        logsigma     = sigma_phase)
    }
     
     parameters      <-   list(
      log_a        = tmp_a, 
      log_b        = tmp_b ,
      beta         = beta_0,
      lambda       = lambda_0,
      epsi_s       = 0,
      logsigma     = log(.9),
      skipFit      = 0)

   # replace with startVal 
   if(!is.null(startVal)){
    if(!any(names(startVal)%in%names(parameters))) 
      stop("startVal names do not match parameters")
      tmpc<-which(names(startVal)%in%names(parameters))
      for(ii in 1:length(tmpc))
        parameters[[ names(startVal)[ii] ]]<-startVal[[ii]]
   }
     # print(parameters)
     # print(estparams)
    
    
    
    maplist <- makeMap(param=parameters,estpar=estparams)
    parameters$beta   <- rs_dat$beta_0[,1]*0
    parameters$lambda <- rs_dat$lambda_0[,1]*0

 
    return(list(parameters = parameters, 
                rs_dat     = rs_dat, 
                maplist    = maplist,
                estparams  = estparams, 
                phases     = phases))
}









