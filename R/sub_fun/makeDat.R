#' makeDat set up data list and map for futR
#'
#' This function prepares the ceattle estimated data for futR
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @weblink 
#' @param rec_years   required 'years' that match SSB and Rec vectors
#' @param SSB         required vector of SSB (usually previous yr if Rec = age 1)
#' @param Rec         required vector of Rec (usually previous yr if Rec = age 1)
#' @param sdRec       optional sd of rec vector 
#' @param sdSSB       optional sd of SSB vector a
#' @param sigMethod   if set to 2 it will not estimate logsigma
#' @param tMethod     apply 0-1 transformation to gamma for ricker model; if 1  = complementary loglog, if 2 = logit 
#' @param typeIN      default is 3; Recruitment model formulation based Climate enhanced Rec~SSB functions from Holsman et al. 2019; 1 = Linear with biomass (y-1), 2 =  linear, 3 = Beverton holt, 4 = Ricker
#' @param estparams   TRUE / FALSE for if parm will be estimated (default is init vals) c(log_a = TRUE, log_b=TRUE, rs_parm=FALSE, logsigma=TRUE)
#' @param covars      default is NULL; data.frame of covariates (rows) of z-score scaled environmental covariates to evaluate in Rec for each year (columns). Must match  Rec and SSB row index
#' @param covars_sd   default is NULL; data.frame of sd for covariates (rows) of z-score scaled environmental covariates to evaluate in Rec for each year (columns). Must match Rec and SSB row index
#' @param startVal    default is NULL; starting values for parameters can be a subset list but names must match
#' @param phases      default is 1 for all parms except sigma (phase 2), but can be specified.
#' @param fityrs      default is NULL; which years to fit (can be all or subset for out of sample analyses)
#' @param REcompile   default is TRUE; if true will recompile futR.cpp in TMB TRUE/FALSE
#' @export            list with :
#' *  parameters      list of parameters and starting vlaues
#' *  rs_dat          data.frame data.frame of 'years', 'SSB', 'R_obs' and optionally, 'sdSSB','sdR_obs'
#' *  maplist         list with map values for calling tmb
#' *  estparams       vecctor of TRUE/ FALSE for each parameter 
#' *  phases          vector of phases for each parameter
#' *  inputs          List of inputs to makeDat()
#' @examples
#'  
makeDat <- function(
  rec_years,
  Rec,
  SSB,
  sdSSB = NULL,
  sdRec = NULL,
  typeIN,
  tauIN =  1,
  gammaIN = NULL,
  estparams  = c(
    log_a        = TRUE, 
    log_b        = TRUE, 
    #logit_tau    = FALSE,
    gamma        = FALSE,
    beta         = FALSE,
    lambda       = TRUE,
    epsi_s       = FALSE,
    logsigma     = TRUE),
  covars    = NULL,
  covars_sd = NULL,
  startVal  = NULL,
  phases    = NULL,
  fityrs    = NULL ,
  REcompile = TRUE,
  Eat_cov   = NULL,
  tMethod   = 1,
  sigMethod = NULL) {

    inputs<-list(
      dataIN     = data.frame(rec_years,SSB,Rec,sdSSB,sdRec), 
      typeIN     = typeIN,
      tauIN      = tauIN,
      estparams  = estparams,
      covars     = covars,
      covars_sd  = covars_sd,
      startVal   = startVal,
      phases     = phases,
      fityrs     = fityrs ,
      #rationIN   = rationIN,
      REcompile  = REcompile,
      #Eat_cov    = Eat_cov,
      tMethod    = tMethod,
      sigMethod  = sigMethod)

    allyrs       <-  rec_years
    # if(is.null(rationIN))
    #   rationIN<-rep(0,length(allyrs))

    ix<-1:length(allyrs)
    if(!is.null(fityrs))  
      ix           <-  which(allyrs%in%fityrs)
  
    if(is.null(typeIN)) typeIN <-3
    
    rs_dat<-list()
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
    rs_dat$tMethod         <- 
    rs_dat$tau             <-  NA
    
    rs_dat$rectype         <-  typeIN
    rs_dat$tMethod         <-  tMethod
    rs_dat$years           <-  allyrs[ix]
    rs_dat$S_obs           <-  SSB[ix]  #Rec[y]~  S_obs[y-1]
    rs_dat$R_obs           <-  Rec[ix]  #Rec[y]~  S_obs[y-1]
    #rs_dat$Ration_scaled   <-  rationIN[ix]
    if(!is.null(covars)) 
      rs_dat$rs_cov        <-  covars[,ix]
    rs_dat$sdrs_cov        <-  rs_dat$rs_cov*0
    if(!is.null(covars_sd)) 
      rs_dat$sdrs_cov      <-  covars_sd[,ix]
    rs_dat$tau             <-  tauIN
    
    if(is.null(sdSSB))  { rs_dat$sdS <- rep(0,length(ix))}else{rs_dat$sdS <- sdSSB}
    if(is.null(sdRec))  { rs_dat$sdR <- rep(0,length(ix))}else{rs_dat$sdR <- sdRec}
    if(is.null(covars)){
      rs_dat$rs_cov    <-   t(data.frame(NOcov=rep(0,length(ix))))
      rs_dat$ncov      <-   NULL
      estparams["beta"] <- FALSE
      estparams["lambda"] <- FALSE
    }else{
      rs_dat$ncov      <-   dim(covars)[1]
      if(rs_dat$ncov==1){
        rs_dat$rs_cov    <-    t(as.matrix(rs_dat$rs_cov))
        rs_dat$sdrs_cov  <-    t(as.matrix(rs_dat$sdrs_cov))
      }
    }
    
    gamma <- -2
    if(!is.null(gammaIN)) gamma <- gammaIN
    
    rs_dat$cov_type  <-    rep(0,length(ix))
    rs_dat$tMethod           <-   tMethod
    rs_dat$cov_type          <-    rep(0,length(ix))
    rs_dat$tMethod           <-    tMethod
    
    if(sigMethod!=2&estparams['logsigma']){
      rs_dat$sigMethod           <-   sigMethod
      estparams['logsigma']      <-   TRUE
    }else{
      rs_dat$sigMethod           <-   sigMethod
      estparams['logsigma']      <-   FALSE
    }
    
    
    # if(!is.null(Eat_cov))
    #   rs_dat$cov_type[Eat_cov]   <-   1
    
    rs_dat$nyrs            <-  dim(rs_dat$rs_cov)[2]
    rs_dat$ncov            <-  dim(rs_dat$rs_cov)[1]

  #___________________________________________
  # 4.1 set initial conditions:
  
    ln_mn_rec     <-  mean(log(Rec[ix]))
    ln_mn_SSB     <-  mean(log(SSB[ix]))
    

    if(typeIN==1){
      #Linear
        tmp_a       =  .8*ln_mn_rec
        tmp_b       =  1.0/log(ln_mn_SSB)
        tmpg        =  0
        aphase      =  1
        bphase      = -4
        gphase      = -4
        sigma_phase =  2
     }

    
    if(typeIN==2){
      #BH
        tmp_a       =   1.2
        tmp_b       =  -ln_mn_rec
        tmpg        =   0
        bphase      =   1
        aphase      =   1
        gphase      =  -4
        sigma_phase =   2
    }  
    
    if(typeIN==3){
      #Ricker
      tmp_a       =   1.2
      tmp_b       =  -ln_mn_rec
      tmpg        =   gamma
      # if(tMethod ==1)
      #   tmpg  = 1-exp(-exp(gamma)) 
      # if(tMethod ==2)
      #   tmpg  = exp(gamma)/(1+exp(gamma))
      bphase      =   1
      aphase      =   1
      gphase      =   2
      sigma_phase =   2
    }
    
    if(typeIN==4){
      tmp_a       =  .8*ln_mn_rec
      tmp_b       =   -ln_mn_rec
      tmpg        =   1
      bphase      =   1
      aphase      =   1
      gphase      =  -4
      sigma_phase =   2
    }

    if(is.null(phases)){
      phases      <-   c(
        log_a        = aphase, 
        log_b        = bphase ,
        gamma       = gphase,
        #logit_tau    = sigma_phase,
        beta         = aphase+1,
        lambda       = aphase+1,
        epsi_s       = aphase+2,
        logsigma     = sigma_phase)
    }
     
     parameters      <-   list(
      log_a        = tmp_a, 
      log_b        = tmp_b ,
      #logit_tau    = -10,
      gamma        = tmpg ,
      beta         = rep(0,rs_dat$ncov),
      lambda       = rep(0,rs_dat$ncov),
      #epsi_s       = rep(0,rs_dat$nyrs),
      epsi_s       = 0,
      logsigma     = log(.9))


   if(!is.null(startVal)){
    if(!any(names(startVal)%in%names(parameters))) stop("startVal names do not match parameters")
      tmpc<-which(names(startVal)%in%names(parameters))
      for(ii in 1:length(tmpc))
        parameters[[ names(startVal)[ii] ]]<-startVal[[ii]]
   }
    maplist<-makeMap(param=parameters,estpar=estparams)
    return(list(parameters = parameters, rs_dat = rs_dat, maplist= maplist,estparams = estparams, phases=phases,inputs=inputs))
}









