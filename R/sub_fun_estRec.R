#'
#' estRec estimate recruitment with options for phasing
#'
#' This function prepares the ceattle estimated data for futR
#' For more information contact author Kirstin Holsman (kirstin.holsman at noaa.gov)
#'  
#' @param datlistIN   list of input data from makeDat() function: "parameters" "rs_dat"     "maplist"    "estparams"  "phases"  
#' @param returnAll   return all phases? TRUE/FALSE
#' @param quiet       print out phases
#' 
#' @returns a list with final (final phase) dlistfinal, dlist (all phases), phase
#' 
#' @export            
#' 
#' @examples
#' mod<-estRec(dataINUSE  =  dataIN_admb)
#' 
estRec<-function(version="futR",
                 datlistIN     =  datlist,
                 returnAll     =  FALSE, 
                 quiet         =  FALSE,
                 REcompile     =  TRUE
){
  
  parameters   <-  datlistIN$parameters 
  rs_dat       <-  datlistIN$rs_dat 
  maplist      <-  datlistIN$maplist 
  estparams    <-  datlistIN$estparams 
  phases       <-  datlistIN$phases 
  inputs       <-  datlistIN$inputs
  
  # 1. Compile model   # recall that parameters must be in the same order
  if(REcompile){
    require(TMB)  # precompile()
    compile(paste0(version, ".cpp"),"-o0 -g")
    
    source("lldbsource.r")       
    # gdbsource2('futR.r',interactive=TRUE) #- change to lldb
    # R --debugger=lldb
    # (lldb) run
    # breakpoint set --file fut_R.cpp --line 2
  }
  
  dyn.load(dynlib(version))
  
  # 2. set up phases (optional)
  estparams   <-  param_ctl$estparams
  phases      <-  param_ctl$phases
  sigMethod   <-  param_ctl$sigMethod
  phase       <-  dlist  <-  list()
  if(is.null(phases))
    phases<-rep(1,length(estparams))
  maxphase<-max(phases)
  
  # 3. loop through phases{
  for(ph in 1:maxphase){
    dlist<-makeDat(
      dataIN     =  inputs$dataIN,
      typeIN     =  inputs$typeIN,
      estparams  =  inputs$estparams,
      covarsIN   =  inputs$covarsIN,
      startVal   =  inputs$startVal,
      phases     =  inputs$phases,
      fityrs     =  inputs$fityrs,
      rationIN   =  inputs$rationIN,
      REcompile  =  inputs$REcompile,
      Eat_cov    =  inputs$Eat_cov,
      sigMethod  =  inputs$sigMethod)
    
    tmp_estparams             <-  estparams
    dlist                     <-  datlistIN
    tmp_estparams[phases>ph]  <-  FALSE
    dlist$estparams           <-  tmp_estparams
    
    
    if(!quiet) print(paste0("* phase ",ph))
    phase[[ph]]<-runmod(dlistIN=dlist,version=version,recompile=FALSE)
    inputs$startVal<-vec2list(phase[[ph]]$fit$par)
    
    #phase[[1]]<-runmod(dlistIN=dlist) 
  }
  
  if(maxphase>=2){
    if(!quiet) print("* phase 2")
    dlist[[2]]<- makeDat(
      sp        = 1,
      fityrs    = fityrsIN ,
      startVal  = vec2list(phase[[1]]$fit$par),
      dataIN    = dataINUSE, 
      rationIN  = NULL,
      typeIN    = 4, 
      REcompile = TRUE,
      Eat_cov   = Eat_covIN,
      covarsIN  = sub_cov,
      include_sdR = FALSE,
      estparams  = c(
        log_a = TRUE, 
        log_b=TRUE, 
        #logit_tau= TRUE,
        rs_parm=FALSE, 
        epsi_s = FALSE,
        logsigma=FALSE)
    )
    
    phase[[2]]<-runmod(dlistIN=dlist[[2]])
  }
  
  if(maxphase>=3){
    if(!quiet) print("* phase 3")
    dlist[[3]]<- makeDat(
      sp        = 1,
      fityrs    = fityrsIN ,
      startVal  = vec2list(phase[[2]]$fit$par),
      dataIN    = dataINUSE, 
      rationIN  = NULL, 
      typeIN    = 4, 
      REcompile = TRUE,
      Eat_cov   = Eat_covIN,
      covarsIN  = sub_cov,
      include_sdR = FALSE,
      estparams  = c(
        log_a = TRUE, 
        log_b=TRUE, 
        #logit_tau= TRUE,
        rs_parm=TRUE, 
        epsi_s = FALSE,
        logsigma=FALSE)
    )
    
    phase[[3]]<-runmod(dlistIN=dlist[[3]])
  }
  
  
  if (maxphase>=4){
    # FIX not yet working right
    if(!quiet) print("* phase 4")
    dlist[[4]]<- makeDat(
      sp        = 1,
      fityrs    = fityrsIN ,
      dataIN    = dataINUSE, 
      rationIN  = NULL,
      startVal  = vec2list(phase[[3]]$fit$par),
      typeIN    = 4, 
      REcompile = TRUE,
      Eat_cov   = Eat_covIN,
      covarsIN  = sub_cov,
      include_sdR = TRUE,
      estparams  = c(
        log_a = TRUE, 
        log_b=TRUE, 
        #logit_tau= TRUE,
        rs_parm=TRUE, 
        epsi_s = TRUE,
        logsigma=FALSE)
    )
    
    phase[[4]]<-runmod(dlistIN= dlist[[4]])
  }
  if (maxphase>=5){
    if(!quiet) print("* phase 5")
    dlist[[5]]<- makeDat(
      sp        = 1,
      fityrs    = fityrsIN ,
      dataIN    = dataINUSE, 
      rationIN  = NULL,
      startVal  = vec2list(phase[[4]]$fit$par),
      typeIN    = 4, 
      REcompile = TRUE,
      Eat_cov   = Eat_covIN,
      covarsIN  = sub_cov,
      include_sdR = TRUE,
      estparams  = c(
        log_a = TRUE, 
        log_b=TRUE, 
        #logit_tau= TRUE,
        rs_parm=TRUE, 
        sdS=FALSE,
        sdR=TRUE,
        logsigma=TRUE)
    )
    
    phase[[5]]<-runmod(dlistIN=dlist4)
  }
  if(returnAll)
    return(list(final=phase[[maxphase]],dlistfinal=dlist[[maxphase]],phase=phase,dlist=dlist))
  else
    return(list(final=phase[[maxphase]],dlistfinal=dlist[[maxphase]] ))
}

