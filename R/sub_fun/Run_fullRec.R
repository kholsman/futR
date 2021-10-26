
#-------------------------------------
# run the model and extract the results
#-------------------------------------
#' runmod will run the futR.cpp recruitment model
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' 
#' @weblink 
#' @param h which hypothesis to run
#' 
#' @export            returns 
#' 
#' @email             For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @examples
#' Run_fullRec(s=1,h=1,tt=1,i=1)
#' 

Run_fullRec<-function(s,h,tt,i){
  
    print(paste("model ", i, "of", nmod,"hypothesis =", hypotheses[h],", rectype =", rectypes[tt],", species =", s))
    pnum      <-  which(full.mat[i,]==1)
    tmpenv    <-  t(covs[,pnum])
    nm        <-  paste0("RS",s,"_",rownames(full.mat)[i],"_",rectypes[tt])
    
    if(hypotheses[h]!="null"){
      covtxt  <-  paste(names(covs)[pnum],collapse=", ")
      txt     <-  paste0("RS",s,"_",rownames(full.mat)[i])
      
    }else{
      pnum    <-  "null"
      tmpenv  <-  NULL
      covtxt  <-  "null"
      txt     <-  paste0("RS",s,"_null")
      nm      <-  paste0("RS",s,"_null_",rectypes[tt])
    }
    covset    <-  paste0("c(",paste(pnum,collapse=", "),")")
    
    
    datlist   <-  makeDat(
      tauIN     = 1,
      sigMethod = mode, 
      tMethod   = 1,
      estparams = estparamsLIST[[h]],
      typeIN    = tt,
      dataIN    = rec,
      covarsIN  = tmpenv)
    
    # fit the model:
    futR_mods   <-  fitRec(dlistIN=datlist,silentIN = silentINN,version="futR",recompile=F,simulate=simIT)
    
    # get number of parms:
    npar        <-  length(which(!is.na(as.numeric(melt(futR_mods$model$env$map)$value))))
    
    # if(rectypes[tt]=="BLM")
    #   npar <- npar-1
    
    if(simIT){
      require(reshape)
      df           <- melt(futR_mods$sim)
      colnames(df) <- c("parameter","itr","estimate")
      df$sp        <-  s
      df$model     <-  nm
      df$modelnum  <-  i
      df$CE        <-  TRUE
      df$npar      <-  npar
      df$n         <-  dim(rec)[1]
      df$rsType    <-  rectypes[tt]
      df$hypoth    <-  hypotheses[h]
      df$covset    <-  covset
      df$covtxt    <-  covtxt
    }
    
    tmppar2   <-  futR_mods$fit$par
    tmppar2   <-  futR_mods$pred
    tmp       <-  futR_mods$pred
    
    futr_tmp <- data.frame(
      cov_mod   =  as.character(rownames(full.mat)[i]),
      sp        =  s,
      recTypenm =  factor(rectypes[tt],levels=rectypes),
      recType   =  tt,
      hypoth    =  h,
      hypothnm  = factor(hypotheses[tt],levels=hypotheses),
      model     =  "futR",
      i         =  i,
      log_a     =  tmp[tmp$def=="log_a",]$pred,
      log_b     =  tmp[tmp$def=="log_b",]$pred,
      logsigma  =  tmp[tmp$def=="logsigma",]$pred)
    
    futr_tmp                    <-  cbind(futr_tmp,matrix(0,1,ncov),matrix(0,1,ncov))
    colnames(futr_tmp)[-(1:11)] <-  c(paste0("beta[",1:ncov,"]"),paste0("lambda[",1:ncov,"]"))
    
    
    futr_tmpSE <- data.frame(
      cov_mod   =  as.character(rownames(full.mat)[i]),
      sp        =  s,
      recTypenm =  factor(rectypes[tt],levels=rectypes),
      recType   =  tt,
      hypoth    =  h,
      hypothnm  = factor(hypotheses[tt],levels=hypotheses),
      model     =  "futR",
      i         =  i,
      log_a     =  tmp[tmp$def=="log_a",]$pred.se,
      log_b     =  tmp[tmp$def=="log_b",]$pred.se,
      logsigma  =  tmp[tmp$def=="logsigma",]$pred.se)
    
    futr_tmpSE                    <-  cbind(futr_tmpSE,matrix(0,1,ncov),matrix(0,1,ncov))
    colnames(futr_tmpSE)[-(1:11)] <-  c(paste0("beta[",1:ncov,"]"),paste0("lambda[",1:ncov,"]"))
    
    #OKO FIX HERE - needs to be beta and Lambda!
    
    if(length(grep("beta",tmp$def))!=0){
      cc  <- grep("beta",tmp$def)
      futr_tmp[paste0("beta[",pnum,"]")]    <-  tmp[ cc,]$pred
      futr_tmpSE[paste0("beta[",pnum,"]")]  <-  tmp[ cc,]$pred.se
    }
    if(length(grep("lambda",tmp$def))!=0){
      cc  <- grep("lambda",tmp$def)
      futr_tmp[paste0("lambda[",pnum,"]")]    <-  tmp[ cc,]$pred
      futr_tmpSE[paste0("lambda[",pnum,"]")]  <-  tmp[ cc,]$pred.se
    }
    
    rm(tmp)
    
    pred           <-  futR_mods$pred
    pred$sp        <-  s
    pred$model     <-  nm
    pred$modelnum  <-  i
    pred$CE        <-  TRUE
    pred$rsType    <-  rectypes[tt]
    pred$hypoth    <-  hypotheses[h]
    pred$npar      <-  npar
    pred$n         <-  dim(rec)[1]
    pred$covset    <-  covset
    pred$covtxt    <-  covtxt
    Rhat           <-  pred[pred$def=="R_hat",]
    Shat           <-  pred[pred$def=="S_hat",]
    datIN          <-  futR_mods$model$env$data
    lmm            <-  lm(log(datIN$R_obs)~log(Rhat$pred))
    
    tmpaic<- 
      data.frame(
        LL    =   futR_mods$fit$objective,
        npar  =   npar,
        LnDet =   futR_mods$LnDet,
        n     =   dim(rec)[1],
        rsType =  rectypes[tt],
        hypoth =  hypotheses[h],
        name  =   nm,
        sp    =   s,
        CE    =   TRUE,
        conv  =   futR_mods$fit$convergence,
        R2    =   as.numeric(summary(lmm)["adj.r.squared"]),
        covset    =  covset,
        covtxt    =  covtxt)
     
    return(list(pred=pred, tmpaic=tmpaic, dfOut=df,futr_tmp=futr_tmp,futr_tmpSE=futr_tmpSE))
    rm(list=c("lmm","datIN","Shat","Rhat","pred","df","futr_tmp","futr_tmpSE"))
   
}