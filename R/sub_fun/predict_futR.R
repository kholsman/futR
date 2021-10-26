#   newdata<-data.frame(SSB=seq(1,7e6,1e6),sdR=0)
# object<-  mod$final
# dlistIN<-mod$dlistfinal
predict_futR<-function(object,se.fit=FALSE,newdata=NULL,simulate=TRUE,nitr=100){
  if(!is.null(newdata)){
    fl       <-  c("SSB","Robs","Ration_scaled","years","rs_cov")
    fl       <-  names(object$datlist$rs_dat)
    fl2      <-  c("epsi_s")
    dlnew    <-  final$datlist
    nd       <-  newdata
    datlist  <-  object$datlist
    fl_n     <-  data.frame(r=rep(0,length(fl)),c=0)
    
    # get dims for each fl object:
    for(i in 1:length(fl))
      fl_n[i,]<-dim(as.matrix(dlnew$rs_dat[[i]]))
    
    # if there is new data:
    if(!is.null(nd)){
      if(!any(names(nd)%in%c(fl,fl2))) stop("newdata names do not match obj data")
      ml   <-  which(fl%in%names(nd))
      ny   <-  length(nd$SSB)
      nc   <-  dlnew$rs_dat$ncov
      
      # set new nyrs(nobs):
      dlnew$rs_dat$nyrs  <-  ny
      
      # set other data inputs to mean vals:
      for(i in (1:(length(fl)))[-ml]){
        nm<-fl[i]
        if( sum(fl_n[i,]) == 2){
          # if scalar
          if(!nm%in%names(nd))
            dlnew$rs_dat[[nm]]  <-  dlnew$rs_dat[[nm]]
        }else{
          
          if(fl_n[i,1] > 1 & fl_n[i,2] > 1){
            # if matrix
            if(!nm%in%names(nd))
              dlnew$rs_dat[[nm]]<-matrix(apply(dlnew$rs_dat[[nm]],1,mean),nc,ny)
          }else{
            # if vector
            if(!nm%in%names(nd))
              dlnew$rs_dat[[nm]]<-rep(mean(dlnew$rs_dat[[nm]]),ny)
          }   
        }
        
      }
      
      # now fill in the rest with the new data:
      for(nm in fl[ml])
        dlnew$rs_dat[[nm]]<-newdata[[nm]]
    }
    
    # # now for parameters: DEFUNCT 
    
    # fl<-fl2
    # # if there is new data:
    # if(!is.null(nd)){
    #   ml    <-  which(fl%in%names(nd))
    #   for(nm in fl)
    #     if(!nm%in%names(nd))
    #       dlnew$parameters[[nm]]<-dlnew$parameters[[nm]]
    #      # dlnew$parameters[[nm]]<-rep(mean(dlnew$parameters[[nm]]),ny)
    
    #    for(nm in fl[ml])
    #     dlnew$parameters[[nm]]<-newdata[[nm]]
    # }
    
    # re make the datlist with the new data:
    tt<- makeDat(
      dataIN    = dlnew$rs_dat,
      typeIN    = datlist$inputs$typeIN ,
      tauIN     = datlist$inputs$tauIN,
      estparams = dlnew$estparams,
      covarsIN  = dlnew$rs_dat$rs_cov,
      startVal  = dlnew$parameters,
      phases    = NULL,
      fityrs    = dlnew$rs_dat$years,
      rationIN  = dlnew$rs_dat$Ration_scaled,
      REcompile = FALSE,
      Eat_cov   = FALSE,
      sigMethod = datlist$inputs$sigMethod)
    
    dlnew$maplist<-tt$maplist
    
  }else{
    dlnew   <-  final$datlist
    
  }
  newObj  <-  MakeADFun(
    data=dlnew$rs_dat,
    parameters=dlnew$parameters,
    DLL='futR',
    checkParameterOrder=TRUE,
    hessian = TRUE,
    map=dlnew$maplist,
    profile = NULL, # TODO: Optionally "beta"
    silent = TRUE)
  
  oldPar  <-  object$mle
  newObj$fn(oldPar)    ## call once to update internal structures
  lp      <-  newObj$env$last.par
  obj     <-  object$model
  #obj$fn(oldPar) 
  
  stopifnot(identical(oldPar,lp))   
  newObj$simulate(par=object$mle,complete=TRUE)$R_hat
  
  if (!se.fit) {
    pred   <- unlist(newObj$report(lp))
  } else {
    H      <- with(newObj,optimHess(oldPar,newObj$fn,newObj$gr))
    sdr    <- sdreport(newObj,oldPar,hessian.fixed=H,getReportCovariance=F)
    sdrsum <- TMB:::summary.sdreport(sdr, "report")
    pred   <- sdrsum[,"Estimate"]
    se     <- sdrsum[,"Std. Error"]
  }
  if(simulate){
    sim <- replicate(nitr, {
      new<-newObj$simulate(par=oldPar, complete=TRUE)
      unlist(new)
    })
    if (!se.fit) return(list(pred=data.frame(def=names(pred),fit=pred,se.fit=NA),sim=sim)) else return(list(pred=data.frame(def=names(pred),fit=pred,se.fit=se),sim=sim))
    
  }else{
    if (!se.fit) return(data.frame(def=names(pred),fit=pred,se.fit=NA)) else return(data.frame(def=names(pred),fit=pred,se.fit=se))
  }
  
  
}
