#'
#'
#' profile_SSB.R
#' 
#' profile the SSB to get partial effects
#' 
#' 

profile_ssb <- function(
    cov_nm_list = c("Spring_temp_surface5m","Summer_temp_surface5m","Winter_temp_bottom5m"),
    nameIN      = "hindcast 2020",
    newdataIN   = hind%>%filter(year==2020)%>%data.frame(),
    T2          = T,
    SD_range   = 7,
    steps      = 100,
    sim_nitrIN = 100,
    simulateIN = TRUE,
    hind_yearsIN = 1980:2023,
    modIN    = AIC_summry[[s]]$topRicker_R2,
    recIN    = AIC_summry[[s]]$topRicker_R2_rec_fit,
    phasesIN = phases){
  
  tmpmod <- get_model(
    mod        = modIN,
    recIN      = recIN,
    simulateIN = FALSE,
    futR_fldr  = NULL,
    sim_nitrIN = sim_nitrIN,
    phasesIN   = phasesIN)
  
  # get objects from the fitted model
  obj <- tmpmod$object$model
  ee <- environment(obj$fn)
  lp <- ee$last.par.best  ## used in $report() call below
  dd <- ee$data           ## data object
  obj$report(lp)
  nc <- dim(dd$rs_cov)  
  
  mu <- as.numeric(apply(dd$rs_cov,1,mean))
  sd <- as.numeric(apply(dd$rs_cov,1,sd))
  muSSB <-  mean(obj$report()$S_obs)
  sdSSB <-  sd(obj$report()$S_obs)
  j <- jj <-jjj <- 0
  covar_n <- covar_n_new <- covar_n2 <- newD <- rep(0,nc[1])
  
  for(cov_nm in cov_nm_list){
    if(any(rownames(dd$rs_cov)==cov_nm)){
      j <- j+1
      covar_n[j]  <- which(rownames(dd$rs_cov)==cov_nm)
    }else{
      stop(paste("covar names do not match",cov_nm))
    }
    if(any(newdataIN$var==cov_nm)){
      jjj <- jjj+1
      covar_n_new[jjj]  <- which( (newdataIN$var==cov_nm) )
      newD[covar_n[j]]  <- newdataIN[covar_n_new[jjj],]$val_use_scaled
    }else{
      stop(paste("covar_n_new names do not match",cov_nm))
    }
   
    if(T2){
      if(any(rownames(dd$rs_cov)==paste0(cov_nm,"x2"))){
        jj <- jj+1
        covar_n2[jj]       <-  which(rownames(dd$rs_cov)==paste0(cov_nm,"x2")) 
        newD[covar_n2[jj]] <-  newdataIN[covar_n_new[jjj],]$val_use_scaled^2
      }
      
    }
    
  }
  if(T2)
    if(jj==0)
      stop(("no covar names match x2"))
  
  tt      <- seq(0,SD_range*sdSSB+muSSB,round(muSSB/steps))
  ny      <- length(tt)-1
  tt      <- tt[1:ny]
  newdata <- dd

  # profile across SSB:
  newdata$S_obs     <- tt 
  newdata$R_obs    <- rep(mean(obj$report()$R_obs),ny)
  newdata$nyrs      <- ny
  newdata$rs_cov    <- matrix(newD,nc[1], ny,byrow=F)
  rownames(newdata$rs_cov) <- rownames(dd$rs_cov)
  newdata$beta_0    <- matrix( dd$beta_0[,1],nc[1], ny,byrow=F)
  newdata$lambda_0  <- matrix( dd$lambda_0[,1],nc[1], ny,byrow=F)
  
  newdata$sdrs_cov <- newdata$rs_cov*0
  newdata$cov_type <- 
    newdata$sdR      <- 
    newdata$sdS      <- rep(0,ny)
  newdata$years    <- (1:ny)+1980
  
  parms  <- ee$last.par.best
  newpar <- tmpmod$dlist$parameters
  for(i in 1:length(newpar)){
    if(any(names(parms)==names(newpar)[i]))
      newpar[[i]] <- as.numeric(parms[which(names(parms)==names(newpar)[i])])
  }
  report_out <- predict_futR(modelIN=tmpmod$object, newdata=newdata, newpar=newpar,src_fldr  = "src/TMB", futR_fldr = NULL)
  report_out$SSB  <- tt
  report_out$year <- newdata$years
  bb <- data.frame(SSB = report_out$S_obs,R_hat = report_out$R_hat,muSSB=muSSB,sdSSB=sdSSB,name=nameIN,type = "profile", year = report_out$year)

  if(simulateIN){
    #tmpmod$simdat <- array(sim_nitr)
    sim <- replicate(sim_nitrIN, {
      simdata <- tmpmod$object$model$env$simulate(par = tmpmod$object$model$par, complete=TRUE)
      obj2    <- MakeADFun(data       = simdata,
                           parameters = tmpmod$dlist$parameters,
                           DLL        = tmpmod$object$model$env$DLL,
                           checkParameterOrder=TRUE,
                           hessian    = TRUE,
                           map        = tmpmod$dlist$maplist,
                           profile    = NULL, # TODO: Optionally "beta"
                           silent     = TRUE)
      list(par =nlminb(obj2$par, obj2$fn, obj2$gr)$par,model=obj2)
    })
    for(jj in 1:sim_nitrIN){
      print(jj)
      tmpmod_sim <- list(model=sim[,jj]$model)
      parms      <- sim[,jj]$par
      newpar     <- tmpmod$dlist$parameters
      for(i in 1:length(newpar))
        if(any(names(parms)==names(newpar)[i]))
          newpar[[i]] <- as.numeric(parms[which(names(parms)==names(newpar)[i])])
      
      report_out <- suppressWarnings(predict_futR(modelIN=list(model=sim[,jj]$model), 
                                                  newdata=newdata, 
                                                  newpar=newpar,src_fldr  = "src/TMB"))
      
      if(jj==1){
        
        sim_df    <- data.frame(estimate=as.vector(sim[,jj]$par), parameter=names(sim[,jj]$par),itr=jj)
        R_hat_sim <- data.frame(cov = newdata$rs_cov[covar_n,],R_hat = report_out$R_hat,itr=jj)
        
      }else{
        sim_df <- rbind(sim_df,data.frame(estimate=as.vector(sim[,jj]$par),
                                          parameter=names(sim[,jj]$par),
                                          itr=jj))
        R_hat_sim <- rbind(R_hat_sim,data.frame(cov = newdata$rs_cov[covar_n,],R_hat = report_out$R_hat,itr=jj))
      }
      
    }
    aa<-R_hat_sim%>%group_by(cov)%>%
      summarize(lower=quantile(R_hat,probs = .05),
                upper=quantile(R_hat,probs = .95),
                med=quantile(R_hat,probs = .5))%>%data.frame()
    
    p <- ggplot(aa)+
      geom_ribbon(aes(x=cov,ymin=lower,ymax=upper,fill="a_sim"),alpha=.4)+
      geom_line(aes(x=cov,y=med,color="a_sim"),size=1)+
      ylab("recruitment")+
      xlab(cov_nm)+
      geom_line(data=bb,aes(x=cov,y=R_hat,color="b_main"),size=1)+
      theme_minimal()+scale_color_viridis_d(end=.7)+scale_fill_viridis_d(begin=.3,end=.8)
    
    return(list(sim_df=sim_df,R_hat_sim=R_hat_sim, poly=aa,main = bb,p=p))
    
  }else{
    
    return(main=bb)
  }
  
  
}




