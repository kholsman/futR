#'profile_covars.R
#' 
#' profile the covariates to get partial effects
#' 
#' @param cov_nm covariate name
#' @param T2  T/F squared term for the parameter
#' @param SD_range    range of SD
#' @param sim_nitrIN  number of random simulations
#' @param simulateIN  simulate random iterations sim_nitrIN
#' @param modIN     rec model e.g., AIC_summry[[s]]$topRicker_R2
#' @param recIN     recruitment fit from the rec model, e.g., AIC_summry[[s]]$topRicker_R2_rec_fit
#' @param phasesIN  phases for each parm
#' @importFrom dplyr %>%
#' @importFrom ggplot2 aes
#' @importFrom stats nlminb
#' @importFrom stats cov
#' @importFrom stats quantile
#' 
#' @returns list(sim_df=sim_df,R_hat_sim=R_hat_sim, poly=aa,main = bb,p=p)
#' 
#' @export
#' 
#' 
#' 

profile_covars <- function(
    cov_nm     = "Spring_temp_surface5m",
    T2         = T,
    SD_range   =  3,
    sim_nitrIN = 100,
    simulateIN = TRUE,
    modIN,
    recIN,
    phasesIN = phases){
    

  tmpmod <- get_model(
    mod        = modIN,
    recIN      = recIN,
    simulateIN = FALSE,
    futR_fldr =  NULL,
    sim_nitrIN = sim_nitrIN,
    phasesIN   = phasesIN)

  # get objects from the fitted model
  obj <- tmpmod$object$model
  ee <- environment(obj$fn)
  lp <- ee$last.par.best  ## used in $report() call below
  dd <- ee$data           ## data object
  obj$report(lp)
  
  mu <- as.numeric(apply(dd$rs_cov,1,mean))
  sd <- as.numeric(apply(dd$rs_cov,1,sd))
  if(any(rownames(dd$rs_cov)==cov_nm)){
    covar_n  <- which(rownames(dd$rs_cov)==cov_nm)
  }else{
    stop(paste("covar names do not match",cov_nm))
  }
  
  if(T2){
    if(any(rownames(dd$rs_cov)==paste0(cov_nm,"x2"))){
      covar_n2 <-  which(rownames(dd$rs_cov)==paste0(cov_nm,"x2")) 
    }else{
      stop(paste0("covar names do not match ",cov_nm,"x2"))
    }
  }
  tt <- seq(-SD_range,SD_range,.1)*sd[covar_n]+mu[covar_n]
  nc <- dim(dd$rs_cov)
  
  ny <- length(tt)-1
  tt <- tt[1:ny]
  
  newdata           <- dd
  newdata$nyrs      <- ny
  newdata$rs_cov    <- matrix(mu,nc[1], ny,byrow=F)
  rownames(newdata$rs_cov) <- rownames(dd$rs_cov)
  newdata$beta_0    <- matrix( dd$beta_0[,1],nc[1], ny,byrow=F)
  newdata$lambda_0  <- matrix( dd$lambda_0[,1],nc[1], ny,byrow=F)
  
  newdata$rs_cov[covar_n,] <- tt
  if(T2) newdata$rs_cov[covar_n2,] <- tt*tt
  
  newdata$sdrs_cov <-newdata$rs_cov*0
  newdata$cov_type <- 
  newdata$sdR      <- 
  newdata$sdS      <- rep(0,ny)
  
  newdata$years    <- (1:ny)+1980
  newdata$S_obs    <- rep(mean(obj$report()$S_obs),ny)
  newdata$R_obs    <- rep(mean(obj$report()$R_obs),ny)
  
  parms <- ee$last.par.best
  newpar <- tmpmod$dlist$parameters
  for(i in 1:length(newpar)){
    if(any(names(parms)==names(newpar)[i]))
      newpar[[i]] <- as.numeric(parms[which(names(parms)==names(newpar)[i])])
  }
  report_out <- predict_futR(modelIN=tmpmod$object, newdata=newdata, newpar=newpar,src_fldr  = "src/TMB")
  report_out$cov <- tt
  bb <- data.frame(cov = report_out$cov,R_hat = report_out$R_hat,name=cov_nm)
  
  if(simulateIN){
    #tmpmod$simdat <- array(sim_nitr)
    sim <- replicate(sim_nitrIN, {
      simdata <- tmpmod$object$model$env$simulate(par = tmpmod$object$model$par, complete=TRUE)
      obj2    <- TMB::MakeADFun(data       = simdata,
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
    aa <- R_hat_sim%>%
      dplyr::group_by(cov)%>%
      dplyr::summarize(lower=quantile(R_hat,probs = .05),
                upper=quantile(R_hat,probs = .95),
                med=quantile(R_hat,probs = .5))%>%data.frame()
   
    p <- ggplot2::ggplot(aa)+
      ggplot2::geom_ribbon(aes(x=cov,ymin=lower,ymax=upper,fill="a_sim"),alpha=.4)+
      ggplot2::geom_line(aes(x=cov,y=med,color="a_sim"),size=1)+
      ggplot2:: ylab("recruitment")+
      ggplot2::xlab(cov_nm)+
      ggplot2::geom_line(data=bb,aes(x=cov,y=R_hat,color="b_main"),size=1)+
      ggplot2::theme_minimal()+scale_color_viridis_d(end=.7)+
      ggplot2::scale_fill_viridis_d(begin=.3,end=.8)
    
    return(list(sim_df=sim_df,R_hat_sim=R_hat_sim, poly=aa,main = bb,p=p))
    
  }else{
    
    return(main=bb)
  }
  
  
}
  
  
  
  
