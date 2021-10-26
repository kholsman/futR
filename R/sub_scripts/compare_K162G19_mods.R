# Kirstin Holsman & Grant Adams
# kirstin.holsman@noaa.gov
# CEfutR

# 1. Load data
# 2. fit futR to ceattle outputs
# 3. compare futR to ceattle
# 4. create projection files (.dat for ceattle, or futdat for Rceattle)
# 5. Project under HCRS
  
  #___________________________________________
  # 1. Set things up
  #___________________________________________
  
    rm(list=ls())
    
    # 1.0. define directories
    main      <-  path.expand("~/GitHub_new/futR")
    mod_fldr  <-  path.expand("~/Documents/D_AFSC_Files/Manuscripts/00_ROMS_CEATTLE_RS/data/runs/aclim_00_JunV2_2019_2")
    rec_fldr  <-  file.path(mod_fldr,"Recruitment_files")
    
    setwd(main)
    file.copy(from= file.path(rec_fldr,"Rec_figs/recruitment.Rdata"),to="data/recruitment.Rdata",overwrite = TRUE)

    # 1.1. load libraries
    library( Rceattle )
    library( TMB )
    library(reshape)
    library(dplyr)
    library(ggplot2)
    library(lattice)
    
    # 1.2 input switches
    run_est    <- FALSE
    REcompile  <- TRUE
    bias_corr  <- TRUE # bias correct sigma sensu Ludwig and Walters 1981 and Porch and Lauretta 2016
    recruitAge <- 1    # what age is recruitment estimated in the model?
    fityrsIN   <- 1980:2017
    Eat_covIN  <- c(18, 20, 26, 28)
    Eat_covIN  <- -99

    # 1.3 load data
    
   # load("data/aclim_cmip5_BT.Rdata")  # Bottom Temp object called allDat
    load("data/env_covars.Rdata")  # Bottom Temp object
    load("data/NPZcov.Rdata")      # key scaled indices from ROMSNPZ (by year, not yet lagged)
    load("data/recruitment.Rdata")
    load( "data/ms_run.Rdata" )
    rec_yrs  <- yrs
    
    # 1.4 load functions
    source("R/futRfun.R")
    
    GET_ESTM2 <- function(mode, s = 1){
      plotdat     <-  smry_RS[[s]]
      Rectxt      <-  Rec_2$TopR2.txt[[s]]
      Rectxt_aic  <-  Rec_2$TopAIC.txt[[s]]
      rec_path    <-  rec_fldr
      estnm       <-  file.path(mod_fldr,"results/ceattle_est.std")
      rsdat       <-  rs_data
      covv        <-  covuse.all
      return(list(covv = covv,  Rectxt = Rectxt, Rectxt_aic = Rectxt_aic,rsdat=rsdat,estnm=estnm,rec_path=rec_path,plotdat=plotdat))
    }
    
    # 1.5 set up mode for future R
    mode = 1  # 1 = estimate simga, similar to ceattle approach
    
    # this will generate warnings - they can be ignored if "0" is returned
    compile('futR.cpp') 
    
    # 4.1 set up covars for epsi_s # epsi_svironmental variabile from the age class year (assuming Jan rec, and summer = age0)
    
    #years        <-  ms_run$data_list$styr-1+(1:(ms_run$data_list$nyrs))[-1]
    ALLcovars    <-  ROMS_NPZ_covars[["aclim_hindcast.use"]]
    rownames(ALLcovars)  <- ALLcovars$year
    covars       <-  ALLcovars[which(ALLcovars$year%in%fityrsIN)-recruitAge,]
    topmod       <-  which(names(covars)%in%c("BottomTemp","ColdPool","fallZavg","springZavg"))
    topmod       <-  which(names(covars)%in%c("Btm_sal_tmp", "ColdPool" , "fallZavg" , "springZavg" ))
    modc         <-  which((AICtable[[1]]$model)=="Btm_sal_tmp_ColdPool_fallZavg_springZavg_") 

    PAR<-data.frame(
      phases = c(
        log_a        = 1, 
        log_b        = 1 ,
        #logit_tau    = 1,
        gamma        = 1,
        beta         = 1,
        lambda       = 1,
        epsi_s       = 1,
        logsigma     = 1),
      estparams  = c(
        log_a        = TRUE, 
        log_b        = TRUE, 
        #logit_tau    = FALSE,
        gamma        = FALSE,
        beta         = FALSE,
        lambda       = TRUE,
        epsi_s       = FALSE,
        logsigma     = TRUE))
    
    cc<-c(6,9,20,28)
    TopR2_Ricker.txt
    covars[,topmod]
    cc<-topmod-1
    
    tmp                   <-   as.matrix((t(covars[,-1][,cc])))
    sub_cov               <-   matrix(0,dim(tmp)[1],dim(tmp)[2])
    
    for(rr in 1:dim(tmp)[1])
      sub_cov[rr,]        <-   (tmp[rr,])
    colnames(sub_cov)     <-   covars[,1]
    nm                    <-   rownames(tmp)
    
    env_covars<-sub_cov[c(1,2),]
    
    # set up recruitment input data file Years, SSB, Robs, sdSSB, sdRobs
    rec_dat  <-  rec_dat_tmb  <-  list()
    for(sp in 1:3){
      est        <-  ms_run
      smry_RSIN  <-  smry_RS[[sp]]
      sdlist     <-  as.list(est$sdrep,"Std. Error")
      estlist    <-  as.list(est$sdreNAp,"Estimate")
      mpRsd      <-  sdlist$rec_dev[sp,-1]
      #mpRsd    <-  mpRsd*0
      
      dataIN_tmb <-data.frame(
        years  = est$data_list$styr-1+(1:(est$data_list$nyrs))[-1],
        S_obs    = est$quantities$biomassSSB[sp,-length(est$quantities$biomassSSB[sp,])],
        R_obs   = est$quantities$R[sp,-1],
        sdS     = 0,
        sdR     = mpRsd)
      
      
      dataIN_admb<-data.frame(
        years   = smry_RSIN$Rec_years,
        S_obs   = smry_RSIN$SSB.y.1.,
        R_obs   = smry_RSIN$R_obs,
        sdS     = 0,
        sdR     = mpRsd) 
      
      rec_dat[[sp]]<- dataIN_admb
      rec_dat_tmb[[sp]]<- dataIN_tmb
      if(sp==1)
        ration_tmb<- data.frame(sp1=as.numeric(scale(apply(est$quantities$ration2Age[-sp,-1,],3,sum)))[-nyrs],sp2=NA,sp3=NA)
      if(sp>1)
        ration_tmb[,sp]    <-  (as.numeric(scale(apply(est$quantities$ration2Age[-sp,-1,],3,sum)))[-nyrs])
    }
    
    # save(env_covars,file="data/env_covars.Rdata")
    # save(rec_dat,file="data/rec_dat.Rdata")
    # save(ration_tmb,file="data/ration_tmb.Rdata")
    
    #___________________________________________
    # 4. COMPILE REC MODEL
    #___________________________________________
    # General example
    
    rec     <-  rec_dat[[1]]
    env     <-  env_covars
    env[1,]    <-  as.numeric(scale(env_covars[1,]))
    env[2,]    <-  as.numeric(scale(env_covars[2,]))
    ration  <-  ration_tmb[,1]
    version <-  "futR"
    
    tau_set<-c(0,.01,.1,0.5,.8,1,1.2,1.5,2,5,10)
    
    
    # 4.2 Set up data
    PAR$phases
    PAR$estparams
    estparams  = c(
      log_a        = TRUE,
      log_b        = TRUE,
      #logit_tau     = TRUE,
      gamma        = FALSE,
      beta         = FALSE,
      lambda       = TRUE,
      epsi_s       = FALSE,
      logsigma     = TRUE)
    #estparams<- PAR$estparams
    rec_noerr<-rec
    rec_noerr$sdRobs<-0
    
    #__________________________________________________________________
    # rectype:
    # 1 = Linear with biomass ( y-1 )  BLM
    # 2 = Beverton Holt
    # 3 = Ricker
    # 4 = exponential  
    
    # sigMethod:
    # 1 = estimate simga
    # 2 = unbiased sigma empirically calculated (rather than biased estimate from MLE; sensu Ludwig and Walters):
    # 3 = as in 1 but with specified measurement error for rec (indep of random effects on S)
    # 4 = as in 1 but with specified measurement error for rec and SSB   
    
    # tMethod:
    # 1 = cloglog link (g = 1-exp(-exp(gamma)))
    # 2 = logit link (g = exp(gamma)/(1+exp(gamma)))
    
    # environmental effects:
    # beta = effects on effective number of spawners
    # lamba = effects on post-spawning success (e.g., age 0 survival)
    
    
    #__________________________________________________________________
    Rec_2  <- list(TopR2.txt=TopR2.txt,TopAIC.txt=TopAIC.txt)
   
   
    
    
    #__________________________________________________________________
    # get parm vals from ceattle model:
    rsdat     <-  
      Rectxt    <-
      Rectxt_aic <- 
      estnm     <-  
      rec_path  <-
      splist    <-
      rs_aic    <-
      rs        <-  list()
    
    
    for(ss in 1:nspp){
      splist[[ss]]    <-  c("pollock","Pacific cod","arrowtooth")[ss]
      ptt             <-  GET_ESTM2(mode=mode,s=ss)
      covv            <-  ptt$covv
      rsdat[[ss]]     <-  ptt$rsdat
      Rectxt[[ss]]    <-  ptt$Rectxt
      Rectxt_aic[[ss]]    <-  ptt$Rectxt_aic
      estnm[[ss]]     <-  ptt$estnm
      rec_path[[ss]]  <-  ptt$rec_path
      if(ss==1) plotdat  <-  data.frame(Species=ss,ptt$plotdat)
      if(ss!=1) plotdat  <-  merge(plotdat,data.frame(Species=ss,ptt$plotdat),
                                   by=c("Species","spp","Rec_years","SSB.y.1.","R_obs","mnRec","mnRS_Rec" ,"TopR2_Rec","TopAIC_Rec","mnRecwB"),all=TRUE)
      
      dirlist        <-  dir(file.path(rec_path[[ss]],"RS_fits"))
      #________________________________
      # for top aic model
      rsm       <-  NA
      if(length(grep("LM",Rectxt_aic[[ss]]))!=0){
        rsm  <- "LM"
        if(length(grep("BLM",Rectxt_aic[[ss]]))!=0)
          rsm  <- "BLM"
      }
       
      cc             <-  grep(Rectxt_aic[[ss]],dirlist)
      flnm           <-  dirlist[cc[grep(".std",dirlist[cc])]]
      
      if(is.na(rsm)){
        flnm           <-  flnm[-grep("LM",flnm)]
        flnm           <-  flnm[grep(paste0(Rectxt_aic[[ss]],".std"),flnm)]
      }
     
      rs_aic[[ss]]              <-  read.csv(file.path(rec_path[[ss]],"RS_fits",flnm),header=F, skip=1,sep="")
      colnames( rs_aic[[ss]])   <-  c("index","name","value","std")
      rs_aic[[ss]]$var_nm<-""
      
      #________________________________
      # now for RS function (top R2)
      rsm       <-  NA
      if(length(grep("LM",Rectxt[[ss]]))!=0){
        rsm  <- "LM"
        if(length(grep("BLM",Rectxt[[ss]]))!=0)
          rsm  <- "BLM"
      }
      
      cc             <-  grep(Rectxt[[ss]],dirlist)
      flnm           <-  dirlist[cc[grep(".std",dirlist[cc])]]
      
      if(is.na(rsm)){
        flnm           <-  flnm[-grep("LM",flnm)]
        flnm           <-  flnm[grep(paste0(Rectxt[[ss]],".std"),flnm)]
      }
      
      rs[[ss]]              <-  read.csv(file.path(rec_path[[ss]],"RS_fits",flnm),header=F, skip=1,sep="")
      colnames( rs[[ss]])   <-  c("index","name","value","std")
      
    }
    
    # First run with no observation error on SSB:
    # tauIN = 1; pre-specified tau if =0 then no ranom effects are included
    estparams$gamma <- TRUE
    datlist<-makeDat(tauIN=1,sigMethod=mode, tMethod=1,estparams=estparams,typeIN=4,dataIN=rec,
                     covarsIN  = env)
    mm1<-fitRec(dlistIN=datlist,version="futR",recompile=T,simulate=TRUE)
    
    df <- data.frame(estimate=as.vector(mm1$sim), parameter=names( mm1$mle)[row(mm1$sim)])
    densityplot( ~ estimate | parameter, data=df, layout=c(5,1),ylim=c(0,10))
    
    datlist$rs_dat$tau<-0.000001
    mm1_t0<-fitRec(dlistIN=datlist,version="futR",recompile=F)
    
    get_par_vals<-function(fl,ss){
      dirlist        <-  dir(file.path(rec_path[[ss]],"RS_fits"))
      
      
      rsm       <-  NA
      if(length(grep("LM",fl))!=0){
        rsm  <- "LM"
        if(length(grep("BLM",fl))!=0)
          rsm  <- "BLM"
      }
      
      cc             <-  grep(fl,dirlist)
      flnm           <-  dirlist[cc[grep(".std",dirlist[cc])]]
      
      if(is.na(rsm)){
        flnm           <-  flnm[-grep("LM",flnm)]
        flnm           <-  flnm[grep(paste0(fl,".std"),flnm)]
      }
      
      rs_out              <-  read.csv(file.path(rec_path[[ss]],"RS_fits",flnm),header=F, skip=1,sep="")
      colnames( rs_out)   <-  c("index","name","value","std")
      return(rs_out)
      
    }    
    
  #___________________________________________________________
  # RUN futR models
  #___________________________________________________________ 
      
    # for each rec model create a matrix of ncov X nyears:
    # rectype:
        # 1 = Linear with biomass ( y-1 )  BLM
        # 2 = Linear LM
        # 3 = Beverton Holt 
        # 4 = ricker 
    rectypes  <- c("BLM","LM","BH","Ricker")
    futR_mods <- list()
    nmod      <- dim(full.mat)[1]  # number of covariate models
    ncov      <- dim(full.mat)[2]
    nspp      <- 3  # number of species
    simIT     <- FALSE
    silentINN <- TRUE
    
    for(s in 1:3){
      rec     <-  rec_dat[[s]]
      print(paste("--------------- species ", s,"--------------------"))
      futR_mods[[s]]<-list()
      for(tt in 1:length(rectypes)){
        print(paste("--------------- rectype ", rectypes[tt],"--------------------"))
        futR_mods[[s]][[tt]]<-list()
        for(i in 1:nmod){
         
          
          print(paste("model ", i, "of", nmod,"rectype", rectypes[tt],"species ", s))
          pnum     <-  which(full.mat[i,]==1)
          tmpenv   <-  t(covs[,pnum])
          tmpenv2  <-  t(covars[,-1][,pnum])
          if(length(pnum)>1){
            covtxt   <-  paste(names(covs[,pnum]),collapse=", ")
          }else{
            covtxt   <-  names(covs)[pnum]
          }
          covset   <-  paste0("c(",paste(pnum,collapse=", "),")")
          
          if(any((tmpenv2-tmpenv)>0.0001))
            message("covs and covars don't match (covars different than those used to fit ceattle rec)")
         
          txt         <- paste0("RS",s,"_",rownames(full.mat)[i])
          tmppar      <- get_par_vals(fl=txt,ss=s)
          ceattle_tmp <- data.frame(
            cov_mod   =  as.character(rownames(full.mat)[i]),
            sp        =  sp,
            recTypenm =  factor(rectypes[tt],levels=rectypes),
            recType   =  tt,
            model     =  "futR",
            i         =  i,
            log_a     =  tmppar[tmppar$name=='log_aa_c',]$value,
            log_b     =  tmppar[tmppar$name=='log_bb_c',]$value,
            logsigma  =  tmppar[tmppar$name=='logsigma',]$value)
          ceattle_tmp <-cbind(ceattle_tmp,matrix(0,1,ncov))
          colnames(ceattle_tmp)[-(1:9)]<-paste0("rs_parm[",1:ncov,"]")
          
          ceattle_tmpSE <- data.frame(
            cov_mod   =  as.character(rownames(full.mat)[i]),
            sp        =  sp,
            recTypenm =  factor(rectypes[tt],levels=rectypes),
            recType   =  tt,
            model     =  "futR",
            i         =  i,
            log_a     =  tmppar[tmppar$name=='log_aa_c',]$std,
            log_b     =  tmppar[tmppar$name=='log_bb_c',]$std,
            logsigma  =  tmppar[tmppar$name=='logsigma',]$std)
          ceattle_tmpSE <-cbind(ceattle_tmpSE,matrix(0,1,ncov))
          colnames(ceattle_tmpSE)[-(1:9)]<-paste0("rs_parm[",1:ncov,"]")
          
           
          if(length(grep("rs_parm",tmppar$name))!=0){
            cc <- grep("rs_parm",tmppar$name)
            ceattle_tmp[as.character(tmppar$name[cc])]   <-  tmppar$value[cc]
            ceattle_tmpSE[as.character(tmppar$name[cc])] <-  tmppar$std[cc]
            txt1  <- as.character(tmppar$name[cc])
          }
          
          
          datlist<-makeDat(tauIN=1,sigMethod=mode, 
                           estparams=estparams,typeIN=tt,
                           dataIN=rec,
                           covarsIN  = tmpenv2)
          
          futR_mods[[s]][[tt]][[i]]<-fitRec(dlistIN=datlist,
                                            silentIN = silentINN,
                                            version="futR",
                                            recompile=F,
                                            simulate=simIT)
          
          npar           <-  length(which(!is.na(as.numeric(melt(futR_mods[[s]][[tt]][[i]]$model$env$map)$value))))
          if(rectypes[tt]=="LM")
            npar<-npar-1
          if(simIT){
            require(reshape)
            df <- melt(futR_mods[[s]][[tt]][[i]]$sim); colnames(df)<-c("parameter","itr","estimate")
            df$sp        <-  s
            df$model     <-  nm
            df$modelnum  <-  i
            df$CE        <-  TRUE
            df$npar      <-  npar
            df$n         <-  dim(rec)[1]
            df$rsType    <-  rectypes[tt]
            df$covset    <-  covset
            df$covtxt    <-  covtxt
          }
          
          
          tmppar2<- futR_mods[[s]][[tt]][[i]]$fit$par
          tmppar2<- futR_mods[[s]][[tt]][[i]]$pred
          futr_tmp <- data.frame(
            cov_mod   =  as.character(rownames(full.mat)[i]),
            sp        =  sp,
            recTypenm =  factor(rectypes[tt],levels=rectypes),
            recType   =  tt,
            model     =  "futR",
            i         =  i,
            log_a     =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="log_a",]$pred,
            log_b     =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="log_b",]$pred,
            logsigma  =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="logsigma",]$pred)
          futr_tmp <-cbind(futr_tmp,matrix(0,1,ncov))
          colnames(futr_tmp)[-(1:9)]<-paste0("rs_parm[",1:ncov,"]")
          
          futr_tmpSE <- data.frame(
            cov_mod   =  as.character(rownames(full.mat)[i]),
            sp        =  sp,
            recTypenm =  factor(rectypes[tt],levels=rectypes),
            recType   =  tt,
            model     =  "futR",
            i         =  i,
            log_a     =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="log_a",]$pred.se,
            log_b     =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="log_b",]$pred.se,
            logsigma  =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="logsigma",]$pred.se)
          futr_tmpSE <-cbind(futr_tmpSE,matrix(0,1,ncov))
          colnames(futr_tmpSE)[-(1:9)]<-paste0("rs_parm[",1:ncov,"]")
          
          
          if(any(!(paste0("rs_parm[",pnum,"]")%in%txt1)))
            message( "parm numbers don't match!")
          
          if(length(grep("rs_parm",futR_mods[[s]][[tt]][[i]]$pred$def))!=0){
            
            cc  <- grep("rs_parm",futR_mods[[s]][[tt]][[i]]$pred$def)
            futr_tmp[paste0("rs_parm[",pnum,"]")]    <-  futR_mods[[s]][[tt]][[i]]$pred[ cc,]$pred
            futr_tmpSE[paste0("rs_parm[",pnum,"]")]  <-  futR_mods[[s]][[tt]][[i]]$pred[ cc,]$pred.se
            
          }
          
          nm <- paste0("RS",s,"_",rownames(full.mat)[i],"_",rectypes[tt])
          if(rectypes[tt]=="Ricker")
            nm <- paste0("RS",s,"_",rownames(full.mat)[i])
          
          
          pred           <-  futR_mods[[s]][[tt]][[i]]$pred
          pred$sp        <-  s
          pred$model     <-  nm
          pred$modelnum  <-  i
          pred$CE        <-  TRUE
          pred$rsType    <-  rectypes[tt]
          pred$npar      <-  npar
          pred$n         <- dim(rec)[1]
          pred$covset    <-  covset
          pred$covtxt    <-  covtxt
          Rhat   <-  pred[pred$def=="R_hat",]
          Shat   <-  pred[pred$def=="S_hat",]
          datIN  <-  futR_mods[[s]][[tt]][[i]]$model$env$data
          lmm    <- lm(datIN$Robs~Rhat$pred)
          
          tmpaic<- 
              data.frame(
                LL    =   futR_mods[[s]][[tt]][[i]]$fit$objective,
                npar  =   npar,
                LnDet =   futR_mods[[s]][[tt]][[i]]$LnDet,
                n     =   dim(rec)[1],
                rsType =  rectypes[tt],
                name  =   nm,
                sp    =   s,
                CE    =   TRUE,
                conv  =   futR_mods[[s]][[tt]][[i]]$fit$convergence,
                R2    =   as.numeric(summary(lmm)["adj.r.squared"]),
                covset    =  covset,
                covtxt    =  covtxt)
          
          if(s==1&i==1&tt==1){
            predOut  <-  pred
            aicOut   <-  tmpaic
            if(simIT) dfOut    <-  df
            par_mat_ceattle    <- ceattle_tmp
            par_mat_futR       <- futr_tmp
            par_mat_ceattleSE  <- ceattle_tmpSE
            par_mat_futRSE     <- futr_tmpSE
          }else{
            predOut  <-  rbind(predOut,pred)
            aicOut   <-  rbind(aicOut,tmpaic)
            if(simIT) dfOut    <-  rbind(dfOut,df)
            par_mat_ceattle    <- rbind(par_mat_ceattle,ceattle_tmp)
            par_mat_futR       <- rbind(par_mat_futR,futr_tmp)
            par_mat_ceattleSE  <- rbind(par_mat_ceattleSE,ceattle_tmpSE)
            par_mat_futRSE     <- rbind(par_mat_futRSE,futr_tmpSE)
          }
          
        }
        
        names(futR_mods[[s]][[tt]])<-paste0("RS",s,"_",rownames(full.mat),"_",rectypes[tt])
      }
      names(futR_mods[[s]])<-rectypes
    }

  
  #____________________________________________
  #now run no covar versions and null
  #____________________________________________
      i  <- i + 1
      for(s in 1:3){
        print(paste("--------------- species ", s,"--------------------"))
        for(tt in 1:length(rectypes)){
          print(paste("--------------- rectype ", rectypes[tt],"--------------------"))
            
            covtxt   <-  ""
            covset   <-  "nocov"
            
            if(any((tmpenv2-tmpenv)>0.0001))
              message("covs and covars don't match (covars different than those used to fit ceattle rec)")
            txt<-paste0("RS",s,"_")
           
            ceattle_tmp <- data.frame(
              cov_mod   =  as.character(rownames(full.mat)[i]),
              sp        =  sp,
              recTypenm =  factor(rectypes[tt],levels=rectypes),
              recType   =  tt,
              model     =  "futR",
              i         =  i,
              log_a     =  tmppar[tmppar$name=='log_aa_c',]$value,
              log_b     =  tmppar[tmppar$name=='log_bb_c',]$value,
              logsigma  =  tmppar[tmppar$name=='logsigma',]$value)
            ceattle_tmp <-cbind(ceattle_tmp,matrix(0,1,ncov))
            colnames(ceattle_tmp)[-(1:9)]<-paste0("rs_parm[",1:ncov,"]")
            
            ceattle_tmpSE <- data.frame(
              cov_mod   =  as.character(rownames(full.mat)[i]),
              sp        =  sp,
              recTypenm =  factor(rectypes[tt],levels=rectypes),
              recType   =  tt,
              model     =  "futR",
              i         =  i,
              log_a     =  tmppar[tmppar$name=='log_aa_c',]$std,
              log_b     =  tmppar[tmppar$name=='log_bb_c',]$std,
              logsigma  =  tmppar[tmppar$name=='logsigma',]$std)
            ceattle_tmpSE <-cbind(ceattle_tmpSE,matrix(0,1,ncov))
            colnames(ceattle_tmpSE)[-(1:9)]<-paste0("rs_parm[",1:ncov,"]")
            
            
            if(length(grep("rs_parm",tmppar$name))!=0){
              cc <- grep("rs_parm",tmppar$name)
              ceattle_tmp[as.character(tmppar$name[cc])]   <-  tmppar$value[cc]
              ceattle_tmpSE[as.character(tmppar$name[cc])] <-  tmppar$std[cc]
              txt1  <- as.character(tmppar$name[cc])
            }
            
            datlist<-makeDat(tauIN=1,sigMethod=mode, 
                             estparams=estparams,typeIN=tt,
                             dataIN=rec,
                             covarsIN  = NULL)
            
            futR_mods[[s]][[tt]][[i]]<-fitRec(dlistIN=datlist,
                                              silentIN = silentINN,
                                              version="futR",
                                              recompile=F,
                                              simulate=simIT)
            
            npar           <-  length(which(!is.na(as.numeric(melt(futR_mods[[s]][[tt]][[i]]$model$env$map)$value))))
            
            if(rectypes[tt]=="LM")
              npar<-npar-1
            if(simIT){
              require(reshape)
              df <- melt(futR_mods[[s]][[tt]][[i]]$sim)
              colnames(df)<-c("parameter","itr","estimate")
              df$sp        <-  s
              df$model     <-  nm
              df$modelnum  <-  i
              df$CE        <-  TRUE
              df$npar      <-  npar
              df$n         <-  dim(rec)[1]
              df$rsType    <-  rectypes[tt]
              df$covset    <-  covset
              df$covtxt    <-  covtxt
            }
            
            
            tmppar2<- futR_mods[[s]][[tt]][[i]]$fit$par
            tmppar2<- futR_mods[[s]][[tt]][[i]]$pred
            futr_tmp <- data.frame(
              cov_mod   =  as.character(rownames(full.mat)[i]),
              sp        =  sp,
              recTypenm =  factor(rectypes[tt],levels=rectypes),
              recType   =  tt,
              model     =  "futR",
              i         =  i,
              log_a     =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="log_a",]$pred,
              log_b     =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="log_b",]$pred,
              logsigma  =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="logsigma",]$pred)
            futr_tmp <-cbind(futr_tmp,matrix(0,1,ncov))
            colnames(futr_tmp)[-(1:9)]<-paste0("rs_parm[",1:ncov,"]")
            
            futr_tmpSE <- data.frame(
              cov_mod   =  as.character(rownames(full.mat)[i]),
              sp        =  sp,
              recTypenm =  factor(rectypes[tt],levels=rectypes),
              recType   =  tt,
              model     =  "futR",
              i         =  i,
              log_a     =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="log_a",]$pred.se,
              log_b     =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="log_b",]$pred.se,
              logsigma  =  futR_mods[[s]][[tt]][[i]]$pred[futR_mods[[s]][[tt]][[i]]$pred$def=="logsigma",]$pred.se)
            futr_tmpSE <-cbind(futr_tmpSE,matrix(0,1,ncov))
            colnames(futr_tmpSE)[-(1:9)]<-paste0("rs_parm[",1:ncov,"]")
            
            
            if(any(!(paste0("rs_parm[",pnum,"]")%in%txt1)))
              message( "parm numbers don't match!")
            
            if(length(grep("rs_parm",futR_mods[[s]][[tt]][[i]]$pred$def))!=0){
              
              cc  <- grep("rs_parm",futR_mods[[s]][[tt]][[i]]$pred$def)
              futr_tmp[paste0("rs_parm[",pnum,"]")]    <-  futR_mods[[s]][[tt]][[i]]$pred[ cc,]$pred
              futr_tmpSE[paste0("rs_parm[",pnum,"]")]  <-  futR_mods[[s]][[tt]][[i]]$pred[ cc,]$pred.se
              
            }
            
            nm <- paste0("RS",s,"_",rectypes[tt])
            
            if(rectypes[tt]=="Ricker")
              nm <- paste0("RS",s,"_")
            
            # tmp1 <- data.frame(LL = LL, 
            #                    npar = npar, 
            #                    lab = 1:nn, 
            #                    R2 = R2)
            pred           <-  futR_mods[[s]][[tt]][[i]]$pred
            pred$sp        <-  s
            pred$model     <-  nm
            pred$modelnum  <-  i
            pred$CE        <-  TRUE
            pred$rsType    <-  rectypes[tt]
            pred$npar      <-  npar
            pred$n         <- dim(rec)[1]
            pred$covset    <-  covset
            pred$covtxt    <-  covtxt
            Rhat   <-  pred[pred$def=="R_hat",]
            Shat   <-  pred[pred$def=="S_hat",]
            datIN  <-  futR_mods[[s]][[tt]][[i]]$model$env$data
            lmm    <- lm(datIN$Robs~Rhat$pred)
            
            tmpaic<- 
              data.frame(
                LL    =   futR_mods[[s]][[tt]][[i]]$fit$objective,
                npar  =   npar,
                LnDet =   futR_mods[[s]][[tt]][[i]]$LnDet,
                n     =   dim(rec)[1],
                rsType =  rectypes[tt],
                name  =   nm,
                sp    =   s,
                CE    =   TRUE,
                conv  =   futR_mods[[s]][[tt]][[i]]$fit$convergence,
                R2    =   as.numeric(summary(lmm)["adj.r.squared"]),
                covset    =  covset,
                covtxt    =  covtxt)
            
           
              predOut  <-  rbind(predOut,pred)
              aicOut   <-  rbind(aicOut,tmpaic)
              if(simIT) dfOut    <-  rbind(dfOut,df)
              par_mat_ceattle    <- rbind(par_mat_ceattle,ceattle_tmp)
              par_mat_futR       <- rbind(par_mat_futR,futr_tmp)
              par_mat_ceattleSE  <- rbind(par_mat_ceattleSE,ceattle_tmpSE)
              par_mat_futRSE     <- rbind(par_mat_futRSE,futr_tmpSE)
            
          
          names(futR_mods[[s]][[tt]])[i]<-paste0("RS",s,"_",rectypes[tt])
        }
      }
      
  #____________________________________________
  # run aic analysis:
  #____________________________________________
    source("R/AIC_selection.R")
    AICtable_futR<-list()
  
    for(s in 1:3){
      aicsub <- as_tibble(aicOut)%>%filter(sp==s )
      AICtable_futR[[s]]<- AICselection(
        LL   = aicsub$LL, 
        npar = aicsub$npar, 
        n    = aicsub$n, 
        LnDet   = -aicsub$LnDet,
        mnames1 = aicsub$name, 
        covnm   = aicsub$covtxt, 
        rsType  = aicsub$rsType,
        CE      = aicsub$CE,
        R2      = aicsub$R2, 
        type2   = 2) 
      
    } 
  
    s<-1
    head(AICtable[[s]])
    head(AICtable_futR[[s]])
  
    head(AICtable[[s]][-grep("LM",AICtable[[s]]$names),])
    head(AICtable_futR[[s]][-grep("LM",AICtable_futR[[s]]$names),])

  #___________
  # model average the parm results:
  #
  PARMAT<-as_tibble(par_mat_futR)
  PARMAT%>%filter(recTypenm=="Ricker")
  # do it better for top AIC models using simIT = True
    
  
  #____________________________________________
  # confirm that input parms yield the same answer for both models:
  #____________________________________________
  i  <-  157
  s  <-  1

  # Ricker

  tt  <-  4
  m   <-  futR_mods[[s]][[tt]][i]
  m[[1]]$model$env$par
      
  # BLM
  
  # LM
  
  # BevHolt


  #______________
  # create projections





