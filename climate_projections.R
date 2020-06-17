# Kirstin Holsman
# kirstin.holsman@noaa.gov
# futR() for use with CEATTLE (holsman et al. 2016,2020) and Rceattle() (adams et al. in prep)
# based on the Climate-Enhanced, Age-based model with Temperature-specific 
# Trophic Linkages and Energetics (CEATTLE)
# model epsi_s for strategic management advice in the EBS 

# 1. Load data
# 2. run estimation mode (or load est models)
# 3. Compile the model
# 4. Run recruitment estimates
# 5. Project under HCRS
 
  # rm(list=ls())
   
  #___________________________________________
  # 1. Set things up
  #___________________________________________
  
  # 1.0. define directories
    main<-path.expand("~/GitHub/futR")
    setwd(main)
     
  # 1.1. load libraries
    library( Rceattle )
    library( TMB )
    library(reshape)
    library(dplyr)
    library(ggplot2)
    library(lattice)
    
  # 1.2 input switches
    run_est    <- FALSE  # run the ceattle estimation model?
    REcompile  <- TRUE   # recompile futR 
    bias_corr  <- TRUE   # bias correct sigma sensu Ludwig and Walters 1981 and Porch and Lauretta 2016
    recruitAge <- 1      # what age (lag) is recruitment estimated in the model?
    fityrsIN   <- 1980:2017  # what years of the recruiment data do you want to fit too (allows for leave one out analysis)
    Eat_covIN  <- c(18, 20, 26, 28)   # specific to CEATTLE
    Eat_covIN  <- -99    # specific to CEATTLE; nulled here using -99
    
  # 1.3 load data
    load("data/env_covars.Rdata")
    load("data/rec_dat.Rdata")
    load("data/rec_dat_tmb.Rdata")
    load("data/ration_tmb.Rdata")
  
  # 1.4 load functions
    source("R/futRfun.R")
  
  #___________________________________________
  # 2. run Rceattle() estimation mode
  #___________________________________________

    if( run_est ){   
      # Rceattle based on Grant Adams TMB formulation of the EBS CEATTLE model
        data(BS2017SS)
        data(BS2017MS)

        # 2.1 Run single species model est
        ss_run <- Rceattle(
          data_list    = BS2017SS,
          inits        = NULL, # Initial parameters = 0
          file_name    = NULL, # Don't save
          debug        = 0, # Estimate
          random_rec   = FALSE, # No random recruitment
          msmMode      = 0, # Single species mode
          avgnMode     = 0,
          silent       = TRUE
          )

        # 2.2 Run multispecies model est
        ms_run <- Rceattle(
          data_list   = BS2017MS,
          inits       = ss_run$estimated_params, # Initial parameters from ss run
          file_name   = NULL, # Don't save
          debug       = 0, # Estimate
          random_rec  = FALSE, # No random recruitment
          niter       = 10, # Number of iterations around predation/pop dy functions
          msmMode     = 1, # Multi-species holsman mode
          avgnMode    = 0 # Use average N
        )

        save( ss_run,file   =   "data/ss_run.Rdata" )
        save( ms_run,file   =   "data/ms_run.Rdata" )
    }else{

        load( "data/ss_run.Rdata" )
        load( "data/ms_run.Rdata" )

    }

  #___________________________________________
  # 3. Compile futR	
  #___________________________________________
  
  compile('futR.cpp') # this will generate warnings - they can be ignored if "0" is returned

 
  # 3.1 Set up model data list and map and covars for epsi_s 
    # epsi_s environmental variabile from the age class year (assuming Jan rec, and summer = age0)
  
    PAR<-data.frame(
      phases = c(
          log_a        = 1, 
          log_b        = 1 ,
          #logit_tau    = 1,
          rs_parm      = 1,
          epsi_s       = 1,
          logsigma     = 1),
      estparams  = c(
          log_a        = TRUE, 
          log_b        = TRUE, 
          #logit_tau    = FALSE,
          rs_parm      = TRUE,
          epsi_s       = FALSE,
          logsigma     = TRUE))
  
   
  
      # For Kir : skip
      if(1==10){
          load("data/aclim_cmip5_BT.Rdata")  # Bottom Temp object called "allDat"
          load("data/NPZcov.Rdata")
          load("data/smry_RS.Rdata")
          topmod  <-  colnames(NPZcov)[-1]
          subCOV1 <-  NPZcov[NPZcov$year%in%(yearsIN-1),]  # off set by 1 year Rec in year y ~f(covars in y-1)
          subCOV2 <-  allDat[allDat$t%in%(yearsIN-1),]     # off set by 1 year Rec in year y ~f(covars in y-1)
          
          covars                <-  subCOV1[,-1]
          rownames(covars)      <-  subCOV1[,1]
          tmp                   <-  as.matrix((t(covars)))
          sub_cov               <-  matrix(0,dim(tmp)[1],dim(tmp)[2])  #empty matrix
          
          for(rr in 1:dim(tmp)[1])
            sub_cov[rr,]        <-   (tmp[rr,])
          colnames(sub_cov)     <-   subCOV1[,1]
          nm                    <-   rownames(tmp)
            env_covars            <-   sub_cov[c(1,2),]
        
        # set up recruitment input data file Years, SSB, Robs, sdSSB, sdRobs
            rec_dat  <-  rec_dat_tmb  <-  list()
            for(sp in 1:3){
                est        <-  ms_run
                smry_RSIN  <-  smry_RS[[sp]]
                sdlist     <-  as.list(est$sdrep,"Std. Error")
                estlist    <-  as.list(est$sdreNAp,"Estimate")
                mpRsd      <-  sdlist$rec_dev[sp,-1]
                
                dataIN_tmb <-data.frame(
                        years  = est$data_list$styr-1+(1:(est$data_list$nyrs))[-1],
                        SSB    = est$quantities$biomassSSB[sp,-length(est$quantities$biomassSSB[sp,])],
                        Robs   = est$quantities$R[sp,-1],
                        sdSSB  = 0,
                        sdRobs = mpRsd)
          
          
                dataIN_admb<-data.frame(
                   years  = smry_RSIN$Rec_years,
                   SSB    = smry_RSIN$SSB.y.1.,
                   Robs   = smry_RSIN$R_obs,
                   sdSSB  = 0,
                   sdRobs = mpRsd) 
          
                rec_dat[[sp]]     <-  dataIN_admb
                rec_dat_tmb[[sp]] <-  dataIN_tmb
                
                if(sp==1)
                  ration_tmb<- data.frame(sp1=as.numeric(scale(apply(est$quantities$ration2Age[-sp,-1,],3,sum)))[-length(smry_RSIN$Rec_years)],sp2=NA,sp3=NA)
                if(sp>1)
                  ration_tmb[,sp]    <-  (as.numeric(scale(apply(est$quantities$ration2Age[-sp,-1,],3,sum)))[-length(smry_RSIN$Rec_years)])
              }
          
              save(env_covars,file =  "data/env_covars.Rdata")
              save(rec_dat,file    =  "data/rec_dat.Rdata")
              save(rec_dat_tmb,file    =  "data/rec_dat_tmb.Rdata")
              save(ration_tmb,file =  "data/ration_tmb.Rdata")
      }
  
  
    rec        <-  rec_dat[[1]]
    env        <-  env_covars
    env[1,]    <-  as.numeric(scale(env_covars[1,]))
    env[2,]    <-  as.numeric(scale(env_covars[2,]))
    ration     <-  ration_tmb[,1]
    version    <-  "futR"  # which TMB model to use

    tau_set<-c(0,.01,.1,0.5,.8,1,1.2,1.5,2,5,10)
  
  # 3.2 Set up data
  
    PAR$phases
    PAR$estparams
    
    # which parameters to estimate with futR?
    estparams  = c(
          log_a        = TRUE, 
          log_b        = TRUE, 
          #logit_tau     = TRUE,
          rs_parm      = TRUE,
          epsi_s       = FALSE,
          logsigma     = TRUE)
  
    rec_noerr<-rec
    rec_noerr$sdRobs<-0

    #___________________________________________
    # 4. RUN FUTR()
    #___________________________________________
    
    #__________________________________________________________________
    # rectype/ typeIN:
        # 1 = Linear with biomass ( y-1 )  BLM
        # 2 = Linear LM
        # 3 = Beverton Holt 
        # 4 = ricker  
    
    # sigMethod:
        # 1 = estimate simga
        # 2 = unbiased sigma empirically calculated (rather than biased estimate from MLE; sensu Ludwig and Walters):
        # 3 = as in 1 but with specified measurement error for rec (indep of random effects on S)
        # 4 = as in 1 but with specified measurement error for rec and SSB    
    #__________________________________________________________________

    
    # RICKER MODELS (typeIN = 4)
    
    # 4.1 First run with no observation error on SSB:
         # makeDat will make the input values, data, and phases for the model:
          datlist  <-  makeDat(
                    tauIN      =  1,
                    sigMethod  =  1, 
                    estparams  =  estparams,
                    typeIN     =  4,
                    dataIN     =  rec,
                    covarsIN   =  env)
          
          # run the model
          mm1                 <-  runmod(dlistIN=datlist,version="futR",recompile=T,simulate=TRUE)
          df1                 <-  data.frame(estimate=as.vector(mm1$sim), parameter=names( mm1$mle)[row(mm1$sim)])
          densityplot( ~ estimate | parameter, data=df1, layout=c(5,1),ylim=c(0,10))
          
          #now change tau
          datlist$rs_dat$tau  <-  0.000001
          
          # re-run the model with tau 
          mm1_t0              <-  runmod(dlistIN=datlist,version="futR",recompile=F,simulate=TRUE)
          df1_t0              <-  data.frame(estimate=as.vector(mm1_t0$sim), parameter=names( mm1_t0$mle)[row(mm1_t0$sim)])
          densityplot( ~ estimate | parameter, data=df1_t0, layout=c(5,1),ylim=c(0,10))
     
          
          
          
    # 4.2 Now run the model with unbiased sigma (sigMethod = 2):
          datlist  <-  makeDat(
                    tauIN      =  1,
                    sigMethod  =  2, 
                    estparams  =  estparams,
                    typeIN     =  4,
                    dataIN     =  rec,
                    covarsIN   = env)
         
          # run the model
          mm2                 <-  runmod(dlistIN=datlist,version="futR",recompile=F,simulate=TRUE)
          df2                 <-  data.frame(estimate=as.vector(mm2$sim), parameter=names( mm2$mle)[row(mm2$sim)])
          densityplot( ~ estimate | parameter, data=df2, layout=c(5,1),ylim=c(0,10))
          
          #now change tau
          datlist$rs_dat$tau<-0.000001
          
          # re-run the model with tau 
          mm2_t0              <-  runmod(dlistIN=datlist,version="futR",recompile=F,simulate=TRUE)
          df2_t0              <-  data.frame(estimate=as.vector(mm2_t0$sim), parameter=names( mm2_t0$mle)[row(mm2_t0$sim)])
          densityplot( ~ estimate | parameter, data=df2_t0, layout=c(5,1),ylim=c(0,10))
          
          
          
          
    # 4.3 Now run the model estimating sigma; specify sdR (random effects on SSB, random effects on Rec = sig+sdR)
          datlist   <-  makeDat(
                      tauIN      =  1,
                      sigMethod  =  3, 
                      estparams  =  estparams,
                      typeIN     =  4,
                      dataIN     =  rec,
                      covarsIN   =  env)
          
          # run the model
          mm3                 <-  runmod(dlistIN=datlist,version="futR",recompile=F,simulate=TRUE)
          df3                 <-  data.frame(estimate=as.vector(mm3$sim), parameter=names( mm3$mle)[row(mm3$sim)])
          densityplot( ~ estimate | parameter, data=df3, layout=c(5,1),ylim=c(0,10))
          
          #now change tau
          datlist$rs_dat$tau<-0.000001
          
          # re-run the model with tau 
          mm3_t0              <-  runmod(dlistIN=datlist,version="futR",recompile=F,simulate=TRUE)
          df3_t0              <-  data.frame(estimate=as.vector(mm3_t0$sim), parameter=names( mm3_t0$mle)[row(mm3_t0$sim)])
          densityplot( ~ estimate | parameter, data=df3_t0, layout=c(5,1),ylim=c(0,10))
          
    # 4.4 Now run the model estimating sigma; specify sdR (random effects on SSB, random effects on Rec = sig+sdR)
          datlist   <-  makeDat(
                    tauIN        =  1,
                    sigMethod    =  4, 
                    estparams    =  estparams,
                    typeIN       =  4,
                    dataIN       =  rec,
                    covarsIN     =  env)
        
          # run the model
          mm4                 <-  runmod(dlistIN=datlist,version="futR",recompile=F,simulate=TRUE)
          df4                 <-  data.frame(estimate=as.vector(mm4$sim), parameter=names( mm4$mle)[row(mm4$sim)])
          densityplot( ~ estimate | parameter, data=df4, layout=c(5,1),ylim=c(0,10))
          
          #now change tau
          datlist$rs_dat$tau<-0.000001
          
          # re-run the model with tau 
          mm4_t0              <-  runmod(dlistIN=datlist,version="futR",recompile=F,simulate=TRUE)
          df4_t0              <-  data.frame(estimate=as.vector(mm4_t0$sim), parameter=names( mm4_t0$mle)[row(mm4_t0$sim)])
          densityplot( ~ estimate | parameter, data=df4_t0, layout=c(5,1),ylim=c(0,10))
    
          
        mle_set<- data.frame(rbind(
                  mm1$mle,
                  mm1_t0$mle,
                  c(mm2$mle,NA),
                  c(mm2_t0$mle,NA),
                  mm3$mle,
                  mm3_t0$mle,
                  mm4$mle,
                  mm4_t0$mle))
      #_______________________________________________
      ## NOW FOR EACH SPP:
      #_______________________________________________

        # 4.5 Estimate recruitment
        tmplist<-list()
        env     <-  env_covars
          env[1,]    <-  as.numeric(scale(env_covars[1,]))
          env[2,]    <-  as.numeric(scale(env_covars[2,]))
        ration  <-  ration_tmb[,1]
        version <-  "futR"

        tau_set<-c(0,.01,.1,0.5,.8,1,1.2,1.5,2,5,10)
      
        # 4.2 Set up data and fit to env covariates
        PAR$phases
        PAR$estparams
        estparams  = c(
              log_a        = TRUE, 
              log_b        = TRUE, 
              #logit_tau     = TRUE,
              rs_parm      = TRUE,
              epsi_s       = FALSE,
              logsigma     = TRUE)
      
        noSDR   <- FALSE
        mod<-list()
        sp<-1

        for(sp in 1:3){ 
            rec     <-  rec_dat[[sp]]
            if(noSDR) rec$sdRobs<-0
            tmpp<-list()
        
            datlist  <  -makeDat(tauIN=1,sigMethod=1, estparams=estparams,typeIN=4,dataIN=rec,covarsIN  = env)
              mm1<-runmod(dlistIN=datlist,version="futR",recompile=T)
              tmpp$mm1<-mm1;tmpp$mm1$datlist<-datlist
        
              datlist$rs_dat$tau<-0.000001
              mm1_t0<-runmod(dlistIN=datlist,version="futR",recompile=F)
              tmpp$mm1_t0<-mm1_t0;tmpp$mm1_t0$datlist<-datlist
        
            # unbiased sigma
              # tauIN = 1; estparms$epsi_s = FALSE; estparms$sigma = TRUE; rciker
            datlist<-makeDat(tauIN=1,sigMethod=2, estparams=estparams,typeIN=4,dataIN=rec,covarsIN  = env)
              mm2<-runmod(dlistIN=datlist,version="futR",recompile=F)
              tmpp$mm2<-mm2;tmpp$mm2$datlist<-datlist
            datlist$rs_dat$tau<-0.000001
              mm2_t0<-runmod(dlistIN=datlist,version="futR",recompile=F)
               tmpp$mm2_t0<-mm2_t0;tmpp$mm2_t0$datlist<-datlist
        
            # est sigma; specify sdR (random effects on SSB, random effects on Rec = sig+sdR)
              # tauIN = 1; estparms$epsi_s = FALSE; estparms$sigma = TRUE; rciker
               estparams$epsi_s<-TRUE
            datlist<-makeDat(tauIN=1,sigMethod=3, estparams=estparams,typeIN=4,dataIN=rec,covarsIN  = env)
              mm3<-runmod(dlistIN=datlist,version="futR",recompile=F)
               tmpp$mm3<-mm3;tmpp$mm3$datlist<-datlist
        
            # est sigma; specify sdR (random effects on SSB, random effects on Rec = sig+sdR)
              # tauIN = 1; estparms$epsi_s = FALSE; estparms$sigma = TRUE; rciker
            datlist<-makeDat(tauIN=.5,sigMethod=4, estparams=estparams,
              typeIN=4,dataIN=rec,covarsIN  = env)
              mm4<-runmod(dlistIN=datlist,version="futR",recompile=F)
              tmpp$mm4<-mm4;tmpp$mm4$datlist<-datlist
        
        
                mod[[sp]]<-tmpp
                rm(tmpp);rm(mm1);rm(mm1_t0);rm(mm2);rm(mm2_t0);rm(mm3);rm(mm4)
          cat(paste("----------------------\nspecies ",sp, "complete\n----------------------"))
        }


       # 4.6 (plotit)
        par(mfrow=c(3,1))
        mode<-"mm2"
        mode<-1
        xlimm<-c(9e6,6e5,5e5)
        ylimm<-c(1.2e8,2e6,6e5)
        cexx<-.5
        units<-c(1e5,1e4,1e4)
  

        for(mode in 1:length(mod[[1]])){
            for(recIN in 1:2){
                 if(recIN==1) {
                   rec_datUSE  <- rec_dat
                   model      <- "ADMB"
                 }else{
                   rec_datUSE  <- rec_dat_tmb
                   model      <- "TMB"
                 }
                for(sp in 1:3){
                   if(sp==1&recIN==1&mode==1){
                     rec  <-  data.frame(rec_datUSE[[sp]],sp=sp,mode=names(mod[[sp]])[mode],mod=model)
                   }else{
                     rec  <-  rbind(rec,data.frame(rec_datUSE[[sp]],sp=sp,mode=names(mod[[sp]])[mode],mod=model))
                   }
                   final  <-  mod[[sp]][[mode]]
                   rc     <-  final$datlist$rs_dat$rs_cov
                   
                   # 4.4 predict Rec from new data - plot RS relationship under +1 and -1 covars
                   Snew         <-  list(SSB=seq(1,max(rec_datUSE[[sp]]$SSB)*2,units[sp]))
                   Snew$rs_cov  <-  matrix(1,length(rc[,1]),length(Snew[[1]]))
                   nd1          <-  predict_futR(object =  final,newdata=Snew,se.fit=TRUE, simulate=TRUE,nitr=1000)
                   
                   Snew$rs_cov<-matrix(-1,length(rc[,1]),length(Snew[[1]]))
                   nd     <-  predict_futR(object =  final,newdata=Snew,se.fit=TRUE, simulate=TRUE,nitr=1000)
                   dlist  <-  data.frame(years=final$datlist$rs_dat$years,Robs=final$datlist$rs_dat$Robs,SSB=final$datlist$rs_dat$SSB)
                   
                   dat    <-  nd$sim
                   out0   <-  data.frame(
                     melt(apply(dat[grep("R_hat",rownames(dat)),],1,quantile,probs=c(.95,.5,.05))),
                     melt(apply(dat[grep("S_hat",rownames(dat)),],1,quantile,probs=c(.95,.5,.05))),
                     sp = sp,type="-1SD",mode=names(mod[[sp]])[mode],mod=model,SSB=rep(Snew$SSB,each=3))
                   
                   dat    <-  nd1$sim
                   out1   <-  data.frame(
                     melt(apply(dat[grep("R_hat",rownames(dat)),],1,quantile,probs=c(.95,.5,.05))),
                     melt(apply(dat[grep("S_hat",rownames(dat)),],1,quantile,probs=c(.95,.5,.05))),
                     sp = sp,type="+1SD",mode=names(mod[[sp]])[mode],mod=model,SSB=rep(Snew$SSB,each=3))
                   if(sp==1&recIN==1&mode==1){
                     out  <-  rbind(out0,out1)
                   } else{
                     out  <-  rbind(out,rbind(out0,out1))
                   }
                }
            }  
        }

        rec$Rec  <- rec$Robs
        rec$mod  <- factor(rec$mod,levels=c("ADMB","TMB"))
        out$mod  <- factor(out$mod,levels=c("ADMB","TMB"))
        rec$mode <- factor(rec$mode,levels=names(mod[[sp]]))
        out$mode <- factor(out$mode,levels=names(mod[[sp]]))
        out$type <- factor(out$type,levels=c("-1SD","+1SD"))
        
        if(!dir.exists("Figs")) dir.create("Figs")
        source("R/plot_rec.R")
        
        graphics.off()
        HH<-7;WW<-6
        dev.new(height=HH,width=WW)
        for(modeIN in 1:length(mod[[1]])){
          dev.new(height=HH,width=WW)
          plotCompareMod(H=HH,W=WW,mode=modeIN,outIN=out,recIN=rec)
          ggsave(file=paste0("Figs/compare_",names(mod[[sp]])[modeIN],".jpg"),device="jpg",width=WW,height=HH)
        }
  
        
      #___________________________________________
      # 5. Project forward (under development)
      #___________________________________________
      
      # 5.1 Run model   # recall that parameters mepsi_s in the same order
    
      # project forward epsi_snew data and HCR
      for(sp in 1:3){
        for(mode in 1:length(mod[[sp]])){
          pred   <-  mod[[1]][[1]]$pred
          Rhat   <-  pred[pred$def=="R_hat",]
          Shat   <-  pred[pred$def=="S_hat",]
          datIN  <-  mod[[1]][[1]]$datlist$inputs$dataIN
          lmm    <-  lm(datIN$Robs~Rhat$pred)
          
          data.frame(
            LL    =   mod[[1]][[1]]$fit$objective,
            npar  =   length(mod[[1]][[1]]$fit$par),
            model =   factor(names(mod[[1]])[mode],levels=names(mod[[1]])),
            conv  =   mod[[1]][[1]]$fit$convergence,
            R2    =   as.numeric(summary(lmm)["adj.r.squared"]))
          
         # predict_futR()
        }
      }

        
        

