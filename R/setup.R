# ----------------------------------------
# setup.R
# subset of Holsman et al. 2020 Nature Comm.
# kirstin.holsman@noaa.gov
# updated 2020
# ----------------------------------------
  
    # set up directory paths:
    #-------------------------------------------
    fun_dir         <- file.path(main,"R/sub_Fun")
    in_dir          <- file.path(main,"data/in")
    out_dir         <- file.path(main,"data/out")  
    
    # switches and options:
    #-------------------------------------------
    version         <-  "src/futR"  # which TMB model to use
    update.figs     <-  FALSE   # set to true to re-save figs
    update.inputs   <-  FALSE   
    update.outputs  <-  FALSE   # overwrite the existing Rdatafiles in data/out
    status          <-  TRUE   # print progress
    scaleIN         <-  1      # controls the ratio (relative scaling of window)
    dpiIN           <-  600    # dpi for figures (set to lower res for smaller file size- these will be about 3.5 MB)
    
    run_est    <- FALSE  # run the ceattle estimation model?
    REcompile  <- TRUE   # recompile futR 
    bias_corr  <- TRUE   # bias correct sigma sensu Ludwig and Walters 1981 and Porch and Lauretta 2016
    recruitAge <- 1      # what age (lag) is recruitment estimated in the model?
    fityrsIN   <- 1980:2017  # what years of the recruiment data do you want to fit too (allows for leave one out analysis)
    Eat_covIN  <- c(18, 20, 26, 28)   # specific to CEATTLE
    Eat_covIN  <- -99    # specific to CEATTLE; nulled here using -99
    
    
    # Set up model data list and map and covars for epsi_s 
    #-------------------------------------------
    # epsi_s environmental variabile from the age class year (assuming Jan rec, and summer = age0)
    
    PAR <-data.frame(
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
      
      save(env_covars,file     =  file.path(in_dir,"env_covars.Rdata"))
      save(rec_dat,file        =  file.path(in_dir,"rec_dat.Rdata"))
      save(rec_dat_tmb,file    =  file.path(in_dir,"rec_dat_tmb.Rdata"))
      save(ration_tmb,file     =  file.path(in_dir,"ration_tmb.Rdata"))
    }
    
    
    rec        <-  rec_dat[[1]]
    env        <-  env_covars
    env[1,]    <-  as.numeric(scale(env_covars[1,]))
    env[2,]    <-  as.numeric(scale(env_covars[2,]))
    ration     <-  ration_tmb[,1]
    
    
    tau_set<-c(0,.01,.1,0.5,.8,1,1.2,1.5,2,5,10)
    
    # 3.2 Set up data
    
    PAR$phases
    PAR$estparams
    
    # which parameters to estimate with futR?
    estparams  = c(
      log_a        = TRUE, 
      log_b        = TRUE, 
      #logit_tau     = TRUE,
      beta         = FALSE,
      lambda       = TRUE,
      epsi_s       = FALSE,
      logsigma     = TRUE)
    
    rec_noerr<-rec
    rec_noerr$sdRobs<-0

   # Species stuff: (used for plotting and manipulating data)
   #-------------------------------------------
    sppINFO<-list(
    plk=list(abv="plk",
               guildIN="Walleye pollock",
               plotSPP="walleye pollock",
               bin2=c(seq(0,300,10),1000),
               binJvAD=c(0,40,1000),
               splistIN="W. Pollock",doNEBS=T,plotIT=T),
      pcod=list(abv="pcod",
                guildIN="Pacific cod",
                plotSPP="Pacific cod",
                bin2=c(seq(0,300,10),1000),
                binJvAD=c(0,40,1000),
                splistIN="P. Cod",doNEBS=T,plotIT=T),
      atf=list(abv="atf",
               guildIN="Arrowtooth or Kamchatka",
               plotSPP="arrowtooth flounder",
               bin2=c(seq(0,300,10),1000),
               binJvAD=c(0,40,1000),
               splistIN=c("Arrowtooth","Arrow or Kam", "Kamchat fl"),doNEBS=F,plotIT=T)
    )
    
    # These switches for KHolsman during simulation updates:
    #-------------------------------------------
    retroFL         <-  "data/in/raw/retro_data2018_long_ext_bcs.dat"
    futFL           <-  "data/in/raw/proj_data2018_long_ext_bcs.dat"
    fldr_nm         <-  "aclim_00_JunV2_2019"  # folder with the CEATTLE assessment runs
    UpdateMCMC      <-  1      # update MCMC? 1 = TRUE, 0 = FALSE
    readdat         <-  FALSE  # re-read in new data?
    update.simlist  <-  FALSE  # only TRUE when re-running CEATTLE simulations 
    update.romsnpz  <-  FALSE  # only TRUE when re-running CEATTLE simulations
    
    
    # Plotting stuff:
    #-------------------------------------------
    # The width of figures, when printed, 
    # will usually be 5.5 cm (2.25 inches or 1 column) 
    # or 12.0 cm (4.75 inches or 2 columns). 

    # set up color palettes
    plt     <- c("Zissou1","Darjeeling1","Darjeeling2","FantasticFox1")
    blues   <- RColorBrewer::brewer.pal(5, "Blues")
    BG      <- RColorBrewer::brewer.pal(9, "GnBu")  #5
    Ornjazz <- RColorBrewer::brewer.pal(5, "Oranges")
    YGB     <- (RColorBrewer::brewer.pal(5, "YlGnBu"))
    bg      <- colorRampPalette(BG)
    YlGnBu  <- colorRampPalette(YGB[-1])
    blu     <- colorRampPalette(blues[-1])
    night   <- colorRampPalette(colors()[c(653,47,474,72,491,477)])
    dawn    <- colorRampPalette(c(colors()[c(477,491,72,474,47,653)],"orange","red"))
    orng    <- colorRampPalette(Ornjazz[1:5])
    colIN1   <- colorRampPalette(c(wes_palette(n=5, name=plt[1])[1:5]))
    col4     <- colorRampPalette(c(colIN1(6),"maroon"))
    
    col_in  <- colorRampPalette(colors()[c(459,122,73)])
    col_in  <- colorRampPalette(colors()[c(408,44,73)])
    col_in2 <- colorRampPalette(c("orange","red"))
    wes     <- colorRampPalette(c(wes_palette(n=5, name=plt[1])[1:5]))
    col1    <- colorRampPalette(colors()[c(280,320)])
    #col2   <- colorRampPalette(colors()[c(70,491)])
    col2    <- colorRampPalette(colors()[c(114,491)])
    col2    <- colorRampPalette(c(wes(7)[c(3,1)],col2(3)[3]))
    col3    <- colorRampPalette(c(wes(7)[4:7]))
    
    # set the color scheme
    coll_use         <-  c(colors()[320],col2(8)[c(2,5,7)],col3(8)[c(2,5,7)])
    
    # Set up plotting stuff:
    #------------------------------------     
   
    probbs        <-  c(.1,.25,.5,.75,.9)
    alphaAll      <-  ceiling(rep(100/(length(probbs)/2),length(probbs)/2))
    c1            <-  col2(5)
    c2            <-  col2(6)[seq(2,6,2)]
    c3            <-  col3(6)[seq(2,6,2)]
    collIn        <-  rep(NA,13)
    collIn[1:2]   <-  col2(2)[2]
    ltyall                  <-  rep(1,13)  
    ltyall[1:2]             <-  1
    lwdall                  <-  rep(1,13)  
    lwdall[1:2]             <-  2
 

