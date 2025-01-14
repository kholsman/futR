rm(list=ls()); setwd("D:/GitHub_cloud/futR")
#___________________________________________
# 1. Set things up
#___________________________________________
# rm(list=ls()) ; dir()

# load data, packages, setup, etc.
source("R/01_make.R")

#___________________________________________
# 2. Compile futR	(first time through - can skip this step after )
#___________________________________________

recompile_model <- FALSE # to recompile the model set to TRUE

if(recompile_model){
  wd0 <- getwd()
  setwd("src")
  recompile('futR')
  setwd(wd0)
}
# this will generate warnings - they can be ignored if "0" is returned

# read in the data and create a datlist
datlist <- makefutR_data("data/in/futR_Inputs.xlsx" )

# recruitment data:
datlist$rs_dat$R_obs
datlist$rs_dat$S_obs

# covar data:
datlist$rs_dat$rs_cov
# rec        <-  rec_dat[[1]]
# env        <-  env_covars

# which parameters to estimate with futR?
datlist$estparams

# starting values? 
datlist$parameters

# parameter map:
datlist$maplist

# which phases to estimate in (not yet coded up)
datlist$phases


# set some global values for the demo below:
estparams <-  datlist$estparams[1:6]
rec       <-  data.frame(readxl::read_xlsx("data/in/futR_Inputs.xlsx" , sheet = "rec_data" ))



datlist2  <-  makeDat(
  rectype    =  4,
  tauIN      =  0,
  sigMethod  =  1, # (default, no random effects)
  estparams  =  estparams,
  estMode    =  1,
  rec_years  =  rec$rec_year,
  Rec        =  rec$Robs,
  SSB        =  rec$SSB,
  sdSSB      =  rec$sdSSB,
  sdRec      =  rec$sdRobs,
  covars     =  NULL,
  covars_sd  =  NULL)

datlist2$estparams["beta"] <- datlist$estparams["lambda"] <- TRUE

Rec4 <-  mm <-runmod(dlistIN   = datlist2,
                     version   = 'futR',
                     recompile = FALSE,
                     simulate  = TRUE,
                     sim_nitr  = 1000)

# summarize results
dfR4    <-  data.frame(model = "Rec 4",
                       estimate  = as.vector(mm$sim),
                       parameter = names( mm$mle)[row(mm$sim)])
r4_fit  <- getFit(mm, nm = "recType = 4")

rm(mm)


#-----------------------------------------------
# Estimate covariate effects:
#-----------------------------------------------
# read in the data and create a datlist (rather than hand code it)
datlist <- makefutR_data("data/in/futR_Inputs.xlsx" )

datlist$rs_dat$rs_cov
datlist$rs_dat$sdrs_cov

estparams3 <- estparams
estparams3["beta"] <- estparams3["lambda"] <- TRUE

datlist3  <-  makeDat(
  rectype    =  4,
  tauIN      =  0,
  sigMethod  =  1, # (default, no random effects)
  estparams  =  estparams3,
  estMode    =  1,
  rec_years  =  rec$rec_year,
  Rec        =  rec$Robs,
  SSB        =  rec$SSB,
  sdSSB      =  rec$sdSSB,
  sdRec      =  rec$sdRobs,
  beta_0     =  datlist$rs_dat$beta_0,
  lambda_0   =  datlist$rs_dat$lambda_0,
  covars     =  datlist$rs_dat$rs_cov,
  covars_sd  =  datlist$rs_dat$sdrs_cov)

# set the rectype to 4
datlist3$rs_dat$rectype <- 4


# run the basic model
Rec4_covar <-  mm <-runmod(dlistIN   = datlist,
                           version   = 'futR',
                           recompile = FALSE,
                           simulate  = TRUE,
                           sim_nitr  = 1000)

# summarize results
dfR4_c    <-  data.frame(model = "Rec 4 cov2",
                         estimate  = as.vector(mm$sim),
                         parameter = names( mm$mle)[row(mm$sim)])
df      <- rbind(dfR4,dfR4_c)
r4_fitc  <- getFit(mm, nm = "recType = 4 cov2")
rec_fit <- rbind(r4_fit,r4_fitc)
rm(mm)

#plot_par_pdf(df)
plot_rs(rec_fit)

#-----------------------------------------------
# compare to just Coldpool (second row)
#-----------------------------------------------

# Estimate covariate effects:
estparams3["beta"] <- estparams3["lambda"] <- TRUE

datlist4  <-  makeDat(
  rectype    =  4,
  tauIN      =  0,
  sigMethod  =  1, # (default, no random effects)
  estparams  =  estparams3,
  estMode    =  1,
  rec_years  =  rec$rec_year,
  Rec        =  rec$Robs,
  SSB        =  rec$SSB,
  sdSSB      =  rec$sdSSB,
  sdRec      =  rec$sdRobs,
  beta_0     =  datlist$rs_dat$beta_0,
  lambda_0   =  datlist$rs_dat$lambda_0,
  covars     =  matrix(datlist$rs_dat$rs_cov[2,],1,dim(datlist$rs_dat$rs_cov)[2]),
  covars_sd  =  matrix(datlist$rs_dat$rs_cov[2,],1,dim(datlist$rs_dat$sdrs_cov)[2]))

# set the rectype to 4
datlist4$rs_dat$rectype <- 4


# run the basic model
Rec5_covar <-  mm <-runmod(dlistIN   = datlist,
                           version   = 'futR',
                           recompile = FALSE,
                           simulate  = TRUE,
                           sim_nitr  = 1000)

# summarize results
dfR5_c    <-  data.frame(model = "Rec 5 cov2",
                         estimate  = as.vector(mm$sim),
                         parameter = names( mm$mle)[row(mm$sim)])
df       <- rbind(df,dfR5_c)
r5_fitc  <- getFit(mm, nm = "recType = 5 cov2")
rec_fit  <- rbind(rec_fit,r5_fitc)
rm(mm)

#plot_par_pdf(df)
plot_rs(rec_fit)


