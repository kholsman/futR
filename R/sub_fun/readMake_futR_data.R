#'
#'
#'
#'read_futR_data.R
#'
#'

readMake_futR_data <- function(fn= "data/in/futR_Inputs.xlsx"){
  
  tmp            <- data.frame(readxl::read_xlsx(fn, sheet = "rec_data" ))
  ec             <- data.frame(readxl::read_xlsx(fn, sheet = "covar_val" ,skip=3)) 
  ec_sd          <- data.frame(readxl::read_xlsx(fn, sheet = "covar_val_sd")) 
  beta_0IN       <- data.frame(readxl::read_xlsx(fn, sheet = "covar_val", n_max =2)[1,-1])
  lambda_0IN     <- data.frame(readxl::read_xlsx(fn, sheet = "covar_val", n_max =2 )[2,-1])
  estparams      <- data.frame(readxl::read_xlsx(fn, sheet = "estparams" ))
  switches       <- tibble::deframe(readxl::read_xlsx(fn, sheet = "switches"))
  rec_dat        <- tmp[,1:5]
  fityears       <- tmp$rec_year[tmp$fit==1]
  
  ec             <- ec[ec$year%in%rec_dat$rec_year,]
  tt             <- ec[,-1]
  rownames(tt)   <- as.character(ec[,1])
  ec <- t(tt)
  ec_sd          <- ec_sd[ec_sd$year%in%rec_dat$rec_year,]
  tt             <- ec_sd[,-1]
  rownames(tt)   <- as.character(ec_sd[,1])
  ec_sd <- t(tt)
  
  datlist  <-  makeDat(
    rectype    =  as.numeric(switches["rectype"]),
    tauIN      =  as.numeric(switches["tauIN"]),
    sigMethod  =  as.numeric(switches["sigMethod"]), # (default, no random effects)
    estparams  =  tibble::deframe(estparams[,c(1,3)]),
    estMode    =  1,
    rec_years  =  rec_dat$rec_year,  # recruitment year
    fityrs     =  fityears ,         # fit the model to the whole set of rec_years
    Rec        =  rec_dat$Robs,      # recruitment index
    SSB        =  rec_dat$SSB,       # spawning index (usually the year before SSB)
    sdSSB      =  rec_dat$sdSSB,
    sdRec      =  rec_dat$sdRobs,
    covars    = ec,                # (n_cov, n_yrs) matrix
    covars_sd = ec_sd,
    beta_0    = matrix(as.numeric(beta_0IN),dim(ec)[1],dim(ec)[2],byrow=F),
    lambda_0  = matrix(as.numeric(lambda_0IN),dim(ec)[1],dim(ec)[2],byrow=F),
    #startVal  = tibble::deframe(estparams[,c(1,4)]),
    startVal  = NULL,
    phases    = tibble::deframe(estparams[,c(1,2)]))
  
  return(datlist = datlist)
}