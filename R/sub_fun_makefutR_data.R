#' makefutR_data.R
#' 
#' reads an xlxs input file and creates the datlist for input into runmod()
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#'  
#' @param fn path to data input sheet; default is 'data/in/futR_Inputs.xlsx'
#'
#' @return datlist  list of input data and parameters used for the runmod() dependency
#' * datlist$parameters  list of parameters for the TMB futR model:
#'    [] datlist$parameters$log_a    scalar starting value for log_a (intercept)
#'    [] datlist$parameters$log_b    scalar starting value for log_b (slope) 
#'    [] datlist$parameters$beta     vector of starting values for lambda (env. effects on pre-spawning)   
#'    [] datlist$parameters$lambda   vector of starting values for beta (env. effects on post-spawning)   
#'    [] datlist$parameters$epsi_s   starting value for epsi_s (error around S estimates) 
#'    [] datlist$parameters$logsigma starting value for logsigma (process error)
#'    [] datlist$parameters$skipFit  starting value for skipFit (switch used for projecting)
#' * datlist$rs_dat list of data to read into the model
#'    --- Observation Error inputs
#'    [] tau        scalar for observation and process error
#'    [] sigMethod  Integer (switch): which method for estimating sigma
#'    --- Covariate inputs
#'    [] ncov       Integer number of environmental covariates
#'    [] cov_type   (defunct) vector of covariate type if 1 then run ration
#'    [] nyrs       Integer number of years in the data
#'    [] lambda_0   matrix (ncov, nyrs), binary (0,1) matrix of covar. inclusion in post-spawning effects
#'    [] beta_0     matrix (ncov, nyrs), binary (0,1) matrix of covar. inclusion in pre-spawning effects
#'    [] sdrs_cov   matrix (ncov, nyrs) of the stdev of environmental covariates
#'    [] rs_cov     matrix (ncov, nyrs) of the of environmental covariates
#'    --- RS data inputs
#'    [] sdR        vector (nyrs) of stdev of recruitment in a given year
#'    [] sdS        vector (nyrs) of stdev of spawnerindex in a given year
#'    [] R_obs      vector (nyrs) of recruitment in a given year
#'    [] S_obs      vector (nyrs) of spawner index in a given year
#'    [] years      vector (nyrs) recruiment year
#'    --- Switches
#'    [] rectype    Integer (switch) which type of recruitment model: 1 = loglinear, 2 = loglinear ~f(S), 3 = Beverton-holt, 4 = Ricker
#'    [] estMode    Integer (switch) should the model run in estimation mode (y=1,no=0)?
#' * datlist$maplist   list of the map for the TMB model
#' * datlist$estparams vector of T/F of which parameters will be estimated
#' * datlist$phases    (defunct) phases for estimating parameters
#' * datlist$inputs    archive of intputs for this function
#' @examples
#' datlist <- readMake_futR_data("data/in/futR_Inputs.xlsx" )
#' mm      <- runmod(datlist   = datlist, version   = 'futR',recompile = TRUE,simulate  = TRUE,sim_nitr  = 1000)  
#' @export
#' 

makefutR_data <- function(fn = "data/in/futR_Inputs.xlsx",export_all=F){
  
  tmp            <- data.frame(readxl::read_xlsx(fn, sheet = "rec_data" ))
  ec             <- data.frame(readxl::read_xlsx(fn, sheet = "covar_val" )) 
  ec_sd          <- data.frame(readxl::read_xlsx(fn, sheet = "covar_val_sd")) 
  beta_0IN       <- data.frame(readxl::read_xlsx(fn, sheet = "pre_post_spawning_effects", n_max =2)[1,-1])
  lambda_0IN     <- data.frame(readxl::read_xlsx(fn, sheet = "pre_post_spawning_effects", n_max =2 )[2,-1])
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
  
  datIN <-list()
  datIN$ec      <- ec
  datIN$ec_sd   <- ec_sd
  
  datIN$beta_0IN   <- beta_0IN
  datIN$lambda_0IN <- lambda_0IN
  datIN$switches   <- switches
  datIN$rec_dat    <- rec_dat
  datIN$fityears   <- fityears
  
  
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
   # beta_0    = matrix(as.numeric(beta_0IN),dim(ec)[1],dim(ec)[2],byrow=F),
  #  lambda_0  = matrix(as.numeric(lambda_0IN),dim(ec)[1],dim(ec)[2],byrow=F),
    beta_0     = beta_0IN,
    lambda_0   = lambda_0IN,
    #startVal  = tibble::deframe(estparams[,c(1,4)]),
    startVal  = NULL,
    phases    = tibble::deframe(estparams[,c(1,2)]))
  if(!export_all)
    return(datlist = datlist)
  if(export_all)
    return(list(datlist = datlist,datIN = datIN))
}