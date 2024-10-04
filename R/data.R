#' @title datlist
#'
#' @description A data input file for the futR model
#'
#' @format A list with the following objects
#' \describe{
#'   \item{parameters}{a list with the starting values for each model parameter}
#'   \item{parameters$loga}{starting value for the log of the RS intercept}
#'   \item{parameters$logb}{starting value for the log of the RS slope}
#'   \item{parameters$lambda}{starting value for the env co-variate parms}
#'   \item{parameters$epsi_s}{starting value for the espsi parm}
#'   \item{parameters$logsigma}{starting value for the log of the error term}
#'   \item{parameters$skipFit}{a 0=NO 1=YES switch for if the fit should be skipped (simulate)}
#'   \item{rs_dat}{a list with each of the data objects needed to run the model}
#'   \item{rs_dat$tau}{observation error effect 0 = no observation error; >0 tau value used in sigMethod 2-5 (random effects on SSB), }
#'   \item{rs_dat$sigMethod}{index of 1-5, 
#'        1. No observation error (𝜏 = 0)
#'        2. estimate sigma, random effects on SSB if 𝜏 > 0, 𝜏 input
#'        3. as in 1 but with unbiased sigma estimate, 𝜏 input
#'        4. as in 1 but with defined measurement error for rec (indep of random effects on Spawners/SSB)
#'        5. as in 1 but with defined measurement error for rec and Spawners/SSB)}
#'   \item{rs_dat$ncov}{number of environmental covariates}
#'   \item{rs_dat$cov_type}{depreciated covariate type if 1 then run ration cov_type.allocate( 1,ncov,"cov_type" );}
#'   \item{rs_dat$nyrs}{number of recruitment years}
#'   \item{rs_dat$sdrs_cov}{optional standard deviation of the environmental covariates}
#'   \item{rs_dat$rs_cov}{matrix of covariate values with rows for each covariate and columns for each year}
#'   \item{rs_dat$sdR}{recruitment standard deviation (optional)}
#'   \item{rs_dat$sdS}{spawning biomass standard deviation (optional)}
#'   \item{rs_dat$R_obs}{recruitment}
#'   \item{rs_dat$S_obs}{spawning biomass}
#'   \item{rs_dat$years}{years}
#'   \item{rs_dat$rectype}{index of 1-4
#'          1. linear model (not Spawning biomass effect)
#'          2. Linear with biomass ( y-1 ) as a predictor
#'          3. Beverton Holt
#'          4. Ricker model}
#'   \item{rs_dat$estMode}{0/1 for estimate (0) or skip nll (1)}
#'   \item{rs_dat$beta_0}{0,1 matrix of inclusion of prespwaning covariate effects: beta_0.allocate( 1,ncov,1,nyrs,"beta_0" );}
#'   \item{rs_dat$lambda_0}{0,1 matrix of inclusion of postspwaning covariate effects: beta_0.allocate( 1,ncov,1,nyrs,"beta_0" );}
#'   \item{maplist}{a TMB map list of the parameters to be estimated}
#'   \item{estparams}{a list of T/F estimate each par }
#'   \item{phases}{phases for estimating each par (currently in devel)}
#' }
#' @source K. Holsman et al. 2024
"datlist"

#' @title rec
#'
#' @description 
#'
#' @format A dataframe with the following columns
#' \describe{
#'   \item{rec_year}{ year of recruitment}
#'   \item{SSB}{Spawning biomass in the previous year}
#'   \item{Robs}{Recruitment in the year}
#'   \item{sdSSB}{Error in the Spawning biomass estimate}
#'   \item{sdRobs}{Error in the recruitment estimate}
#'   \item{fit}{0/1 no/yes fit in the model?}
#' }
#' @source K. Holsman et al. 2024
"rec"
