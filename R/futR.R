#' futR
#'
#' futR() is a generic Rpackage for fitting recruitment models to indices of 
#' spawning biomass (or numbers) and recruitment with or without climate 
#' covariates. The recruitment model is based on Template Model Builder 
#' (`TMB`) and standard stock-recruitment models modified to include the 
#' effects of environmental covariates (Deriso 1980; Schnute 1985). This 
#' includes Ricker (logistic), Beverton Holt, log-linear, and log-linear 
#' with biomass lagged by year 'y-1'. The model can be fit with and with 
#' out random effects on spawning stock biomass (SSB) and recruitment (R) 
#' (i.e., measurement error on SSB and rec) using the methods of  Porch 
#' and Lauretta (2016) and with the optional unbiased estimate of sigma 
#' (sensu Ludwid and Walters 1981, Porch and Lauretta 2016). Environmental 
#' covariates are optional and can be included to influence pre-spawning 
#' and post-spawning success.
#'
#'
#' @docType package
#'
#' @author Kirstin K. Holsman \email{kirstin.holsman@noaa.gov}
#'
#' @name futR
NULL