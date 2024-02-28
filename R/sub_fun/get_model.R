#'
#'
#' get_model.R
#' 
#' get the model from the meta evaluation 
#' 


get_model <- function(
    mod        ,
    recIN      ,
    simulateIN = FALSE,
    sim_nitrIN = 1000,
    futR_fldr =  NULL,
    phasesIN   = phases){
    
    #--------------------------
    
    betaIN <- lambdaIN <- FALSE
    
    if(mod$type == "postspawn effects")  {
      mod_mat  <- postspawn_mod_mat
      lambdaIN <- TRUE
    }
    if(mod$type == "prespawn effects")  { 
      mod_mat  <- prespawn_mod_mat
      betaIN   <- TRUE
      
    }
    if(mod$type == "pre- and postspawn effects") { 
      mod_mat  <- preNpostspawn_mod_mat
      lambdaIN <-  betaIN   <- TRUE
    }
    
    dlist  <-  makeDat(
      rectype    =  mod$rectype, 
      tauIN      =  0, # no obs error
      sigMethod  =  1,  # no obs error
      estparams = c(
        log_a        = TRUE, 
        log_b        = TRUE, 
        beta         = betaIN,
        lambda       = lambdaIN,
        epsi_s       = FALSE,
        logsigma     = TRUE,
        skipFit      = FALSE),
      estMode    =  1,
      rec_years  =  recIN$year,  # recruitment year
      fityrs     =  recIN$year ,         # fit the model to a subset of years
      Rec        =  recIN$R_obs,      # recruitment index
      SSB        =  recIN$S_hat,       # spawning index (usually the year before SSB)
      sdSSB      =  recIN$S_hat*0,
      sdRec      =  recIN$R_obs*0,
      covars    = t(mod_mat[[mod$modn]]),                # (n_cov, n_yrs) matrix
      covars_sd = t(mod_mat[[mod$modn]])*0,
      beta_0     = rep(1,dim(mod_mat[[mod$modn]])[2])*0,
      lambda_0   = rep(1,dim(mod_mat[[mod$modn]])[2]),
      REcompile = FALSE,
      startVal  = NULL,
      phases    = phasesIN)
    
    tmpdir2 <- getwd()
    tryCatch(
      {
        if(!is.null(futR_fldr) ) setwd(futR_fldr)
        object   <- suppressMessages(suppressWarnings(
          runRecMod(dlistIN   = dlist,
                    version   = 'futR',                           
                    src_fldr  = "src/TMB",
                    recompile = FALSE,
                    simulate  = simulateIN,
                    sim_nitr  = sim_nitrIN))
        )
        object
      },
      error = function(cond) {
        message("Error in runRecMod, error message:")
        message(conditionMessage(cond))
        # Choose a return value in case of error
        NA
      },
      warning = function(cond) {
        message("Warning in runRecMod, warning message:")
        message(conditionMessage(cond))
        # Choose a return value in case of warning
        NULL
      },
      finally = {
        setwd(tmpdir2)
      }
    )
    return(list(object= object, dlist=dlist))
}
    
    
    
    
    
    