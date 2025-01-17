#'predict_futR.R
#'
#'
#' Predict recruitment from new data and or new par
#' 
#' @param modelIN   model output from runRecMod()
#' @param newdata   new predictive data to simulate (can be new S_obs or cov)                    
#' @param src_fldr  src folder path
#' @param futR_fldr NULL
#' @param  newpar   new parameters to run through the mode description
#' 
#' @returns  pred_model$report()
#' 
#' @export
#' 
predict_futR <- function(
                  modelIN = mm,  # model output from runRecMod()
                  newdata,       # new predictive data to simulate (can be new S_obs or cov)                    
                  src_fldr  = "src/TMB",
                  futR_fldr =  NULL,
                  newpar){        # new parameters to run through the mode
  tmpdir<- getwd()
  tryCatch({
    # get objects from the fitted model
    obj <- modelIN$model
    ee  <- environment(obj$fn)
    lp  <- ee$last.par.best  ## used in $report() call below
    dd  <- ee$data           ## data object
    # obj$report(lp)
    if(!is.null(futR_fldr)) setwd(futR_fldr)
    if(!is.null(src_fldr)) setwd(src_fldr)
     dyn.load( dynlib(ee$DLL))
     
      pred_model  <- suppressWarnings(MakeADFun(
        data                 =  newdata,
        parameters           =  newpar,
        DLL                  =  ee$DLL,
        checkParameterOrder  =  TRUE,
        hessian              =  FALSE,
        map                  =  ee$map,
        silent               =  TRUE))
      
      pred_model$report()
     
    },
    error = function(cond) {
      message("Problem with predict, error message:")
      message(conditionMessage(cond))
      # Choose a return value in case of error
      NA
    },
    warning = function(cond) {
      message("Problem with predict, warning message:")
      message(conditionMessage(cond))
      # Choose a return value in case of warning
      NULL
    },
    finally = {
      setwd(tmpdir)
    }
  )

}


