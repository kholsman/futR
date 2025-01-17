#' get the shape of the recruitment curve
#'
#' This function will get the SSB R shape
#' For more information contact author Kirstin Holsman (kirstin.holsman at noaa.gov)
#' 
#' 
#' @export
#' 
#' @param AIC_summryIN AIC summary to input
#' @param simname simulation name
#' @param SD_range range of the SD
#' @param modIN model input, e.g., AIC_summry[[1]]$topRicker_R2
#' @param recfitIN recruitment model fit, e.g., AIC_summry[[1]]$topRicker_R2_rec_fit
#' @param simulateIN simulate T or F
#' @param sim_nitrIN number of simulations
#' @param NewCovsIN covariates to input, e.g, hind%>%dplyr::filter(year==2020)%>%data.frame()
#' @param steps number of steps
#' @param phases Input of phases for each parameter
#' @param futR_path = path for the futR model
#' 
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_line
#' 
#' @returns The sum of `x` and `y`
#'
#' @examples
#' #get_SSB_shape()

get_SSB_shape <- function(
    AIC_summryIN,
    simname = "hindcast 2020",
    SD_range = 7,
    modIN,
    recfitIN,
    simulateIN= F,
    sim_nitrIN = 200,
    NewCovsIN, 
    steps   = 100,
    phases, 
    futR_path = "../../futR"){
  
  stepsIN <- steps
  tmdir   <- getwd()
  setwd(futR_path)
  
  tryCatch({
    AIC_summry <- AIC_summryIN
    covs       <- strsplit(as.character(modIN$name),":")[[1]][2]
    covs       <- strsplit(covs,"_PLUS_")[[1]]
    covs2      <- covs[-grep("x2",covs)]
    
    T2IN <- F
    for(c in 1:length(covs2)){
      nn   <- grep(covs2[c],covs)
      if(length(nn)>1) T2IN <-T
    }
  
      out <- profile_ssb(cov_nm_list = covs2,
        nameIN      = simname,
        newdataIN   = NewCovsIN,
        T2          = T2IN,
        SD_range    = SD_range,
        steps       = stepsIN,
        hind_yearsIN = 1980:2023,
        simulateIN = F,
        modIN      = modIN,
        recIN      = recfitIN,
        phasesIN   = phases)
      
      maxR <- max(out$R_hat)
      minS <- out$muSSB[1]-1.95*out$sdSSB[1]
      maxS <- out$muSSB[1]+1.95*out$sdSSB[1]
      
    p <- ggplot2::ggplot(data=out)+
      ggplot2::geom_line(aes(x=SSB,y=R_hat,color=name)) +
      ggplot2::geom_line(data=out[ which(out$SSB>=minS&out$SS<=maxS),], 
                         aes(x=SSB,y=R_hat,color=name),size=2) +
      ggplot2::theme_minimal()+theme(legend.title = element_blank())
      
    return(list(p=p, out = out))
  },
  error = function(cond) {
    message("error message:")
    message(conditionMessage(cond))
    # Choose a return value in case of error
    NA
  },
  warning = function(cond) {
    
    message("warning message:")
    message(conditionMessage(cond))
    # Choose a return value in case of warning
    NULL
  },
  finally = {
    setwd(tmpdir)
  })
}

