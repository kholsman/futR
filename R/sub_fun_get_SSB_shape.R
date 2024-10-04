#' get_SSB_shape.R
#'
#' get the SSB / R shape
#'
#' @param x A description of the parameter 'x'. The
#'   description can span multiple lines.
#' @param y A description of the parameter 'y'
#' @param AIC_summryIN ,
#' @param simname = "hindcast 2020",
#' @param SD_range = 7,
#' @param modIN = AIC_summry[[1]]$topRicker_R2,
#' @param recfitIN = AIC_summry[[1]]$topRicker_R2_rec_fit,
#' @param simulateIN= F,
#' @param sim_nitrIN = 200,
#' @param NewCovsIN =hind%>%filter(year==2020)%>%data.frame(), 
#' @param steps   = 100,
#' @param futR_path = "../../futR")
#' 
#' @returns The sum of `x` and `y`.
#' @export
#'
#' @examples

get_SSB_shape <- function(
    AIC_summryIN,
    simname = "hindcast 2020",
    SD_range = 7,
    modIN = AIC_summry[[1]]$topRicker_R2,
    recfitIN = AIC_summry[[1]]$topRicker_R2_rec_fit,
    simulateIN= F,
    sim_nitrIN = 200,
    NewCovsIN =hind%>%filter(year==2020)%>%data.frame(), 
    steps   = 100,
    futR_path = "../../futR"){
  
  tmdir <- getwd()
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
        SD_range    =  SD_range,
        steps       = steps,
        hind_yearsIN = 1980:2023,
        simulateIN = F,
        modIN      = modIN,
        recIN      = recfitIN,
        phasesIN   = phases)
      maxR <- max(out$R_hat)
      minS <- out$muSSB[1]-1.95*out$sdSSB[1]
      maxS <- out$muSSB[1]+1.95*out$sdSSB[1]
    p <- ggplot()+
      geom_line(data=out, aes(x=SSB,y=R_hat,color=name)) +
      geom_line(data=out[ which(out$SSB>=minS&out$SS<=maxS),], aes(x=SSB,y=R_hat,color=name),size=2) +
      theme_minimal()+theme(legend.title = element_blank())
      # geom_polygon(aes(
      #   y = c(0,maxR,maxR,0),
      #   x = c(rep(minS,2),rep(maxS,2)), color="avg. 1980:2020"),alpha=.4, fill = "white"  ) + 
      # geom_line(data=out, aes(x=SSB,y=R_hat,color=name)) +
      
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

