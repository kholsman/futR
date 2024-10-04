#' extract the model fits from the runmod() results
#'
#' This function prepares data for plot_rs()
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' 
#' @param mm  output from the model
#' @param nm  model name to assign to the model
#' 
#' @returns maplist a list with map values

#' @export 
#' 
#' @examples
#'  datlist <- readMake_futR_data("data/in/futR_Inputs.xlsx" )
#'  mm      <- runmod(dlistIN   = datlist, version   = 'futR',recompile = FALSE,simulate  = TRUE,sim_nitr  = 1000)  
#'  r_fit  <- getFit(mm, nm = paste0("recType = ",datlist$rs_dat$rectype))
#'  plot_rs(r_fit)
#' 
#' 
getFit <- function(mm, nm = "recType = 2"){
  out <- data.frame(mm$pred %>%
  dplyr::select(def,pred)%>%
    dplyr::group_by(def) %>%
    dplyr::filter(def%in%c("S_hat","S_obs","R_hat","R_obs"))%>%
    dplyr::mutate(id = 1:n()) %>%  
  tidyr::spread(def, pred)%>%
    dplyr::mutate(R_obsIN =  mm$input$R_obs, 
         S_obsIN =  mm$input$S_obs, 
         year    =  mm$input$years,model = nm)%>%
    dplyr::ungroup())
  return(out)
}