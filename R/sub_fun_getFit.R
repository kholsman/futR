#' extract the model fits from the runmod() results
#'
#' This function prepares data for plot_rs()
#' For more information contact author Kirstin Holsman (kirstin.holsman at noaa.gov)
#' 
#' @param mm  output from the model
#' @param nm  model name to assign to the model
#' @importFrom dplyr %>%
#' @returns maplist a list with map values
#' 
#' @export 
#' 
#' @examples
#'  #add
#' 
#' 
getFit <- function(mm, nm = "recType = 2"){
  library("dplyr")
  mm  <- mm 
  out <- data.frame(mm$pred%>%
    dplyr::select(def,pred)%>%
    dplyr::group_by(def)%>%
    dplyr::filter(def%in%c("S_hat","S_obs","R_hat","R_obs"))%>%
    dplyr::mutate(id = 1:n())%>%  
    tidyr::spread(def, pred)%>%
    dplyr::mutate(R_obsIN =  mm$input$R_obs, 
         S_obsIN =  mm$input$S_obs, 
         year    =  mm$input$years,model = nm)%>%
    dplyr::ungroup())
  return(out)
}