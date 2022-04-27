#' extract the model fits from the runmod() results
#'
#' This function prepares data for plot_rs()
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @weblink 
#' @param parameters  list of starting values and parameters for the TMB futR() model
#' @param estpar      vector of T/F of which parameters will be estimated
#' @export            maplist a list with map values
#' @email             For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @examples
#'  datlist <- readMake_futR_data("data/in/futR_Inputs.xlsx" )
#'  mm      <- runmod(dlistIN   = datlist, version   = 'futR',recompile = FALSE,simulate  = TRUE,sim_nitr  = 1000)  
#'  r_fit  <- getFit(mm, nm = paste0("recType = ",datlist$rs_dat$rectype))
#'  plot_rs(r_fit)
#' @export
getFit <- function(mm, nm = "recType = 2"){
  out <- mm$pred%>%
  select(def,pred)%>%
  group_by(def) %>%
  filter(def%in%c("S_hat","S_obs","R_hat","R_obs"))%>%
  mutate(id = 1:n(), model = nm) %>%  
  tidyr::spread(def, pred)%>%
  mutate(R_obsIN =  mm$input$R_obs, 
         S_obsIN =  mm$input$S_obs, 
         year =  mm$input$years) 
  return(out)
}