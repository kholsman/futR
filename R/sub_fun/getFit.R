#'
#'
#'
#'getFit.R
#'
#'


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