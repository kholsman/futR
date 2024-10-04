#' plot_par.R
#' 
#' plot the pdf of each parameter estimated in the model
#' 
#' @email For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @weblink 
#' @param df   dataframe of outputs from runmod(), e.g., mm <- runmod(...)
#' * model     model name to be displayed on the legend
#' * estimate  as.vector(mm$sim)
#' * parameter names( mm$mle)[row(mm$sim)]
#' @return list
#'
#' @examples
#' # read in the data
#' datlist <- datlist2 <-  readMake_futR_data("data/in/futR_Inputs.xlsx" )
#' 
#' # set the sigMethod to 1
#' datlist2$rs_dat$sigMethod <-1
#' mm      <- runmod(dlistIN   = datlist, version   = 'futR',recompile = FALSE,simulate  = TRUE,sim_nitr  = 1000)
#' dfR4_s1 <-  data.frame(model = "sigMethod 1",estimate  = as.vector(mm$sim),parameter = names( mm$mle)[row(mm$sim)])
#' plot_par_pdf(dfR4_s1)
#' 
#' #set the sigMethod to 4 (unbiased sigma)
#' datlist2$rs_dat$sigMethod <- 4
#' mm      <- runmod(dlistIN   = datlist2, version   = 'futR',recompile = FALSE,simulate  = TRUE,sim_nitr  = 1000)
#' dfR4_s3 <-  data.frame(model = "sigMethod 4",estimate  = as.vector(mm$sim),parameter = names( mm$mle)[row(mm$sim)])
#' #now compare plots
#' plot_par_pdf(rbind(dfR4_s1,dfR4_s3))
#' @export
plot_par_pdf <- function(df = df1_t0){
  mu   <- df%>%group_by(model,parameter)%>%summarise(grp.mean=mean(estimate))
  peak <- df%>%group_by(model,parameter)%>%
    count(parameter,round(estimate,1))%>%
    slice(which.max(n))
  names(peak)<- c("model","parameter","freq","n")
  
  p <-ggplot(data=df) +
    geom_density( aes(x=estimate, color=model))+
    facet_wrap(~parameter,scales="free")+
    geom_vline(data=mu,aes(xintercept=grp.mean, color = model), linetype="solid", size=1)+
    geom_vline(data=peak,aes(xintercept=freq, color = model), linetype="dashed", size=.6)+
    theme_minimal()
    # theme_kir_EBM()
  p
}