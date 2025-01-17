#' plot_par.R
#' 
#' plot the pdf of each parameter estimated in the model
#' 
#' For more information contact author Kirstin Holsman (kirstin.holsman at noaa.gov)
#' 
#' @param df   dataframe of outputs from runmod(), e.g., mm <- runmod(...)
#' * model     model name to be displayed on the legend
#' * estimate  as.vector(mm$sim)
#' * parameter names( mm$mle)[row(mm$sim)]
#' @return list
#' @importFrom dplyr %>%
#' @importFrom ggplot2 aes
#'
#' @examples
#' # read in the data
#' 
#' #datlist <- datlist2 <-  
#' #readMake_futR_data("data/in/futR_Inputs.xlsx" )
#' 
#' # set the sigMethod to 1
#' #now compare plots
#' # plot_par_pdf(rbind(dfR4_s1,dfR4_s3))
#' @export
plot_par_pdf <- function(df = df1_t0){
  mu   <- df%>%
    dplyr::group_by(model,parameter)%>%
    dplyr::summarise(grp.mean=mean(estimate))
  peak <- df%>%
    dplyr::group_by(model,parameter)%>%
    dplyr::count(parameter,round(estimate,1))%>%
    dplyr::slice(which.max(n))
  names(peak)<- c("model","parameter","freq","n")
  
  p <-ggplot2::ggplot(data=df) +
    ggplot2::geom_density( aes(x=estimate, color=model))+
    ggplot2::facet_wrap(~parameter,scales="free")+
    ggplot2::geom_vline(data=mu,aes(xintercept=grp.mean, color = model), 
                        linetype="solid", size=1)+
    ggplot2::geom_vline(data=peak,aes(xintercept=freq, color = model), 
                        linetype="dashed", size=.6)+
    ggplot2::theme_minimal()
    # theme_kir_EBM()
  p
}