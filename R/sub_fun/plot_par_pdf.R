#'
#'
#'plot_par_pdf.R
#'



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