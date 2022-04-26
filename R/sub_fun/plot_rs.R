#'
#'
#'plot_rs
#'
#'

plot_rs <- function(rec_df = rec_fit, cbp = "Dark2", obscol = "gray"){
  Rmax <- 1.05* max( max(rec_df$R_obs), max(rec_df$R_hat))
  Smax <- 1.05* max( max(rec_df$S_obs), max(rec_df$S_hat))
  
  prs <- ggplot(data =  rec_df) +
    geom_point (aes(x=S_obs,y = R_obs),color=obscol,size =2)+ 
    xlab("Spawner Index")+
    ylab("Recruitment")+    
    coord_cartesian(ylim = c(0, Rmax),xlim = c(0, Smax))+
    geom_point (aes(x=S_obs,y = R_hat, color = model))+ 
    theme_minimal()+
    # theme_kir_EBM()+
    scale_color_brewer(palette = cbp)

  
  prt <- ggplot(data =  rec_df) +
    geom_point (aes(x=year,y = R_obs),color=obscol,size =2)+ 
    coord_cartesian(ylim = c(0, Rmax))+
    xlab("Year")+
    ylab("Recruitment")+    
    geom_point (aes(x=year,y = R_hat, color = model))+ 
    geom_line (aes(x=year,y = R_hat, color = model))+ 
    theme_minimal()+
    # theme_kir_EBM()+
    scale_color_brewer(palette = cbp)
  #display.brewer.all()
  #col(1)
  library("cowplot")
  plot_grid(prs, prt, #bp + rremove("x.text"), 
            labels = c("A", "B"),
            ncol = 1, nrow = 2)
  
  
}