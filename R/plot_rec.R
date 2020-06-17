### Ploting code for futR

library(wesanderson)
plt<-c("Zissou1","Darjeeling1","Darjeeling2","FantasticFox1")
dawn<-colorRampPalette(c(colors()[c(477,491,72,474,47,653)],"orange","red"))
Ornjazz <- RColorBrewer::brewer.pal(5, "Oranges")
orng<-colorRampPalette(Ornjazz[1:5])
blues <- RColorBrewer::brewer.pal(5, "Blues")
BG <- RColorBrewer::brewer.pal(9, "GnBu")  #5
wes<-colorRampPalette(c(wes_palette(n=5, name=plt[1])[1:5]))
col2<-colorRampPalette(c(wes(7)[c(3,1)],orng(3)[3]))
col3<-colorRampPalette(c(wes(7)[4:7]))
coll_use         <-  c(colors()[320],col2(6)[c(2,3,4)],col3(6)[c(3,4,6)])

# load Kir's theme:
source("~/Documents/D_AFSC_Files/AFSC_code/R_scripts/Misc_scripts/THEMES_GGPLOT.r")



# Plot comparisons of recruitment from ADMB and TMB
plotCompareMod<-function(H=7,W=6,mode=1,outIN=out,recIN=rec){
  dev.new(height=H,width=W)
  captionIN  <- "Comparison of (-1,1) standard dev. of covariates on Ricker RS for each species"
  titleIN    <- paste(names(mod[[1]])[mode]," recruitment estimates ")
  
  outIN   <-  outIN[outIN$mode==names(mod[[1]])[mode],]
  recIN   <-  recIN[recIN$mode==names(mod[[1]])[mode],]
  
  out_tmp <-  as_tibble(outIN)
  out_tmp <-  out_tmp%>%filter(mode==levels(out_tmp$mode)[mode])
  rec_tmp <-  as_tibble(recIN)
  rec_tmp <-  data.frame(rec_tmp%>%filter(mode==levels(rec_tmp$mode)[mode]))
  
  types<-levels(out_tmp$type)
  
  upp <- out_tmp%>%filter(type==types[2],X1=='95%')
  dwn <- out_tmp%>%filter(type==types[2],X1=='5%')
  mid <- out_tmp%>%filter(type==types[2],X1=='50%')
  mid2 <- out_tmp%>%filter(type==types[1],X1=='50%')
  upp2 <- out_tmp%>%filter(type==types[1],X1=='95%')
  dwn2 <- out_tmp%>%filter(type==types[1],X1=='5%')
  out_tmp <- data.frame(out_tmp)
  p <- ggplot(data=rec_tmp,aes(x=SSB,y=Rec))
  p <- p + facet_wrap(sp~mod,scales="free",nrow=3)+ expand_limits(y = 0) 
  p <- p + geom_point(aes(x=SSB,y=Rec),col=coll_use[1])
  
  # p <- p + geom_line(data=mid,aes(x=value.1,y=value),col=coll_use[5],inherit.aes = FALSE)
  # p <- p + geom_ribbon(data=upp,aes(x=dwn$'value.1',ymin=dwn$value,ymax=upp$value),fill=coll_use[5],alpha=.2,inherit.aes = FALSE)
  # p <- p + geom_line(data=mid2,aes(x=value.1,y=value),col=coll_use[3],inherit.aes = FALSE)
  # p <- p + geom_ribbon(data=upp2,aes(x=dwn2$'value.1',ymin=dwn2$value,ymax=upp2$value),fill=coll_use[3],alpha=.2,inherit.aes = FALSE)
  # 
  p <- p + geom_line(data=mid,aes(x=SSB,y=value),col=coll_use[5],inherit.aes = FALSE)
  p <- p + geom_ribbon(data=upp,aes(x=SSB,ymin=dwn$value,ymax=upp$value),fill=coll_use[5],alpha=.2,inherit.aes = FALSE)
  p <- p + geom_line(data=mid2,aes(x=SSB,y=value),col=coll_use[3],inherit.aes = FALSE)
  p <- p + geom_ribbon(data=upp2,aes(x=SSB,ymin=dwn2$value,ymax=upp2$value),fill=coll_use[3],alpha=.2,inherit.aes = FALSE)
  
  p <-  p + theme_light() +
    labs(x=NULL, y=NULL,title=titleIN,
         subtitle=NULL,
         caption=captionIN) +
    theme(plot.subtitle=element_text(margin=margin(b=20))) +
    theme(legend.title=element_blank()) +
    theme(legend.position="right") +
    theme(legend.key.width = unit(.5, "cm")) +
    theme(legend.text=element_text(size=5)) +
    theme(legend.key.size=unit(.01, "cm")) +
    theme(plot.margin=margin(t = 10, r = 10, b = 10, l =10)) 
  
  p<- p+ theme_kir_EBM(sub_title_size=12,
                       sub_title_just="l",
                       axis_title_just = "cm") +
    theme(legend.title=element_blank(),
          legend.background = element_rect(colour = NA),
          legend.key = element_rect(colour = "white", fill = NA)) 
  p
}
