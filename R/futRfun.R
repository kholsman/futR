
dirlist<-dir(file.path(main,"R"))
dirlist<-dirlist[dirlist!="futRfun.R"]

dirlist<-c("estRec.R","fitRec.R","makeDat.R","makeMap.R","plot_rec.R","predict_futR.R","runmod.R","vec2list.R") 

print(paste(dirlist,collapse=", "))
for(fl in dirlist) {
  print(paste("loading",fl))
  source(file.path(main,"R",fl))
}