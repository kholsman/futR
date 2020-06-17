

# rm(list=ls())
# graphics.off()
setwd("/Users/kholsman/GitHub/futR/")
load("rs_dat.Rdata")


######################################
#phase example

# PHase 1 - no covars"
	parameters     <-  list(log_aa_c=log(1),log_bb_c=2,rs_parm=0,logsigma=log(.9))
	require(TMB) 
	# precompile()
	compile('futR.cpp') 

	model<-MakeADFun(data=rs_dat,parameters=parameters,DLL='futR',map=list(
      log_aa_c=factor(1),
      log_bb_c=factor(2),
      rs_parm=factor(NA),
      logsigma=factor(2))) 
	fit<-nlminb(model$env$last.par.best,model$fn,model$gr) 
	rep2<-sdreport(model,getJointPrecision=TRUE)
	out<-summary(rep2)
