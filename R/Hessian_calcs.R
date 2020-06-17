library(r4ss);library("corpcor") #first get these libraries
install.packages("corpcor",repos="https://cloud.r-project.org",dependencies=TRUE)

GET_HESS_AIC<-function(flnm){
	#Based on Jim Thorson code via ingrid
	require("corpcor")

	#source code - use this function.
	get.admb.hes <- function(flnm){
	 wd.old <- getwd(); on.exit(setwd(wd.old))
	 setwd(model.path)
	# filename <- file("admodel.hes", "rb")
	 filename <- file (flnm,"rb")
	 on.exit(close(filename), add=TRUE)
	 num.pars <- readBin(filename, "integer", 1)
	 hes.vec <- readBin(filename, "numeric", num.pars^2)
	 hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
	 hybrid_bounded_flag <- readBin(filename, "integer", 1)
	 scale <- readBin(filename, "numeric", num.pars)
	 result <- list(num.pars=num.pars, hes=hes,
	                hybrid_bounded_flag=hybrid_bounded_flag, scale=scale)
	 return(result)
	}

	#this pulls out the hessian and then transforms it into parameter space
	HESS = get.admb.hes()
	# Calculate Hessian
	cov <- pseudoinverse(HESS$hes)
	#sqrt(diag(HESS$hes))
	scale <- HESS$scale
	cov.bounded <- cov*(scale %o% scale)
	Hess = pseudoinverse(cov.bounded)


	NLL=4759.74 #this is the objective function value
	num.pars=120
	LnDet = determinant(Hess, logarithm=TRUE)$modulus[[1]]
	Ln_Integral = log(2*pi) + (-1/2)*sum(LnDet) + -1*NLL#this is the MARGINAL likelihood
	AIC = -2*Ln_Integral + 2*num.pars
	return(AIC)
}
#15.1c
NLL=3824.53 #this is the objective function value
num.pars=167
AIC=8078.448

#18.0
NLL=3831.72 #this is the objective function value
num.pars=159
AIC=8406.3

#18.1
NLL=3831.48 #this is the objective function value
num.pars=121
AIC=8384.0

#18.2
NLL=3876.08 #this is the objective function value
num.pars=119
AIC=8459.202

#18.3
NLL=3747.05 #this is the objective function value
num.pars=167
AIC=7817.843

#18.4
NLL=3801.56 #this is the objective function value
num.pars=167
AIC=-inf

#18.5
NLL=6473.34 #this is the objective function value
num.pars=167
AIC=-inf

#18.6
NLL=3730.32 #this is the objective function value
num.pars=159
AIC=8208.4


#here num.pars is the number of parameters from the .par file
#I think it is usually okay, except here I actually constrained the number of parameters
#because of the way I estimated selectivity so that is why Grant said I actually had 
#fewer than that.