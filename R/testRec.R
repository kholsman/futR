## Explore Recruitment formulations:


R_hat<-function(N,gamma=-1,b=1,loga=1,beta=.2,lambda=.2,x){
  sapply(N,function(NIn=N,gammaIN=gamma,bIn=b,logaIn=loga,betaIn=beta,lambdaIn=lambda){
    gN <- NIn*exp(sum(lambdaIn*x))
    hx <- exp(sum(betaIn*x))
    return( exp(loga)*gN* ((1-bIn*gammaIN*gN)^(1/gammaIN))*hx)
  }
)}

bseq<-seq(0,1,.01)
plot(bseq,R_hat(N=100,gamma=-1,b=bseq,a=1,beta=.2,lambda=.2,x=bseq),type="l")

gseq<-seq(-1,1,.01)
Nseq <- seq(0,200,1)
xseq <-Nseq*0
BETA <- 0
LAMBDA <-0
# ricker
plot(Nseq,R_hat(N=Nseq,gamma=.3,b=.01,loga=1,beta=BETA,lambda=LAMBDA,x=xseq),type="l")
# Bev holt
lines(Nseq,R_hat(N=Nseq,gamma=-1,b=.1,loga=2,beta=BETA,lambda=LAMBDA,x=xseq),type="l")
#linear
lines(Nseq,R_hat(N=Nseq,gamma=0,b=.1,loga=-1.3,beta=BETA,lambda=LAMBDA,x=xseq),type="l")
#exponential:
lines(Nseq,R_hat(N=Nseq,gamma=1,b=-0.002,loga=-1.3,beta=BETA,lambda=LAMBDA,x=xseq),type="l",col="red")

lines(Nseq,R_hat(N=Nseq,gamma=0,b=0,loga=-2.3,beta=BETA,lambda=LAMBDA,x=xseq),type="l",col="blue")

x_seq2<-datlist[[6]][[1]]$S_obs*0
plot(datlist[[6]][[1]]$S_obs,R_hat(N=datlist[[6]][[1]]$S_obs,gamma=a,b=-12.3076257,loga=1.2,beta=0,lambda=0,x=x_seq2),type="p")

g<- 0.126577
a<-3.320117
b<-4.517166e-06
S_obs<-rec$S_obs
Acovars<-Bcovars<-0

plot(S_obs,R_hat(N=S_obs,gamma=g,b=b,loga=log(a)/10,beta=0,lambda=0,x=x_seq2),type="p")

R_hat<-S_obs*0
for(i in 1:length(S_obs))
R_hat[i]  =  a*S_obs[i]*exp(Acovars)*( 1 - (b*g*S_obs[i]*exp(Acovars))^(1/g)  )*exp(Bcovars)

plot(S_obs,R_hat)

LAMBDA <-0
nobs  <- length(Nseq)
err   <- rnorm(nobs,0,.2)
AC    <- 1.2
xseq <-exp(err+c(0,err[-nobs])*AC)
xseq <- rep(2,nobs)
xseq2 <- rep(-2,nobs)
plot(xseq,type="l")

plot(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=.002,lambda=LAMBDA,x=xseq),type="l",lty=2,col="red")
lines(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0,lambda=LAMBDA,x=xseq),type="l")
lines(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=.002,lambda=LAMBDA,x=xseq2),type="l",lty=2)

plot(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0,lambda=0.002,x=xseq2),type="l",lty=2)
lines(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0,lambda=0,x=xseq),type="l")
lines(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0,lambda=.002,x=xseq),type="l",lty=2,,col="red")

plot(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0.0,lambda=0,x=xseq2),type="l")
lines(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0.002,lambda=.002,x=xseq),type="l",lty=2)
lines(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0.002,lambda=.002,x=xseq2),type="l",lty=2,,col="red")
lines(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0.002,lambda=-.002,x=xseq2),type="l",lty=2,,col="red")
lines(Nseq,R_hat(N=Nseq,gamma=0.12,b=.1,a=30,beta=0.002,lambda=-.002,x=xseq),type="l",lty=2,,col="red")
x<-seq(0,1,.01)
x1<- log(-log(1-x))
plot(x,x1,type="l",ylim=c(-4,4))
x2<-  log((x/(1-x)))
lines(x,x2,lty=2)
