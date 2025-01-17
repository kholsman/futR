#' AICselection
#'
#' AICselection runs AIC information criterion analysis on resulting set of recruitment models
#' 
#' @param LL vector of logliklihoods of each model
#' @param npar number of parameteris in a given model
#' @param n number of data observations (nobs)
#' @param R2 default is NULL, vector of R2 values for each model
#' @param mnames1 names for each model
#' @param include_marginal T/F include mariginal AIC
#' @param type2   if 1 then AIC is calculated, if 2 then AICc (for small sample sizes is calculated)
#' @param covnm    default is NULL, covariate num
#' @param rsType  default is NULL, recruitment type
#' @param hypoth  default is NULL, hypothesis
#' @param LnDet   the natural logarithm of the determinant of a matrix  X;  defualt is NULL
#' @param CE Climate enhanced model? default = NULL
#'    list of input data from makeDat() function: "parameters" "rs_dat"     "maplist"    "estparams"  "phases"  
#' @return   a dataframe of summary statistics for AICc of submodels         
#' 
#' @export
#'
#' @examples
#' #add
#' 
AICselection <- function (LL, 
                          npar, 
                          n, 
                          R2 = NULL,
                          mnames1 = NULL, 
                          include_marginal = FALSE,
                          type2   = 2, 
                          covnm   = NULL,
                          rsType  = NULL,
                          hypoth  = NULL,
                          LnDet   = NULL,
                          CE      = NULL) 
{
  # LL<--1*LL
  # 
  aicfun<-function (npar, LL, n, type = 2) 
  {
    if (type == 1) 
      return((2 * npar - 2 * (LL)))
    if (type == 2) 
      return((2 * npar - 2 * (LL)) + (2 * npar * (npar + 1))/(n - npar - 1))
  }
  
 
  aicfun_marg<-function (npar, LL, n, type = 2,LnDet) 
  {
    #LnDet = determinant(HESS, logarithm=TRUE)$modulus[[1]]
    Ln_Integral <- log(2*pi) + (-1/2)*(LnDet) + LL #this is the MARGINAL likelihood
    
    if (type == 1) 
      return((2 * npar - 2 * (Ln_Integral)))
    if (type == 2) 
      return((2 * npar - 2 * (Ln_Integral)) + (2 * npar * (npar + 1))/(n - npar - 1))
  }
  
  nn         <- length(LL)
  tmp1       <- data.frame(LL = LL, npar = npar, lab = 1:nn, R2 = R2)
  tmp1$name  <- mnames1
  tmp1$aicc  <- tmp1$aicc_marg <- rep(0, nn)
  for (i in 1:nn) {
    tmp1$aicc[i]        <- aicfun(npar[i], -LL[i], n[i], type = type2)
    if(!is.null(LnDet))
      tmp1$aicc_marg[i] <- aicfun_marg(npar=npar[i], LL=-LL[i],n= n[i], type = type2,LnDet=LnDet[i])
   # tmp1$aicc_marg[i] <- GET_HESS_AIC(HESS[[i]], npar = npar[i],NLL = -1 * LL[i])[[1]]
  }
  tmp1$deltaAIC  <- (tmp1$aicc - min(tmp1$aicc, na.rm = T))
  tmp1$AICweight <- exp(-0.5 * tmp1$deltaAIC)
  tmp1$rank      <- rank(tmp1$aicc)
  if(!is.null(covnm))
    tmp1$covnm   = covnm
  if(!is.null(rsType))
    tmp1$rsType   = rsType
  if(!is.null(hypoth))
    tmp1$hypoth   = hypoth
  if(!is.null(CE))
    tmp1$CE   = CE
  
  # re order from smallest AICc to largest
  tmp1           <- tmp1[order(tmp1$aicc), ]
  tmp1$AICw_std  <- tmp1$AICweight/sum(tmp1$AICweight, na.rm = T)
  tmp1$cumlAIC   <- cumsum(tmp1$AICw_std)
  
  cutoff <- which(tmp1$cumlAIC > 0.95)[1]
  if (is.na(cutoff)) 
    cutoff <- nn
  cutoff2 <- which(tmp1$deltaAIC > 2)[1]
  if (is.na(cutoff2)) 
    cutoff2 <- nn
  cutoff4 <- which(tmp1$deltaAIC > 4)[1]
  if (is.na(cutoff4)) 
    cutoff4 <- nn
  
  t1            <- rep("", nn)
  t1[1:cutoff]  <- "o"
  t2            <- rep("", nn)
  t2[1:cutoff2] <- "*"
  t4            <- rep("", nn)
  t4[1:cutoff4] <- "*"
  
 
  tmp1$topSet         <- paste(t2, t4, t1, sep = "")
  if(include_marginal){
    tmp1$aicc_marg[tmp1$aicc_marg==-Inf] <- NA
    tmp1$deltaAIC_marg  <- (tmp1$aicc_marg - min(tmp1$aicc_marg,na.rm=T))
    tmp1$AICweight_marg <- exp(-0.5 * tmp1$deltaAIC_marg)
    tmp1$rank_marg      <- rank(tmp1$aicc_marg)
    tmp1$AICw_std_marg  <- tmp1$AICweight_marg/sum(tmp1$AICweight_marg,na.rm=T)
    tmp1$AICw_std_marg[is.na(tmp1$AICw_std_marg)]<-0
    tmp1$cumlAIC_marg   <- cumsum(tmp1$AICw_std_marg)
    cutoff_marg         <- which(tmp1$cumlAIC_marg > 0.95)[1]
    
    if (is.na(cutoff_marg)) 
      cutoff_marg <- nn
    cutoff2_marg <- which(tmp1$deltaAIC_marg > 2)[1]
    if (is.na(cutoff2_marg)) 
      cutoff2_marg <- nn
    cutoff4_marg <- which(tmp1$deltaAIC_marg > 4)[1]
    if (is.na(cutoff4_marg)) 
      cutoff4_marg <- nn
    t1_marg <- rep("", nn)
    t1_marg[1:cutoff_marg] <- "o"
    t2_marg <- rep("", nn)
    t2_marg[1:cutoff2_marg] <- "*"
    t4_marg <- rep("", nn)
    t4_marg[1:cutoff4_marg] <- "*"
    tmp1$topSet_marg <- paste(t2_marg, t4_marg, t1_marg, sep = "")
    }
  
  return(tmp1)
}