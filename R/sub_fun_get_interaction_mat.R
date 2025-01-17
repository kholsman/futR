#' get_interaction_mat.R
#'
#' Function to get the interaction matrix for env covariates
#'
#'@param covIN    covars (matrix with cols for each covar) e.g., pre_spawning_covs[,-1]
#'@param maxIN    max number of interactions
#'@param ADDTEMP2 add a temperature ^2 term to those covars with "temp" in the name
#'@param cor_cutoff  cut off value for correlation matrix (e.g., 0.5), those with cor > than this value will be excluded
#' get full set of 4 interactions (max_interactions)
#' returns a list with data.frames of all possible combinations
#' @export
#' @importFrom stats cor
#' @importFrom utils combn
#'
#' @examples
#' #get_interaction_mat()
#'
get_interaction_mat<-function(covIN,
                              maxIN,
                              ADDTEMP2=TRUE,
                              cor_cutoff=0.5){
  mod_mat <- list()
  cor_mat <- abs(cor(covIN,covIN))
  for(i in 1:maxIN){
    tt <- combn(covIN,i,simplify=F)
    for(k in 1:length(tt)){
      if(i>1){
        ex <- cor(tt[[k]])
        diag(ex) <- 0
        ex[abs(ex)>cor_cutoff] <- 0
        # }else{
        #   ex = 0
        # }
        names(tt)[k] <-  paste0(names(tt[[k]]),collapse="_PLUS_")
        # if the covars are not related:
        if(  length(which(as.numeric(ex[1,])!=0)) ==(i-1) ){
          mod_mat <- c(mod_mat,tt[k])
        }
        rm(ex)
      }
      if(ADDTEMP2){
        aa <- grep("temp",colnames(tt[[k]]))
        if(length(aa)>0){
          # add temperature squared
          tmp           <- data.frame(tt[[k]][,aa]^2)
          colnames(tmp) <- paste0(colnames(tt[[k]])[aa],"x2")
          ll            <- list(data.frame(  tt[[k]],tmp))
          names(ll)     <- paste0(names(ll[[1]]),collapse="_PLUS_")
          mod_mat       <- c(mod_mat,ll)
          rm(tmp)
        }
        rm(aa)
      }
    }
    rm(tt)
  }
  return(mod_mat)
}
