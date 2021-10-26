#-------------------------------------
# makeMap
#-------------------------------------

#' makeMap will create a map list for reading into the futR function
#'
#' This function prepares the estimated data for futR
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @weblink 
#' @param parameters  parameters for the tmb function
#' @param estpar      a vector of integers indicating phase. NA results in no estimation.
#' @export            maplist a list with map values
#' @email             For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @examples
#' 
#' mod<-estRec(dataINUSE  =  dataIN_admb)
#' 
makeMap<-function(param=parameters,estpar=estparams){
  mapseq  <- 1:length(unlist(param))
  ixx     <- 0
  maplist <- param
  for(inm in names(param)){
    if(inm %in% names(estpar)){
      if(estpar[[inm]]){
        
        ixx             <-  1:length(param[[inm]])+rev(ixx)[1]
        maplist[[inm]]  <- factor(ixx)
        
      }else{
        maplist[[inm]]  <- factor(rep(NA,length(maplist[[inm]])))
      }
    }else{
      maplist[[inm]]    <- factor(rep(NA,length(maplist[[inm]])))
    }
  }
  return(maplist)
  
}
