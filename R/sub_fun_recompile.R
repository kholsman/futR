#' recompile()
#' 
#' recompiles the model
#' 
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#'  
#' @param model (without the '.cpp'); default is "futR"
#'
#' @return list
#'
#' @examples
#' recompile(model="futR")
#' @export
#' 
recompile<-function(model="futR"){
  if(file.exists( paste0(model,'.so')) )
    file.remove( paste0(model,'.so')) 
  if(file.exists( paste0(model,'.o')) )
  file.remove( paste0(model,'.o')) 

 aa <-  TMB::compile( paste0(model,'.cpp')) 
  if(file.exists(paste0(model,'.dll')) ){
    file.copy(from =paste0(model,'.dll'),to =paste0("../",model,'.dll'))
    #file.remove(paste0(model,'.dll'))
  }
 file.copy(from =paste0(model,'.o'),to =paste0("../",model,'.o'))
 file.copy(from =paste0(model,'.so'),to =paste0("../",model,'.so'))
 if(aa != 0) stop(paste(model, "did not compile"))
 print(aa)
}
