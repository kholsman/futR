#'
#'recompile
#'
#'
recompile<-function(model="futR"){
  
  if(file.exists( paste0(model,'.so')) )
    file.remove( paste0(model,'.so')) 
  if(file.exists( paste0(model,'.o')) )
  file.remove( paste0(model,'.o')) 

  TMB::compile( paste0(model,'.cpp')) 
}
