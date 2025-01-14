#' load_functions.R
#'
#' A script to load all the functions in the subfolder
#'
#'  ----------------------------------------
#' load_functions.R
#' subset of futR
#' kirstin.holsman at noaa.gov
#' updated 2024
#' ----------------------------------------
#' 

#load the functions in R folder
#if( length(setdiff("futR", installed.packages()[, 1]))>0 ){
if(1==10){
dirlist<-dir("R")
  for(d in dirlist[grep("sub_fun_",dirlist)])
    source(file.path("R",d))
}
#}

#load_all()
  
 
  
  
  

  
  
  
  