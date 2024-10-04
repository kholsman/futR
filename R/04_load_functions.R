#' load_functions.R
#'
#' A script to load all the functions in the subfolder
#'
#'  ----------------------------------------
#' load_functions.R
#' subset of futR
#' kirstin.holsman@noaa.gov
#' updated 2024
#' ----------------------------------------
#' 
dirlist<-dir("R")
  for(d in dirlist[grep("sub_fun_",dirlist)]) 
    source(file.path("R",d))



  
 
  
  
  

  
  
  
  