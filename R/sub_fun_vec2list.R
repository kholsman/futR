#' convert vector to list using names
#' 
#' @email For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @weblink 
#' @param x Numeric vector.
#'
#' @return list
#'
#' @examples
#' bmi.vals <- rnorm(n = 50, mean = 25, sd = 3)
#' vec2list(bmi.vals)
#'
#' @export
vec2list <- function(x){
  n  <-  unique(names(x))
  outlist<-list()
  for(nx in n)
    outlist[[nx]]<-x[names(x)%in%nx]
  return(outlist)
}

