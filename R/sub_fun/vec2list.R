#-------------------------------------
# convert vector to list using names
#-------------------------------------

vec2list <- function(x){
  n  <-  unique(names(x))
  outlist<-list()
  for(nx in n)
    outlist[[nx]]<-x[names(x)%in%nx]
  return(outlist)
}

