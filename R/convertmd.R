## convert overview doc to read me file
oldr<-getwd()
setwd("/Users/kholsman/GitHub_new/futR/")
library(rmarkdown)
render("futR_description.Rmd", md_document(variant = "markdown_github"),params=F)
file.copy(from="futR_description.md",to="README.md",overwrite=T)
setwd(oldr)
