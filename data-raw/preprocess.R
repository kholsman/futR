# Load raw data from .csv file

#datlist <- makefutR_data(fn=file.path("data-raw","in","futR_Inputs.xlsx"))
load("data-raw/datlist.Rdata")
rec <-  data.frame(readxl::read_xlsx(file.path("data-raw","in","futR_Inputs.xlsx") , sheet = "rec_data" ))

# Apply preprocessing...
# Save the cleaned data in the required R package location

usethis::use_data(datlist, overwrite = T)
usethis::use_data(rec, overwrite = T)