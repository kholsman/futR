# ----------------------------------------
# make.R
# subset of Holsman et al. 2022 
# kirstin.holsman@noaa.gov
# updated 2022
# ----------------------------------------
cat("-- loading packages ...\n")
source("R/02_packages.R")       # loads packages
cat("-- loading setup ...\n")
source("R/03_setup.R")          # load other switches and controls
cat("-- loading functions ...\n")
source("R/04_load_functions.R") # load all functions in the R/sub_fun folder
cat("-- loading data ...\n")
source("R/05_load_data.R")      # load data
make_flag   <- TRUE

#source("R/plan2.R")          # creates the drake plan
cat("-- compiling the model ...\n\n")
wd0 <- getwd()
setwd("src/TMB")
recompile('futR')
setwd(wd0)
# options(clustermq.scheduler = "multicore") # optional parallel computing. Also needs parallelism = "clustermq"
# make(
#   plan2, # defined in R/plan.R
#   verbose = 2
# )

cat("-- make.R is Complete\n")