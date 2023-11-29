# ----------------------------------------
# packages.R
# load or install packages
# subset of Holsman et al. 2020 Nature Comm.
# kirstin.holsman@noaa.gov
# updated 2020
# ----------------------------------------

lib_list <- c(
  # these for reshaping and manipulating data:
    "Matrix", # needed for TMB
    "TMB",
   # "TMBhelper",
    "writexl",
    "viridis",
    "reshape",
    "tidyr",
    "dplyr", 
    "ggplot2",
    "lattice",
    "purrr",
    "readxl",
 #   "Rceattle",
    
  # markdown stuff:
    "knitr",
    "kableExtra",
    
  # These for making plots:
    "cowplot",
    "RColorBrewer",
    "ggplot2", 
    "wesanderson",
    "scales",
    "ggforce",
    "grid",
    "processx",
    "plotly",
    "extrafont",
 
 # These for Rpackages
    "roxygen2"
  )

# Install missing libraries:
missing <- setdiff(lib_list, installed.packages()[, 1])
if (length(missing) > 0) install.packages(missing,repos = "http://cran.us.r-project.org", dependencies = TRUE)

# Load libraries:
for(lib in lib_list)
       eval(parse(text=paste("library(",lib,")")))


