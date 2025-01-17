#' update R package
#'
#' script to updated the Rpackage for futR
#'

# create the package folder
usethis::create_package( "/Users/KKH/Documents/GitHub_mac/futR")
usethis::use_test("test-data")
source("data-raw/preprocess.R")
