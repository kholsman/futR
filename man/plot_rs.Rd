% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sub_fun_plot_rs.R
\name{plot_rs}
\alias{plot_rs}
\title{plot_rs}
\usage{
plot_rs(rec_df = rec_fit, cbp = "Dark2", obscol = "gray")
}
\arguments{
\item{rec_df}{dataframe of recruitment and spawners generated by getFit()}

\item{cbp}{color_brewer palette for the different models; default = "Dark2"}

\item{obscol}{color for the "observed" recruitment; default = "gray"}
}
\value{
list
}
\description{
plot the recruitment model in two panels a) R~f(S) and b) R~f(Year)
For more information contact author Kirstin Holsman (kirstin.holsman at noaa.gov)
}
\examples{
# read in the data
datlist <- datlist2 <-  readMake_futR_data("data/in/futR_Inputs.xlsx" )
mm      <- runmod(dlistIN   = datlist, version   = 'futR',recompile = FALSE,simulate  = TRUE,sim_nitr  = 1000)
r_fit  <- getFit(mm, nm = paste0("recType = ",datlist$rs_dat$rectype))
plot_rs(r_fit)

# set the recType to 2
datlist2$rs_dat$rectype <- 2
mm      <- runmod(dlistIN   = datlist2, version   = 'futR',recompile = FALSE,simulate  = TRUE,sim_nitr  = 1000)
r_fit2  <- getFit(mm, nm = paste0("recType = ",datlist2$rs_dat$rectype))

# now compare plots
plot_rs(rbind(r_fit,r_fit2))
}
