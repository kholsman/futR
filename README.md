R futR
------

futR() is a generic Rpackage for fitting recruitment models to stock assessment estimates of spawning stock biomass and recruitment with or without climate covariates. The recruitment model is based on Template Model Builder (`TMB`) and formultations follow Maunder and Desriso (2011) using a generalized three parameter stock-recruitment model with environmental covariates (Deriso 1980; Schnute 1985). This includes Ricker (logistic), Beverton Holt, log-linear, and log-linear with biomass lagged by year 'y-1'. The model can be fit with and with out random effects on spawning stock biomass (SSB) and recruitment (R) (i.e., measurement error on SSB and rec) using the methods of Porch and Lauretta (2016) and with the optional unbiased estimate of sigma (sensu Ludwid and Walters 1981, Porch and Lauretta 2016). Environmental covariates are optional but can be included as main effects or as interactions.

For more information see Holsman et al. 2020 Climate and trophic controls on groundfish recruitment in Alaska.

Options
-------

### Recruitment formulations (rectype):

1.  Linear (gamma = 0)
2.  Beverton Holt (gamma = -1)
3.  Ricker (0 &lt; gamma &lt;1 ) ; gamma is estimated (tMethod):
    *a. *link = cloglog
    *b. *link = logit
4.  Exponential (gamma=1, b&lt;0)

### observation error options (sigMethod):

1.  No observation error (tau = 0)
2.  estimate sigma, random effects on SSB if tau &gt;0, tau input
3.  unbiased sigma estimate, tau input
4.  as in 1 but with defined measurement error for rec (indep of random effects on Spawners/SSB)
5.  as in 1 but with defined measurement error for rec and Spawners/SSB)

Installing futR()
-----------------

The package can be installed from github using the devtools package:

``` r
install.packages("devtools")
```

The projection package can then be installed to R directly:

``` r
devtools::install_github("kholsman/futR")
```

Fitting Recruitment:
--------------------

The base function for fitting recruitment requires a data.frame of recruitment and spawning biomass:

``` r
library("futR")
```
