# GET_ESTM2

GET_ESTM2 <- function(mode, s = 1){
  plotdat     <-  smry_RS[[s]]
  Rectxt      <-  Rec_2$TopR2.txt[[s]]
  Rectxt_aic  <-  Rec_2$TopAIC.txt[[s]]
  rec_path    <-  rec_fldr
  estnm       <-  file.path(mod_fldr,"results/ceattle_est.std")
  rsdat       <-  rs_data
  covv        <-  covuse.all
  return(list(covv = covv,  Rectxt = Rectxt, Rectxt_aic = Rectxt_aic,rsdat=rsdat,estnm=estnm,rec_path=rec_path,plotdat=plotdat))
}