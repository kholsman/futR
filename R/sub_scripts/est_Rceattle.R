# est_Rceattle.R
# sub script of futR
# Rceattle based on Grant Adams TMB formulation of the EBS CEATTLE model
data(BS2017SS)
data(BS2017MS)

# 2.1 Run single species model est
ss_run <- Rceattle(
  data_list    = BS2017SS,
  inits        = NULL, # Initial parameters = 0
  file_name    = NULL, # Don't save
  debug        = 0, # Estimate
  random_rec   = FALSE, # No random recruitment
  msmMode      = 0, # Single species mode
  avgnMode     = 0,
  silent       = TRUE
)

# 2.2 Run multispecies model est
ms_run <- Rceattle(
  data_list   = BS2017MS,
  inits       = ss_run$estimated_params, # Initial parameters from ss run
  file_name   = NULL, # Don't save
  debug       = 0, # Estimate
  random_rec  = FALSE, # No random recruitment
  niter       = 10, # Number of iterations around predation/pop dy functions
  msmMode     = 1, # Multi-species holsman mode
  avgnMode    = 0 # Use average N
)

save( ss_run,file   =   file.path(in_dir,"ss_run.Rdata" ))
save( ms_run,file   =   file.path(in_dir,"ms_run.Rdata" ))