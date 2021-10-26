# ----------------------------------------
# load_data.R
# subset of futR
# kirstin.holsman@noaa.gov
# updated 2020
# ----------------------------------------
  cat("loading data ...")
  
  # load ggplot theme:
  source("R/sub_scripts/THEMES_GGPLOT.r")

  #_______________________________________
  # 1.3 load data
  load(file.path(in_dir,"env_covars.Rdata"))
  load(file.path(in_dir,"rec_dat.Rdata"))
  load(file.path(in_dir,"rec_dat_tmb.Rdata"))
  load(file.path(in_dir,"ration_tmb.Rdata"))
  
  #___________________________________________
  # 2. run Rceattle() estimation mode
  #___________________________________________
  
  if( run_est ){
    cat("\n re-running R ceattle...")
    
    source("R/sub_scripts/est_Rceattle.R")
           
  }else{
             
    load(file.path(in_dir,"ss_run.Rdata" ))
    load(file.path(in_dir,"ms_run.Rdata" ))
    
  }

    # Scenarios     <-  unique(covariates$Scenario)
    # A1B_n         <-  grep("A1B",Scenarios)
    # bio_n         <-  grep("bio",Scenarios)
    # rcp45_n       <-  grep("rcp45",Scenarios)
    # rcp85_n       <-  grep("rcp85",Scenarios)
    # rcp85NoBio_n  <-  setdiff(rcp85_n,bio_n)
    # plotList      <-  Scenario_set  <- c(1,rcp45_n,rcp85NoBio_n)
    # esnm          <-  list(c(rcp45_n,rcp85NoBio_n))
    # esmlist       <-  list(rcp45_n,rcp85NoBio_n)
    # 
    # simnames  <- Scenarios
    # Years     <- sort(unique(dat_2_5_12$future_year)+start_yr-1)
    # nYrsTot   <- length(Years )
    # riskTypes <- unique(risk12$type)


# subset of downscaled projections used for the paper = Scenario_set
# bio runs are a sensitivity set of runs to evaluate nutrient forcing
# of boundary conditions, not used here bc they are highly similar to 
# non-bio runs (See Kearney et al. 2020 and Hermann et al. 2019 for more info
# A1B not used bc they were AR4 runs and only went to 2040
# print(as.character(Scenarios[Scenario_set]))
# ACLIM Projection simulations
# "###########################################################"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
# [1]  "#  | mn_Hind"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
# [2]  "#  | MIROC_A1B"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
# [3]  "#  | ECHOG_A1B"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
# [4]  "#  | CCCMA_A1B"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
# [5]  "#  | GFDL_rcp45"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
# [6]  "#  | GFDL_rcp85"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
# [7]  "#  | GFDL_rcp85_bio"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
# [8]  "#  | MIROC_rcp45"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
# [9]  "#  | MIROC_rcp85"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
# [10] "#  | CESM_rcp45"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
# [11] "#  | CESM_rcp85"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
# [12] "#  | CESM_rcp85_bio"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#  "###########################################################"  

cat("Load Data Complete")





