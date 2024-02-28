#'
#'
#'
#'get_cov_shape.R
#'
  get_cov_shape <- function(
                      AIC_summryIN,
                      s = 1,
                      SD_range =3,
                      simulateIN=F,
                      sim_nitrIN = 200,
                      futR_path = "../../futR/"){
    
      tmdir <- getwd()
      setwd(futR_path)
      tryCatch({
          AIC_summry <- AIC_summryIN
          modIN      <- AIC_summry[[s]]$topRicker_R2
          recIN      <- AIC_summry[[s]]$topRicker_R2_rec_fit
          covs       <- strsplit(as.character(AIC_summry[[s]]$topRicker_R2$name),":")[[1]][2]
          covs       <- strsplit(covs,"_PLUS_")[[1]]
          covs2      <- covs[-grep("x2",covs)]
          
          for(c in 1:length(covs2)){
            nn   <- grep(covs2[c],covs)
            T2IN <- F
            if(length(nn)>1) T2IN <-T
            
            tmp <- profile_covars(
              cov_nm     = covs2[c],
              T2         = T2IN,
              SD_range   =  SD_range,
              sim_nitrIN = sim_nitrIN,
              simulateIN = simulateIN,
              modIN    = modIN,
              recIN    = recIN,
              phasesIN = phases)
            if(c==1) out <-tmp
            if(c>1)  out <- rbind(out,tmp)
            rm(tmp)
            
          }
          p <- ggplot(out)+geom_line(aes(x=cov,y=R_hat,color=name))+theme_minimal()
          return(list(p=p, out = out))
        },
        error = function(cond) {
            message("error message:")
            message(conditionMessage(cond))
            # Choose a return value in case of error
            NA
          },
        warning = function(cond) {
         
          message("warning message:")
          message(conditionMessage(cond))
          # Choose a return value in case of warning
          NULL
        },
        finally = {
         setwd(tmpdir)
        })
  }
  
  