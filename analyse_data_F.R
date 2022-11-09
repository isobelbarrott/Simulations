task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
seed <- as.numeric(task_id_string)
source("~/simulation/landmarking_functions.R")

mean_TCHDL<-4.066782
sd_TCHDL<-1.227767
mean_SBP<-130.8466
sd_SBP<-15.2515
mean_response_time_F<-59.53185
sd_response_time_F<-11.0672
obs<-40
betaLong_SBP_aux = c(1*0.7092341)
betaLong_TCHDL_aux = c(0.6850294*0.7092341)

for (x_l in c(40,50,60,70,80)){
  set.seed(seed)
  dat<-readRDS(paste0(file=paste0("/rds/user/ib401/hpc-work/dat_simsurv/dat_simsurv_",seed,".RDS")))
  m<-"SBP"
  n_row<-nrow(dat)
  name_response <- paste0("Long_", m)
  dat <- dat[rep(seq_len(n_row), each=62),]
  
  tij_seq_SBP <- rep(c(seq(x_l-10-0.125, x_l+5-0.125,length=61),x_l),n_row)
  dat$tij_SBP <- (tij_seq_SBP-mean_response_time_F)/sd_response_time_F
  
  dat[[paste0("Xij_", m)]] <-
    dat[[paste0("betaLong_",m,"_fixed_intercept")]] +
    dat[[paste0("betaLong_",m,"_random_intercept")]] +
    dat[[paste0("betaLong_",m,"_European")]] * dat[["EthnicityEuropean"]]+
    dat[[paste0("betaLong_",m,"_Pacific")]] * dat[["EthnicityPacific"]]+
    dat[[paste0("betaLong_",m,"_NZMaori")]] * dat[["EthnicityNZMaori"]]+
    dat[[paste0("betaLong_",m,"_Indian")]] * dat[["EthnicityIndian"]]+
    dat[[paste0("betaLong_",m,"_Chinese_other_Asian")]] * dat[["EthnicityChinese_other_Asian"]]+
    dat[[paste0("betaLong_",m,"_Non_smoker")]] * dat[["SmokingNon_smoker"]]+
    dat[[paste0("betaLong_",m,"_Ex_smoker")]] * dat[["SmokingEx_smoker"]]+
    dat[[paste0("betaLong_",m,"_Smoker")]] * dat[["SmokingSmoker"]]+
    dat[[paste0("betaLong_",m,"_No_Diabetes")]] * dat[["DiabetesNo_Diabetes"]]+
    dat[[paste0("betaLong_",m,"_Diabetes")]] * dat[["DiabetesDiabetes"]]+
    dat[[paste0("betaLong_",m,"_NZDep1")]] * dat[["NZDepNZDep1"]]+
    dat[[paste0("betaLong_",m,"_NZDep2")]] * dat[["NZDepNZDep2"]]+
    dat[[paste0("betaLong_",m,"_NZDep3")]] * dat[["NZDepNZDep3"]]+
    dat[[paste0("betaLong_",m,"_NZDep4")]] * dat[["NZDepNZDep4"]]+
    dat[[paste0("betaLong_",m,"_NZDep5")]] * dat[["NZDepNZDep5"]]+
    dat[[paste0("betaLong_",m,"_No_atrial_fibrillation")]] * dat[["Atrial_fibrillationNo_atrial_fibrillation"]]+
    dat[[paste0("betaLong_",m,"_Atrial_fibrillation")]] * dat[["Atrial_fibrillationAtrial_fibrillation"]]+
    dat[[paste0("betaLong_",m,"_No_lipid_med")]] * dat[["Lipid_medNo_lipid_med"]]+
    dat[[paste0("betaLong_",m,"_Lipid_med")]] * dat[["Lipid_medLipid_med"]]+
    dat[[paste0("betaLong_",m,"_No_bp_med")]] * dat[["Bp_medNo_bp_med"]]+
    dat[[paste0("betaLong_",m,"_Bp_med")]] * dat[["Bp_medBp_med"]]+
    dat[[paste0("betaLong_",m,"_No_antithrombotic_med")]] * dat[["Antithrombotic_medNo_antithrombotic_med"]]+
    dat[[paste0("betaLong_",m,"_Antithrombotic_med")]] * dat[["Antithrombotic_medAntithrombotic_med"]]+
    dat[[paste0("betaLong_",m,"_fixed_linear")]] * dat[["tij_SBP"]]+
    dat[[paste0("betaLong_",m,"_random_linear")]] * dat[["tij_SBP"]]
  
  mu <- dat[[paste0("Xij_", m)]]
  sigma <- get(paste0("betaLong_",m,"_aux"))
  dat[[paste0("Yij_", m)]] <- stats::rnorm(length(mu), mu, sigma)
  
  dat$tij_SBP <- tij_seq_SBP
  m<-"TCHDL"
  name_response <- paste0("Long_", m)
  
  dat$tij_TCHDL <- tij_seq_SBP
  dat$tij_TCHDL <- (dat$tij_TCHDL-mean_response_time_F)/sd_response_time_F
  
  dat[[paste0("Xij_", m)]] <-
    dat[[paste0("betaLong_",m,"_fixed_intercept")]] +
    dat[[paste0("betaLong_",m,"_random_intercept")]] +
    dat[[paste0("betaLong_",m,"_European")]] * dat[["EthnicityEuropean"]]+
    dat[[paste0("betaLong_",m,"_Pacific")]] * dat[["EthnicityPacific"]]+
    dat[[paste0("betaLong_",m,"_NZMaori")]] * dat[["EthnicityNZMaori"]]+
    dat[[paste0("betaLong_",m,"_Indian")]] * dat[["EthnicityIndian"]]+
    dat[[paste0("betaLong_",m,"_Chinese_other_Asian")]] * dat[["EthnicityChinese_other_Asian"]]+
    dat[[paste0("betaLong_",m,"_Non_smoker")]] * dat[["SmokingNon_smoker"]]+
    dat[[paste0("betaLong_",m,"_Ex_smoker")]] * dat[["SmokingEx_smoker"]]+
    dat[[paste0("betaLong_",m,"_Smoker")]] * dat[["SmokingSmoker"]]+
    dat[[paste0("betaLong_",m,"_No_Diabetes")]] * dat[["DiabetesNo_Diabetes"]]+
    dat[[paste0("betaLong_",m,"_Diabetes")]] * dat[["DiabetesDiabetes"]]+
    dat[[paste0("betaLong_",m,"_NZDep1")]] * dat[["NZDepNZDep1"]]+
    dat[[paste0("betaLong_",m,"_NZDep2")]] * dat[["NZDepNZDep2"]]+
    dat[[paste0("betaLong_",m,"_NZDep3")]] * dat[["NZDepNZDep3"]]+
    dat[[paste0("betaLong_",m,"_NZDep4")]] * dat[["NZDepNZDep4"]]+
    dat[[paste0("betaLong_",m,"_NZDep5")]] * dat[["NZDepNZDep5"]]+
    dat[[paste0("betaLong_",m,"_No_atrial_fibrillation")]] * dat[["Atrial_fibrillationNo_atrial_fibrillation"]]+
    dat[[paste0("betaLong_",m,"_Atrial_fibrillation")]] * dat[["Atrial_fibrillationAtrial_fibrillation"]]+
    dat[[paste0("betaLong_",m,"_No_lipid_med")]] * dat[["Lipid_medNo_lipid_med"]]+
    dat[[paste0("betaLong_",m,"_Lipid_med")]] * dat[["Lipid_medLipid_med"]]+
    dat[[paste0("betaLong_",m,"_No_bp_med")]] * dat[["Bp_medNo_bp_med"]]+
    dat[[paste0("betaLong_",m,"_Bp_med")]] * dat[["Bp_medBp_med"]]+
    dat[[paste0("betaLong_",m,"_No_antithrombotic_med")]] * dat[["Antithrombotic_medNo_antithrombotic_med"]]+
    dat[[paste0("betaLong_",m,"_Antithrombotic_med")]] * dat[["Antithrombotic_medAntithrombotic_med"]]+
    dat[[paste0("betaLong_",m,"_fixed_linear")]] * dat[["tij_TCHDL"]]+
    dat[[paste0("betaLong_",m,"_random_linear")]] * dat[["tij_TCHDL"]]
  
  mu <- dat[[paste0("Xij_", m)]]
  sigma <- get(paste0("betaLong_",m,"_aux"))
  dat[[paste0("Yij_", m)]] <- stats::rnorm(length(mu), mu, sigma)
  
  dat$tij_TCHDL <- tij_seq_SBP
  dat$Age_start<-dat$Age_start*sd_response_time_F+mean_response_time_F
  dat$eventtime<-dat$Age_start+dat$eventtime
  
  Ethnicity<-dat[,grep("^Ethnicity",colnames(dat))]
  Ethnicity<-apply(Ethnicity,1,function(row){colnames(Ethnicity)[which.max(row)]})
  Ethnicity<-factor(gsub(pattern="^Ethnicity",replacement = "",x = Ethnicity),
                    levels=c("European","Pacific","NZMaori","Indian","Chinese_other_Asian"))
  Smoking<-dat[,grep("^Smoking",colnames(dat))]
  Smoking<-apply(Smoking,1,function(row){colnames(Smoking)[which.max(row)]})
  Smoking<-factor(gsub(pattern="^Smoking",replacement = "",x = Smoking),
                  levels=c("Non_smoker","Ex_smoker","Smoker"))
  Diabetes<-dat[,grep("^Diabetes",colnames(dat))]
  Diabetes<-apply(Diabetes,1,function(row){colnames(Diabetes)[which.max(row)]})
  Diabetes<-factor(gsub(pattern="^Diabetes",replacement = "",x = Diabetes),
                   levels=c("No_Diabetes","Diabetes"))
  NZDep<-dat[,grep("^NZDep",colnames(dat))]
  NZDep<-apply(NZDep,1,function(row){colnames(NZDep)[which.max(row)]})
  NZDep<-factor(gsub(pattern="^NZDep",replacement = "",x = NZDep),
                levels=c("NZDep1","NZDep2","NZDep3","NZDep4","NZDep5"))
  Atrial_fibrillation<-dat[,grep("^Atrial_fibrillation",colnames(dat))]
  Atrial_fibrillation<-apply(Atrial_fibrillation,1,function(row){colnames(Atrial_fibrillation)[which.max(row)]})
  Atrial_fibrillation<-factor(gsub(pattern="^Atrial_fibrillation",replacement = "",x = Atrial_fibrillation),
                              levels=c("No_atrial_fibrillation","Atrial_fibrillation"))
  Lipid_med<-dat[,grep("^Lipid_med",colnames(dat))]
  Lipid_med<-apply(Lipid_med,1,function(row){colnames(Lipid_med)[which.max(row)]})
  Lipid_med<-factor(gsub(pattern="^Lipid_med",replacement = "",x = Lipid_med),levels=c("No_lipid_med","Lipid_med"))
  Bp_med<-dat[,grep("^Bp_med",colnames(dat))]
  Bp_med<-apply(Bp_med,1,function(row){colnames(Bp_med)[which.max(row)]})
  Bp_med<-factor(gsub(pattern="^Bp_med",replacement = "",x = Bp_med),levels=c("No_bp_med","Bp_med"))
  Antithrombotic_med<-dat[,grep("^Antithrombotic_med",colnames(dat))]
  Antithrombotic_med<-apply(Antithrombotic_med,1,function(row){colnames(Antithrombotic_med)[which.max(row)]})
  Antithrombotic_med<-factor(gsub(pattern="^Antithrombotic_med",replacement = "",x = Antithrombotic_med),levels=c("No_antithrombotic_med","Antithrombotic_med"))
  
  dat[[paste0("SBP_slope")]]<-dat[[paste0("betaLong_SBP_random_linear")]]+dat[[paste0("betaLong_SBP_fixed_linear")]]
  dat[[paste0("TCHDL_slope")]]<-dat[[paste0("betaLong_TCHDL_random_linear")]]+dat[[paste0("betaLong_TCHDL_fixed_linear")]]
  
  dat<-data.frame(dat[,c("id","eventtime", "status", "tij_SBP", "Yij_SBP", "Xij_SBP", "tij_TCHDL", "Yij_TCHDL","Xij_TCHDL","Age_start","SBP_slope","TCHDL_slope")],
                  Ethnicity,
                  Smoking,
                  Diabetes,
                  NZDep,
                  Atrial_fibrillation,
                  Lipid_med,
                  Bp_med,
                  Antithrombotic_med
  )
  
  dat[["Yij_SBP"]][which(dat[["Yij_SBP"]]<(0-mean_SBP)/sd_SBP)]<-(0-mean_SBP)/sd_SBP
  dat[["Yij_TCHDL"]][which(dat[["Yij_TCHDL"]]<(0-mean_TCHDL)/sd_TCHDL)]<-(0-mean_TCHDL)/sd_TCHDL
  dat[["Xij_SBP"]][which(dat[["Xij_SBP"]]<(0-mean_SBP)/sd_SBP)]<-(0-mean_SBP)/sd_SBP
  dat[["Xij_TCHDL"]][which(dat[["Xij_TCHDL"]]<(0-mean_TCHDL)/sd_TCHDL)]<-(0-mean_TCHDL)/sd_TCHDL
  
  dat<-dat[dat$Age_start<=x_l,]
  dat <-
    return_ids_with_LOCF(
      data_long = dat,
      individual_id = "id",
      covariates=c("Ethnicity","Smoking","Diabetes","NZDep","Atrial_fibrillation","Lipid_med","Bp_med","Antithrombotic_med","Yij_SBP", "Yij_TCHDL"),
      covariates_time=c(rep("tij_SBP",9), "tij_TCHDL"),
      x_L = x_l
    )
  dat<-dat[dat$eventtime>x_l,]
  dat<-dat[dat$eventtime>=dat$tij_SBP,]
  
  dat_LME_LOCF<-dat[dat$tij_SBP %in% seq(x_l-10-0.125, x_l+5-0.125,by=10/obs),]
  
   dat_LOCF_F_x_l<-fit_LOCF_landmark(data_long=dat_LME_LOCF,
                                     x_L=x_l,
                                     x_hor=x_l+5,
                                     covariates=c("Ethnicity","Smoking","Diabetes","NZDep",
                                                  "Atrial_fibrillation","Lipid_med","Bp_med","Antithrombotic_med","Yij_SBP", "Yij_TCHDL"),
                                     covariates_time=c(rep("tij_SBP",9), "tij_TCHDL"),
                                     individual_id="id",
                                     event_time="eventtime",
                                     event_status="status",
                                     survival_submodel=c("standard_cox")
   )
   saveRDS(dat_LOCF_F_x_l[[1]]$model_survival$coef,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LOCF_coef/dat_LOCF_F_coef_",x_l,"_",obs,"_",seed,".RDS"))
   saveRDS(dat_LOCF_F_x_l[[1]]$prediction_error,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LOCF_prediction_error/dat_LOCF_F_",x_l,"_",obs,"_",seed,".RDS"))
   saveRDS(dat_LOCF_F_x_l[[1]],file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LOCF_F/dat_LOCF_F_",x_l,"_",obs,"/dat_LOCF_F_",x_l,"_",obs,"_",seed,".RDS"))
   rm(dat_LOCF_F_x_l)

  dat_gold<-dat[dat$tij_SBP == x_l,]
  dat_gold_F_x_l<-fit_LOCF_landmark(data_long=dat_gold,
                                    x_L=x_l,
                                    x_hor=x_l+5,
                                    covariates=c("Ethnicity","Smoking","Diabetes","NZDep",
                                                 "Atrial_fibrillation","Lipid_med","Bp_med","Antithrombotic_med","Xij_SBP", "Xij_TCHDL","SBP_slope","TCHDL_slope"),
                                    covariates_time=c(rep("tij_SBP",9),"tij_TCHDL","tij_SBP","tij_TCHDL"),
                                    individual_id="id",
                                    event_time="eventtime",
                                    event_status="status",
                                    survival_submodel=c("standard_cox")
  )
  saveRDS(dat_gold_F_x_l[[1]]$model_survival$coef,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_gold_with_slope_coef/dat_gold_F_coef_",x_l,"_",obs,"_",seed,".RDS"))
  saveRDS(dat_gold_F_x_l[[1]]$prediction_error,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_gold_with_slope_prediction_error/dat_gold_F_",x_l,"_",obs,"_",seed,".RDS"))
  saveRDS(dat_gold_F_x_l[[1]],file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_gold_F_with_slope/dat_gold_F_",x_l,"/dat_gold_F_",x_l,"_",obs,"_",seed,".RDS"))
  rm(dat_gold_F_x_l)

  dat_gold_F_x_l<-fit_LOCF_landmark(data_long=dat_gold,
                                    x_L=x_l,
                                    x_hor=x_l+5,
                                    covariates=c("Ethnicity","Smoking","Diabetes","NZDep",
                                                 "Atrial_fibrillation","Lipid_med","Bp_med","Antithrombotic_med","Xij_SBP", "Xij_TCHDL"),
                                    covariates_time=c(rep("tij_SBP",9),"tij_TCHDL"),
                                    individual_id="id",
                                    event_time="eventtime",
                                    event_status="status",
                                    survival_submodel=c("standard_cox")
  )

  saveRDS(dat_gold_F_x_l[[1]]$model_survival$coef,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_gold_coef/dat_gold_F_coef_",x_l,"_",obs,"_",seed,".RDS"))
  saveRDS(dat_gold_F_x_l[[1]]$prediction_error,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_gold_prediction_error/dat_gold_F_",x_l,"_",obs,"_",seed,".RDS"))
  saveRDS(dat_gold_F_x_l[[1]],file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_gold_F/dat_gold_F_",x_l,"/dat_gold_F_",x_l,"_",obs,"_",seed,".RDS"))
  rm(dat_gold_F_x_l)

  dat_LME_F_x_l<-fit_LME_landmark(data_long=dat_LME_LOCF,
                                    x_L=x_l,
                                    x_hor=x_l+5,
                                    fixed_effects = c("Ethnicity","Smoking","Diabetes","NZDep",
                                                      "Atrial_fibrillation","Lipid_med","Bp_med","Antithrombotic_med"),
                                    random_effects = c("Yij_SBP", "Yij_TCHDL"),
                                    fixed_effects_time=rep("tij_SBP", 8),
                                    random_effects_time=c("tij_SBP", "tij_TCHDL"),
                                    individual_id="id",
                                    event_time="eventtime",
                                    event_status="status",
                                    survival_submodel=c("standard_cox"),
                                    lme_control = nlme::lmeControl(maxIter = 1e8, msMaxIter = 1e8,tolerance = 1e-4,returnObject = TRUE),
                                    standardise_time=TRUE,
                                    standardise_time_mean=mean_response_time_F,
                                    standardise_time_sd=sd_response_time_F,
                                    random_slope_survival = TRUE
   )
   saveRDS(dat_LME_F_x_l[[1]]$model_survival$coef,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LME_with_slope_coef/dat_LME_F_coef_",x_l,"_",obs,"_",seed,".RDS"))
   saveRDS(dat_LME_F_x_l[[1]]$prediction_error,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LME_with_slope_prediction_error/dat_LME_F_",x_l,"_",obs,"_",seed,".RDS"))
   saveRDS(dat_LME_F_x_l[[1]],file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LME_F_with_slope/dat_LME_F_",x_l,"_",obs,"/dat_LME_F_",x_l,"_",obs,"_",seed,".RDS"))
   rm(dat_LME_F_x_l)

  
  dat_LME_F_x_l<-fit_LME_landmark_no_slope(data_long=dat_LME_LOCF,
                                           x_L=x_l,
                                           x_hor=x_l+5,
                                           fixed_effects = c("Ethnicity","Smoking","Diabetes","NZDep",
                                                             "Atrial_fibrillation","Lipid_med","Bp_med","Antithrombotic_med"),
                                           random_effects = c("Yij_SBP", "Yij_TCHDL"),
                                           fixed_effects_time=rep("tij_SBP", 8),
                                           random_effects_time=c("tij_SBP", "tij_TCHDL"),
                                           individual_id="id",
                                           event_time="eventtime",
                                           event_status="status",
                                           survival_submodel=c("standard_cox"),
                                           lme_control = nlme::lmeControl(maxIter = 1e8, msMaxIter = 1e8,tolerance = 1e-4,returnObject = TRUE),
                                           standardise_time=TRUE,
                                           standardise_time_mean=mean_response_time_F,
                                           standardise_time_sd=sd_response_time_F,
                                           random_slope_survival = FALSE,x_L_no_slope=x_l,seed_no_slope=seed,obs_no_slope=obs
  )
  saveRDS(dat_LME_F_x_l[[1]]$model_survival$coef,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LME_coef/dat_LME_F_coef_",x_l,"_",obs,"_",seed,".RDS"))
  saveRDS(dat_LME_F_x_l[[1]]$prediction_error,file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LME_prediction_error/dat_LME_F_",x_l,"_",obs,"_",seed,".RDS"))
  saveRDS(dat_LME_F_x_l[[1]],file=paste0("/rds/project/rds-csoP2nj6Y6Y/ib401/dat_LME_F/dat_LME_F_",x_l,"_",obs,"/dat_LME_F_",x_l,"_",obs,"_",seed,".RDS"))
  rm(dat_LME_F_x_l)
  
}
