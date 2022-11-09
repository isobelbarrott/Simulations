task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
seed <- as.numeric(task_id_string)

library(survival)
library(MASS)
library(simsurv)
library(parallel)
library(riskRegression)
source("~/simulation/landmarking_functions.R")

set.seed(seed)
n<-100000
long_names<-c("SBP","TCHDL")
M<-length(long_names)

mean_TCHDL<-4.066782
sd_TCHDL<-1.227767
mean_SBP<-130.8466
sd_SBP<-15.2515
mean_response_time_F<-59.53185
sd_response_time_F<-11.0672
beta_2_SBP<-(((10/sd_SBP)/(12/sd_response_time_F))-((10/sd_SBP)/(7/sd_response_time_F)))/2
beta_2_TCHDL<-(1-((exp(1.5)-exp(1.35)/sd_TCHDL)/(10/sd_response_time_F)))/2

betaLong_SBP_intercept = c(-0.20609671)
betaLong_TCHDL_intercept = c(-0.40278582)
betaLong_SBP_linear = c(0.23890863)
betaLong_TCHDL_linear = c(0.23890863-0.27378762)
betaLong_SBP_Z = c(
  0,#european
  0.04945855,#pacific
  0.07692593,#maori
  -0.14108632,#indian
  -0.28563515,#Chinese_other_Asian
  0,#smoking0
  0.01401680,#smoking1
  0.04280505,#smoking2
  0,#diabetes0
  0.01566089,#diabetes1
  0,#NZDep1
  0.03135204,#NZDep2
  0.05550907,#NZDep3
  0.06357570,#NZDep4
  0.06625303,#NZDep5
  0,#AF0
  -0.25061181,#AF1
  0,#lipid0
  -0.14311350,#lipid1
  0,#bp0
  0.43614999,#bp1
  0,#anti0
  -0.07350341#anti1
)
betaLong_TCHDL_Z = c(
  0,#european
  0.04945855+0.10419838,#pacific
  0.07692593+0.10003816,#maori
  -0.14108632+0.33833009,#indian
  -0.28563515+0.34065963,#Chinese_other_Asian
  0,#smoking0
  0.01401680+0.01759397,#smoking1
  0.04280505+0.11716144,#smoking2
  0,#diabetes0
  0.01566089+0.11092018,#diabetes1
  0,#NZDep1
  0.03135204+0.01226029,#NZDep2
  0.05550907+0.01885694,#NZDep3
  0.06357570+0.03954735,#NZDep4
  0.06625303+0.06387897,#NZDep5
  0,#AF0
  -0.25061181+0.21241656,#AF1
  0,#lipid0
  -0.14311350-0.32387511,#lipid1
  0,#bp0
  0.43614999-0.35922017,#bp1
  0,#anti0
  -0.07350341+0.05589255#anti1
)
names(betaLong_SBP_Z)<-
  names(betaLong_TCHDL_Z)<-
  c("European","Pacific","NZMaori","Indian","Chinese_other_Asian",
    "Non_smoker","Ex_smoker","Smoker",
    "No_Diabetes","Diabetes",
    "NZDep1","NZDep2","NZDep3","NZDep4","NZDep5",
    "No_atrial_fibrillation","Atrial_fibrillation",
    "No_lipid_med","Lipid_med",
    "No_bp_med","Bp_med",
    "No_antithrombotic_med","Antithrombotic_med"
  )

b_sd = c(0.6487523,0.6948037,0.1620254,0.2143920)
b_rho = matrix(c( 1    , 0.173,-0.112,-0.133,
                  0.173, 1    ,-0.072,-0.269,
                  -0.112,-0.072, 1    , 0.191,
                  -0.133,-0.269, 0.191,  1     ),4,4,byrow=TRUE)

betaEvent_Z = c(0,0.22635722,0.39120690,0.12739327,-0.26949379,
                0,0.16521979,0.60422111,
                0,0.52870299,
                0,0.12739627,0.16064562,0.21822776,0.34931682,
                0,0.79911609,
                0,0.01261087,
                0,0.26236855,
                0,0.28918990
                -0.94760190,
                -0.47114811,
                1.05759181,
                -0.07443217,
                -0.23095285)
names(betaEvent_Z)<-
  c("European","Pacific","NZMaori","Indian","Chinese_other_Asian",
    "Non_smoker","Ex_smoker","Smoker",
    "No_Diabetes","Diabetes",
    "NZDep1","NZDep2","NZDep3","NZDep4","NZDep5",
    "No_atrial_fibrillation","Atrial_fibrillation",
    "No_lipid_med","Lipid_med",
    "No_bp_med","Bp_med",
    "No_antithrombotic_med","Antithrombotic_med",
    "SBP_slope",
    "TCHDL_slope",
    "Age_start",
    "Age_start:SBP",
    "Age_start:Diabetes"
  )

betaEvent_shape<-0.11219950
betaEvent_scale<--5.57556619
betaEvent_lambda<-exp(betaEvent_shape)
betaEvent_gamma<-exp(betaEvent_scale)

betaEvent_SBP_assoc = c(0.36482246)
betaEvent_TCHDL_assoc = c(0.15968587)

jm_hazard<-function(t, x, betas, M = 1,
                    trajectory = "linear", assoc = "etavalue",
                    family = list(gaussian()), grp_assoc = NULL){
  if (t == 0){return(0)}
  
  etavalues<-lapply(long_names,function(m){
    name_response <- paste0("Long_",m)
    etavalue <-
      betas[[name_response]][[paste0("betaLong_",m,"_fixed_intercept")]] +
      betas[[name_response]][[paste0("betaLong_",m,"_European")]] * x[["EthnicityEuropean"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Pacific")]] * x[["EthnicityPacific"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_NZMaori")]] * x[["EthnicityNZMaori"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Indian")]] * x[["EthnicityIndian"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Chinese_other_Asian")]] * x[["EthnicityChinese_other_Asian"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Non_smoker")]] * x[["SmokingNon_smoker"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Ex_smoker")]] * x[["SmokingEx_smoker"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Smoker")]] * x[["SmokingSmoker"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_No_Diabetes")]] * x[["DiabetesNo_Diabetes"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Diabetes")]] * x[["DiabetesDiabetes"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_NZDep1")]] * x[["NZDepNZDep1"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_NZDep2")]] * x[["NZDepNZDep2"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_NZDep3")]] * x[["NZDepNZDep3"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_NZDep4")]] * x[["NZDepNZDep4"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_NZDep5")]] * x[["NZDepNZDep5"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_No_atrial_fibrillation")]] * x[["Atrial_fibrillationNo_atrial_fibrillation"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Atrial_fibrillation")]] * x[["Atrial_fibrillationAtrial_fibrillation"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_No_lipid_med")]] * x[["Lipid_medNo_lipid_med"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Lipid_med")]] * x[["Lipid_medLipid_med"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_No_bp_med")]] * x[["Bp_medNo_bp_med"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Bp_med")]] * x[["Bp_medBp_med"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_No_antithrombotic_med")]] * x[["Antithrombotic_medNo_antithrombotic_med"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_Antithrombotic_med")]] * x[["Antithrombotic_medAntithrombotic_med"]]+
      betas[[name_response]][[paste0("betaLong_",m,"_fixed_linear")]] * (x[["Age_start"]] + t/sd_response_time_F) +
      betas[[name_response]][[paste0("betaLong_",m,"_random_intercept")]] +
      betas[[name_response]][[paste0("betaLong_",m,"_random_linear")]] * (x[["Age_start"]] + t/sd_response_time_F)
    etavalue
  })
  names(etavalues)<-paste0("etavalue_",long_names)
  etavalues[["etavalue_SBP"]][which(etavalues[["etavalue_SBP"]]<(0-mean_SBP)/sd_SBP)]<-(0-mean_SBP)/sd_SBP
  etavalues[["etavalue_TCHDL"]][which(etavalues[["etavalue_TCHDL"]]<(0-mean_TCHDL)/sd_TCHDL)]<-(0-mean_TCHDL)/sd_TCHDL
  
  phi_event <- exp(
    betas[["Event"]][["betaEvent_European"]] * x[["EthnicityEuropean"]] +
      betas[["Event"]][["betaEvent_Pacific"]] * x[["EthnicityPacific"]] +
      betas[["Event"]][["betaEvent_NZMaori"]] * x[["EthnicityNZMaori"]] +
      betas[["Event"]][["betaEvent_Indian"]] * x[["EthnicityIndian"]] +
      betas[["Event"]][["betaEvent_Chinese_other_Asian"]] * x[["EthnicityChinese_other_Asian"]] +
      betas[["Event"]][["betaEvent_Non_smoker"]] * x[["SmokingNon_smoker"]] +
      betas[["Event"]][["betaEvent_Ex_smoker"]] * x[["SmokingEx_smoker"]] +
      betas[["Event"]][["betaEvent_Smoker"]] * x[["SmokingSmoker"]] +
      betas[["Event"]][["betaEvent_No_Diabetes"]] * x[["DiabetesNo_Diabetes"]] +
      betas[["Event"]][["betaEvent_Diabetes"]] * x[["DiabetesDiabetes"]] +
      betas[["Event"]][["betaEvent_NZDep1"]] * x[["NZDepNZDep1"]] +
      betas[["Event"]][["betaEvent_NZDep2"]] * x[["NZDepNZDep2"]] +
      betas[["Event"]][["betaEvent_NZDep3"]] * x[["NZDepNZDep3"]] +
      betas[["Event"]][["betaEvent_NZDep4"]] * x[["NZDepNZDep4"]] +
      betas[["Event"]][["betaEvent_NZDep5"]] * x[["NZDepNZDep5"]] +
      betas[["Event"]][["betaEvent_No_atrial_fibrillation"]] * x[["Atrial_fibrillationNo_atrial_fibrillation"]] +
      betas[["Event"]][["betaEvent_Atrial_fibrillation"]] * x[["Atrial_fibrillationAtrial_fibrillation"]] +
      betas[["Event"]][["betaEvent_No_lipid_med"]] * x[["Lipid_medNo_lipid_med"]] +
      betas[["Event"]][["betaEvent_Lipid_med"]] * x[["Lipid_medLipid_med"]] +
      betas[["Event"]][["betaEvent_No_bp_med"]] * x[["Bp_medNo_bp_med"]] +
      betas[["Event"]][["betaEvent_Bp_med"]] * x[["Bp_medBp_med"]] +
      betas[["Event"]][["betaEvent_No_antithrombotic_med"]] * x[["Antithrombotic_medNo_antithrombotic_med"]] +
      betas[["Event"]][["betaEvent_Antithrombotic_med"]] * x[["Antithrombotic_medAntithrombotic_med"]]+
      betas[["Event"]][["betaEvent_SBP_assoc"]] * etavalues[["etavalue_SBP"]] +
      betas[["Event"]][["betaEvent_TCHDL_assoc"]] * etavalues[["etavalue_TCHDL"]]+
      betas[["Event"]][["betaEvent_SBP_slope"]] * betas[["Long_SBP"]][[paste0("betaLong_SBP_random_linear")]] +
      betas[["Event"]][["betaEvent_TCHDL_slope"]] * betas[["Long_TCHDL"]][[paste0("betaLong_TCHDL_random_linear")]] +
      betas[["Event"]][["betaEvent_Age_start"]] * x[["Age_start"]] +
      betas[["Event"]][["betaEvent_Age_start:SBP"]] * etavalues[["etavalue_SBP"]] * x[["Age_start"]] +
      betas[["Event"]][["betaEvent_Age_start:Diabetes"]] * x[["DiabetesDiabetes"]] * x[["Age_start"]]
  )
  
  h <- phi_event * betaEvent_lambda * betaEvent_gamma * t^(betaEvent_lambda-1)
  if (!length(h) == 1) {
    stop("Bug found: returned hazard should be a scalar.")
  }
  return(h)
}



#############
Ethnicity<-factor(sample(x=c("European","Pacific","NZMaori","Indian","Chinese_other_Asian"),
                         prob = c(117379,26882,28360,16564,22069)/211253,
                         replace = TRUE,
                         size=n),levels=c("European","Pacific","NZMaori","Indian","Chinese_other_Asian"))
Smoking<-factor(sample(x=c("Non_smoker","Ex_smoker","Smoker"),
                       prob=c(151010,36168,24075)/211253,
                       replace=TRUE,
                       size=n),levels=c("Non_smoker","Ex_smoker","Smoker"))
Diabetes<-factor(sample(x=c("No_Diabetes","Diabetes"),
                        prob=c(175835,35418)/211253,
                        replace=TRUE,
                        size=n),levels=c("No_Diabetes","Diabetes"))
NZDep<-factor(sample(x=c("NZDep1","NZDep2","NZDep3","NZDep4","NZDep5"),
                     prob=c(45065,41403,38060,38993,47732)/211253,
                     replace=TRUE,
                     size=n),levels=c("NZDep1","NZDep2","NZDep3","NZDep4","NZDep5"))
Atrial_fibrillation<-factor(sample(x=c("No_atrial_fibrillation","Atrial_fibrillation"),
                                   prob=c(207862,3391)/211253,
                                   replace=TRUE,
                                   size=n),levels=c("No_atrial_fibrillation","Atrial_fibrillation"))
Lipid_med<-factor(sample(x=c("No_lipid_med","Lipid_med"),
                         prob=c(155163,56090)/211253,
                         replace=TRUE,
                         size=n),levels=c("No_lipid_med","Lipid_med"))
Bp_med<-factor(sample(x=c("No_bp_med","Bp_med"),
                      prob=c(173342,37911)/211253,
                      replace=TRUE,
                      size=n),levels=c("No_bp_med","Bp_med"))
Antithrombotic_med<-factor(sample(x=c("No_antithrombotic_med","Antithrombotic_med"),
                                  prob=c(125516,85737)/211253,
                                  replace=TRUE,
                                  size=n),levels=c("No_antithrombotic_med","Antithrombotic_med"))
Age_start<-(runif(n,min=25, max=80)-mean_response_time_F)/sd_response_time_F


Z<-data.frame(Ethnicity,Smoking,Diabetes,NZDep,Atrial_fibrillation,Lipid_med,Bp_med,Antithrombotic_med)
Z<-model.matrix(~-1+.,data=Z,contrasts.arg = lapply(Z,contrasts,contrasts=FALSE))
covs <- data.frame(id = 1:n, Age_start, Z)

b_dd <- MASS::mvrnorm(n = n, mu = rep(0, 4), Sigma = b_rho)
b <- sapply(1:length(b_sd), function(x) b_sd[x] * b_dd[, x])
b<-as.data.frame(b)
colnames(b)<-c("SBP_intercept","TCHDL_intercept","SBP_linear","TCHDL_linear")
betas <- list()
for (m in long_names){
  name_response <- paste0("Long_", m)
  name_fe_intercept<-paste0("betaLong_",m,"_intercept")
  name_fe_linear<-paste0("betaLong_",m,"_linear")
  name_Z<-paste0("betaLong_",m,"_Z")
  betas[[name_response]] <- data.frame(id = 1:n)
  betas[[name_response]][[paste0("betaLong_",m,"_fixed_intercept")]] <- rep(get(name_fe_intercept),n)
  betas[[name_response]][[paste0("betaLong_",m,"_fixed_linear")]] <- rep(get(name_fe_linear),n)
  betas[[name_response]][[paste0("betaLong_",m,"_random_intercept")]] <-  b[[paste0(m,"_intercept")]]
  betas[[name_response]][[paste0("betaLong_",m,"_random_linear")]] <- b[[paste0(m,"_linear")]]
  for (i in names(get(name_Z))){
    betas[[name_response]][[paste0("betaLong_",m,"_",i)]] <- rep(get(name_Z)[[i]], n)
  }
}

betas[["Event"]] <- data.frame(id = 1:n)
for (m in long_names){
  name_assoc<-paste0("betaEvent_",m,"_assoc")
  betas[["Event"]][,name_assoc] <- rep(get(name_assoc), n)
  for (i in names(betaEvent_Z)){
    betas[["Event"]][[paste0("betaEvent_",i)]] <- rep(betaEvent_Z[i], n)
  }
}

interval = c(1e-08, 200)
max_fuptime = 20

ss <- simsurv::simsurv(hazard = jm_hazard, x = covs, betas = betas,
                       idvar = "id",ids = covs$id, maxt = max_fuptime, interval = interval)

dat <- purrr::reduce(list(covs, ss, betas[["Long_SBP"]], betas[["Long_TCHDL"]]),dplyr::left_join,by="id")

saveRDS(dat,file=paste0("/rds/user/ib401/hpc-work/dat_simsurv/dat_simsurv_",seed,".RDS"))
