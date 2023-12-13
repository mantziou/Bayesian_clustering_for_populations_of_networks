# This script contains code for simulations presented in Section 5.3

# source R script with main code for MCMC with SFM extension
setwd("/YOUR_DIRECTORY")
source("MCMC_SFM.R")
source("vec_to_graph.R")

###########################################################################
# ------Tuning and initialisation for Sparse Finite Mixture-------------- #
###########################################################################


# Load simulated data from Section 5.1 (Simulations_Section_5.1.R script)
load("/YOUR_DIRECTORY/Sim_1_1_data.RData")

# function to obtain Jaccard distance between networks
jaccard_mat<-function(adj_lst){
  jaccard_mt<-matrix(rep(0,length(adj_lst)*length(adj_lst)) ,nrow = length(adj_lst),ncol = length(adj_lst))
  for (i in 1:length(adj_lst)){
    for (j in i:length(adj_lst)){
      jaccard_mt[i,j]<-sum(abs(adj_lst[[i]]-adj_lst[[j]]))/(sum(pmax(adj_lst[[i]],adj_lst[[j]])))
      jaccard_mt[j,i]<-jaccard_mt[i,j]
    }
  }
  return(jaccard_mt)
}




# get lists of adjacencies from lists of vectors for regimes 1,2,7,8 of Section 5.1
adj_list_1 <- lapply(Sim_mix_1, function(x) adjac_reform_undir(x,21))
adj_list_2 <- lapply(Sim_mix_2, function(x) adjac_reform_undir(x,21))
adj_list_7 <- lapply(Sim_mix_7, function(x) adjac_reform_undir(x,21))
adj_list_8 <- lapply(Sim_mix_8, function(x) adjac_reform_undir(x,21))

# get jaccard distances for each regime
regimes<-c(1,2,7,8)
for (i in regimes){
  assign(paste("jacc_d_",i,sep = ""),jaccard_mat(get(paste("adj_list_",i,sep=""))))
}


############# Tuning ###################
############## Cmax=10 ################
#########################################

#### Regime 1 ####

# Get membership for networks for init z_init using jaccard and kmemoid
library(kmed)
set.seed(12)
z_init_1<-fastkmed(jacc_d_1,ncluster = 10,iterate = 1000)$cluster

prob_vec_mix_gen_1<-prob_vec_fun_gen(z_init_1,Sim_mix_1,10)

rw_steps_e0<-0.001
e0_init<-.2
a_0<-b_0<-c_0<-d_0<-chi_0<-e_0<-f_0<-rep(0.5,10)
a_e=1
b_e=400
p_init<-q_init<-0.01
theta_init<-.3
pert<-c(0.01,0.001,0.02,0.005,0.03,0.008)
rw_steps<-c(0.01,0.005,0.001,0.0001,0.0005,0.05)

cc_1<-SBM_est(prob_vec_mix_gen_1,10,12,21,2)
c_init_1<-cc_1[[1]]
G_m_init_1<-cc_1[[2]]
C_max<-10
# seed_1<-11 # seed was not used on cluster where run in parallel to other regimes
# 
# set.seed(seed_1)
sim_mix_1_sfm<-MCMC_mix_err_sbm_SFM(C_max,2,500000,0,21,180,Sim_mix_1,prob_vec_mix_gen_1,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,
                                    z_init_1,c_init_1,p_init,q_init,theta_init,e0_init,G_m_init_1,pert,rw_steps,rw_steps_e0)
# save tuning for regime 1
save(C_max,Sim_mix_1,prob_vec_mix_gen_1,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,z_init_1,c_init_1,
     p_init,q_init,theta_init,e0_init,G_m_init_1,pert,rw_steps,rw_steps_e0,seed_1,file="/YOUR_DIRECTORY/Tuning_SFM_reg_1.RData")

#### Regime 2 ####

# Get membership for networks for init z_init using jaccard and kmemoid
library(kmed)
set.seed(12)
z_init_2<-fastkmed(jacc_d_2,ncluster = 10,iterate = 1000)$cluster

prob_vec_mix_gen_2<-prob_vec_fun_gen(z_init_2,Sim_mix_2,10)

cc_2<-SBM_est(prob_vec_mix_gen_2,10,12,21,2)
c_init_2<-cc_2[[1]]
G_m_init_2<-cc_2[[2]]
# seed_2<-1111  # seed was not used on cluster where run in parallel to other regimes
# set.seed(seed_2) 
sim_mix_2_sfm<-MCMC_mix_err_sbm_SFM(C_max,2,500000,0,21,180,Sim_mix_2,prob_vec_mix_gen_2,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,
                                    z_init_2,c_init_2,p_init,q_init,theta_init,e0_init,G_m_init_2,pert,rw_steps,rw_steps_e0)
# save tuning for regime 2
save(C_max,Sim_mix_2,prob_vec_mix_gen_2,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,z_init_2,c_init_2,
     p_init,q_init,theta_init,e0_init,G_m_init_2,pert,rw_steps,rw_steps_e0,seed_2,file="/YOUR_DIRECTORY/Tuning_SFM_reg_2.RData")

#### Regime 7 ####

# Get membership for networks for init z_init using jaccard and kmemoid
library(kmed)
set.seed(12)
z_init_7<-fastkmed(jacc_d_7,ncluster = 10,iterate = 1000)$cluster

prob_vec_mix_gen_7<-prob_vec_fun_gen(z_init_7,Sim_mix_7,10)

cc_7<-SBM_est(prob_vec_mix_gen_7,10,12,21,2)
c_init_7<-cc_7[[1]]
G_m_init_7<-cc_7[[2]]
# seed_7<-11 # seed was not used on cluster where run in parallel to other regimes
# set.seed(seed_7)
sim_mix_7_sfm<-MCMC_mix_err_sbm_SFM(C_max,2,500000,0,21,180,Sim_mix_7,prob_vec_mix_gen_7,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,
                                    z_init_7,c_init_7,p_init,q_init,theta_init,e0_init,G_m_init_7,pert,rw_steps,rw_steps_e0)
# save tuning for regime 7
save(C_max,Sim_mix_7,prob_vec_mix_gen_7,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,z_init_7,c_init_7,
     p_init,q_init,theta_init,e0_init,G_m_init_7,pert,rw_steps,rw_steps_e0,seed_7,file="/YOUR_DIRECTORY/Tuning_SFM_reg_7.RData")

#### Regime 8 ####

# Get membership for networks for init z_init using jaccard and kmemoid
library(kmed)
set.seed(12)
z_init_8<-fastkmed(jacc_d_8,ncluster = 10,iterate = 1000)$cluster

prob_vec_mix_gen_8<-prob_vec_fun_gen(z_init_8,Sim_mix_8,10)

cc_8<-SBM_est(prob_vec_mix_gen_8,10,12,21,2)
c_init_8<-cc_8[[1]]
G_m_init_8<-cc_8[[2]]
# seed_8<-1111  # seed was not used on cluster where run in parallel to other regimes
# set.seed(seed_8)
sim_mix_8_sfm<-MCMC_mix_err_sbm_SFM(C_max,2,500000,0,21,180,Sim_mix_8,prob_vec_mix_gen_8,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,
                                    z_init_8,c_init_8,p_init,q_init,theta_init,e0_init,G_m_init_8,pert,rw_steps,rw_steps_e0)
# save tuning for regime 8
save(C_max,Sim_mix_8,prob_vec_mix_gen_8,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,z_init_8,c_init_8,
     p_init,q_init,theta_init,e0_init,G_m_init_8,pert,rw_steps,rw_steps_e0,seed_8,file="/YOUR_DIRECTORY/Tuning_SFM_reg_8.RData")

# NOTE: MCMC run on cluster in parallel

# Get entropy and purity indices for z

library(NMF) # library for entropy and purity functions
z_tr_mix<-c(rep(1,60),rep(2,60),rep(3,60))


mean_entr_pur<-function(post_draws, true_lab){
  entr<-c()
  pur<-c()
  for (i in 1:nrow(post_draws)){
    print(i)
    entr[i]<-entropy(as.factor(post_draws[i,]),as.factor(true_lab))
    pur[i]<-purity(as.factor(post_draws[i,]),as.factor(true_lab))
  }
  return(c(mean(entr),mean(pur)))
}

cl_entr_summ<-c()
cl_pur_summ<-c()
regimes<-c(1,2,7,8)
ind<-1
for (i in regimes){
  # load results from MCMC
  sfm<-readRDS(file = paste("/YOUR_DIRECTORY/z_job_",i,".rds",sep = ""))
  sfm_lagg<-sfm[-c(1:150000),]
  sfm_lagg<-sfm_lagg[seq(1,350000,50),]
  # get mean entropy/purity
  ep<-mean_entr_pur(sfm_lagg,z_tr_mix)
  cl_entr_summ[ind]<-ep[1]
  cl_pur_summ[ind]<-ep[2]
  ind<-ind+1
}
