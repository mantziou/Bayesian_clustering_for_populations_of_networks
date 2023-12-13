# This script contains code for simulation results and visualisations presented in Section 5.2,
# as well as additional results presented in Supplementary material Section 2.2 corresponding to simulations of Section 5.2 
# Note: ordering of code in this script follows ordering of simulations appearing in Section 5.2

# First, source R scripts to run functions needed for simulating the data and running the MCMC:
setwd("/YOUR_DIRECTORY")
source("simul_fun.R")
source("vec_to_graph.R")
source("label_fun.R")
source("MCMC_mixture_e_sbm_rev.R")
library(igraph)


###########################################
### Simulation studies of Section 5.2 ### 
##########################################
# Cut-off points of accurate inference for various network sizes and sample sizes
# For C=3 clusters 
# For B=2 blocks

#function that returns representatives and their SBM membership of their nodes
gener_clust_repres<-function(w_fixed,theta_fixed,nodes,C,seed_c,seed_a,stp_seed){
  #pairs of nodes for undirected networks
  toy_mat<-matrix(nrow = nodes,ncol = nodes)
  lab<-which(upper.tri(toy_mat),arr.ind = TRUE)
  c_list<-list()
  rep_list<-lapply(1:C,function(x) vector(length = (nodes*(nodes-1))/2))
  for(i in 1:C){
    set.seed(seed_c)
    c<-rmultinom(nodes,1,w_fixed) 
    c_list[[i]]<-which(c!=0,arr.ind = TRUE)[,"row"]
    seed_c<-seed_c+stp_seed
    rep_list[[i]]<-lab_adj(rep_list[[i]],c_list[[i]],lab) #name each entry of representative by the membership of the nodes 
    rep_list[[i]]<-gener_edges_sbm(rep_list[[i]],theta_fixed,seed_a) 
    seed_a<-seed_a+stp_seed
  }
  return(list(c_list,rep_list))
}


#Generate 25-node network representatives 

#SBM paramet for 25-node representatives
lab_theta<-c("11","12","22") #corresponding to B=2
w=c(0.5,0.5)
theta=c(0.9,0.07,0.9)
names(theta)<-lab_theta

repr_25<-gener_clust_repres(w,theta,25,3,1000,10,25)

# Visualise 25-node network representatives Figure 10 in Section 5.2 of main paper, and Figure 1 in Supplement Section 2.2
a1_graph<-graph_reform(repr_25[[2]][[1]],25)
a2_graph<-graph_reform(repr_25[[2]][[2]],25)
a3_graph<-graph_reform(repr_25[[2]][[3]],25)
plot(a1_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_25[[1]][[1]]==1),which(repr_25[[1]][[1]]==2)),vertex.label.cex=2)
plot(a2_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_25[[1]][[2]]==1),which(repr_25[[1]][[2]]==2)),vertex.label.cex=2)
plot(a3_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_25[[1]][[3]]==1),which(repr_25[[1]][[3]]==2)),vertex.label.cex=2)

mtext("25-node cluster representatives",outer = TRUE,cex=1.5,line = -18)
#degree distribution
hist(degree(a1_graph))

#for the 3 25-NODE REPRESENTATIVES, Generate DIFFERENT SIZES of network populations(45,90,135,180,225,270,315)
centr_list_25<-repr_25[[2]]
c_tr_25<-repr_25[[1]]
p<-rep(0.08,3) #false positive
q<-rep(0.08,3) #false negative

## DATA SET 1 ###
seed_sim_mix<-10
for(i in seq(45,315,45)){
  assign(paste("z_tr_mix_",i,sep = ""),c(rep(1,i%/%3),rep(2,i%/%3),rep(3,i%/%3)))
  assign(paste("sim_mix_25node_",i,sep = ""),simulate_data_e(centr_list_25,p,q,i,3,seed_sim_mix))
  assign(paste("prob_vec_mix_25node_",i,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_25node_",i,sep = ""))))
  seed_sim_mix<-seed_sim_mix+100
}

## DATA SET 2,3,4,5 ### (for simulation repetitions, different seeded datasets)
ds_ind<-2 #data set index
seed_sim_mix<-20 
for(k in 1:4){
  for(i in seq(45,315,45)){
    assign(paste("sim_mix_25node_",i,"_dataset_",ds_ind,sep = ""),simulate_data_e_alter_seed(centr_list_25,p,q,i,3,seed_sim_mix))
    assign(paste("prob_vec_mix_25node_",i,"_dataset_",ds_ind,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_25node_",i,"_dataset_",ds_ind,sep = ""))))
    seed_sim_mix<-seed_sim_mix+320 
  }
  ds_ind<-ds_ind+1
}

### Simulate another 5 data sets as above with different seeds after Revisions

## DATA SET 6 ###
seed_sim_mix<-101
for(i in seq(45,315,45)){
  assign(paste("z_tr_mix_",i,sep = ""),c(rep(1,i%/%3),rep(2,i%/%3),rep(3,i%/%3)))
  assign(paste("sim_mix_25node_",i,sep = ""),simulate_data_e(centr_list_25,p,q,i,3,seed_sim_mix))
  assign(paste("prob_vec_mix_25node_",i,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_25node_",i,sep = ""))))
  seed_sim_mix<-seed_sim_mix+100
}

## DATA SET 7,8,9,10 ### 
ds_ind<-2 #data set index
seed_sim_mix<-202
for(k in 1:4){
  for(i in seq(45,315,45)){
    assign(paste("sim_mix_25node_",i,"_dataset_",ds_ind,sep = ""),simulate_data_e_alter_seed(centr_list_25,p,q,i,3,seed_sim_mix))
    assign(paste("prob_vec_mix_25node_",i,"_dataset_",ds_ind,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_25node_",i,"_dataset_",ds_ind,sep = ""))))
    seed_sim_mix<-seed_sim_mix+320 
  }
  ds_ind<-ds_ind+1
}

#Generate representatives of 50-75-100 nodes, EACH GENERATED BY DIFFERENT SBM PARAMETERS TO KEEP FIXED EXPECTED DEGREE with 25-node case
#and then generate corresponding populations

#50 node representatives with RESCALED affinity vec theta (to keep same expected degree, w kept same)
theta_50_resc<-c(22.5/50,1.75/50,22.5/50) #22.5 and 1.75 are the p_fixed resulting from multiplying 25(nodes)*0.9(theta_11 for 25-node) and 25(nodes)*0.07(theta_12 for 25-node) respectively
names(theta_50_resc)<-lab_theta

repr_50_resc<-gener_clust_repres(w,theta_50_resc,50,3,2000,20,27)

# Visualise 50-node network representatives Figure 11 in Section 5.2 of main paper, and Figure 2 in Supplement Section 2.2
a1_graph<-graph_reform(repr_50_resc[[2]][[1]],50)
a2_graph<-graph_reform(repr_50_resc[[2]][[2]],50)
a3_graph<-graph_reform(repr_50_resc[[2]][[3]],50)
plot(a1_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_50_resc[[1]][[1]]==1),which(repr_50_resc[[1]][[1]]==2)),vertex.label.cex=2)
plot(a2_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_50_resc[[1]][[2]]==1),which(repr_50_resc[[1]][[2]]==2)),vertex.label.cex=2)
plot(a3_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_50_resc[[1]][[3]]==1),which(repr_50_resc[[1]][[3]]==2)),vertex.label.cex=2)

mtext("50-node cluster representatives",outer = TRUE,cex=1.5,line = -18)
#degree distribution
hist(degree(a1_graph))

#for the 3 50-NODE REPRESENTATIVES, Generate DIFFERENT SIZES of network populations(45,90,135,180,225,270,315)
centr_list_50_resc<-repr_50_resc[[2]]
c_tr_50_resc<-repr_50_resc[[1]]
p<-rep(0.08,3) #false positive
q<-rep(0.08,3) #false negative

## DATA SET 1 ##
seed_sim_mix<-100
for(i in seq(45,315,45)){
  assign(paste("sim_mix_50node_resc_",i,sep = ""),simulate_data_e(centr_list_50_resc,p,q,i,3,seed_sim_mix))
  assign(paste("prob_vec_mix_50node_resc_",i,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_50node_resc_",i,sep = ""))))
  seed_sim_mix<-seed_sim_mix+100
}

## DATA SET 2,3,4,5 ###
ds_ind<-2 #data set index
seed_sim_mix<-8980
for(k in 1:4){
  for(i in seq(45,315,45)){
    assign(paste("sim_mix_50node_resc_",i,"_dataset_",ds_ind,sep = ""),simulate_data_e_alter_seed(centr_list_50_resc,p,q,i,3,seed_sim_mix))
    assign(paste("prob_vec_mix_50node_resc_",i,"_dataset_",ds_ind,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_50node_resc_",i,"_dataset_",ds_ind,sep = ""))))
    seed_sim_mix<-seed_sim_mix+320 
  }
  ds_ind<-ds_ind+1
}

### Simulate another 5 data sets as above with different seeds after Revisions

## DATA SET 6 ##
seed_sim_mix<-1001
for(i in seq(45,315,45)){
  assign(paste("sim_mix_50node_resc_",i,sep = ""),simulate_data_e(centr_list_50_resc,p,q,i,3,seed_sim_mix))
  assign(paste("prob_vec_mix_50node_resc_",i,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_50node_resc_",i,sep = ""))))
  seed_sim_mix<-seed_sim_mix+100
}

## DATA SET 7,8,9,10 ###
ds_ind<-2 #data set index
seed_sim_mix<-89808
for(k in 1:4){
  for(i in seq(45,315,45)){
    assign(paste("sim_mix_50node_resc_",i,"_dataset_",ds_ind,sep = ""),simulate_data_e_alter_seed(centr_list_50_resc,p,q,i,3,seed_sim_mix))
    assign(paste("prob_vec_mix_50node_resc_",i,"_dataset_",ds_ind,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_50node_resc_",i,"_dataset_",ds_ind,sep = ""))))
    seed_sim_mix<-seed_sim_mix+320 
  }
  ds_ind<-ds_ind+1
}

#75 node representatives with RESCALED affinity vec theta (to keep same expected degree, w kept same)
theta_75_resc<-c(22.5/75,1.75/75,22.5/75) #22.5 and 1.75 are the p_fixed resulting from multiplying 25(nodes)*0.9 and 25(nodes)*0.07 respectively
names(theta_75_resc)<-lab_theta

repr_75_resc<-gener_clust_repres(w,theta_75_resc,75,3,3000,30,28)

# Visualise 75-node network representatives Figure 12 in Section 5.2 of main paper, and Figure 3 in Supplement Section 2.2
library(latex2exp)
a1_graph<-graph_reform(repr_75_resc[[2]][[1]],75)
a2_graph<-graph_reform(repr_75_resc[[2]][[2]],75)
a3_graph<-graph_reform(repr_75_resc[[2]][[3]],75)
plot(a1_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_75_resc[[1]][[1]]==1),which(repr_75_resc[[1]][[1]]==2)),vertex.label.cex=2)
plot(a2_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_75_resc[[1]][[2]]==1),which(repr_75_resc[[1]][[2]]==2)),vertex.label.cex=2)
plot(a3_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_75_resc[[1]][[3]]==1),which(repr_75_resc[[1]][[3]]==2)),vertex.label.cex=2)

title(TeX(sprintf("$G^{*}_{1}")),cex.main=5,line = 0)
mtext("75-node cluster representatives",outer = TRUE,cex=1.5,line = -18)

#degree distribution
hist(degree(a1_graph))

#for the 3 75-NODE REPRESENTATIVES, Generate DIFFERENT SIZES of network populations(45,90,135,180,225,270,315)
centr_list_75_resc<-repr_75_resc[[2]]
c_tr_75_resc<-repr_75_resc[[1]]
p<-rep(0.08,3) #false positive
q<-rep(0.08,3) #false negative

## DATA SET 1 ##
seed_sim_mix<-1000
for(i in seq(45,315,45)){
  assign(paste("sim_mix_75node_resc_",i,sep = ""),simulate_data_e(centr_list_75_resc,p,q,i,3,seed_sim_mix))
  assign(paste("prob_vec_mix_75node_resc_",i,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_75node_resc_",i,sep = ""))))
  seed_sim_mix<-seed_sim_mix+100
}

## DATA SET 2,3,4,5 ###
ds_ind<-2 #data set index
seed_sim_mix<-17940
for(k in 1:4){
  for(i in seq(45,315,45)){
    assign(paste("sim_mix_75node_resc_",i,"_dataset_",ds_ind,sep = ""),simulate_data_e_alter_seed(centr_list_75_resc,p,q,i,3,seed_sim_mix))
    assign(paste("prob_vec_mix_75node_resc_",i,"_dataset_",ds_ind,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_75node_resc_",i,"_dataset_",ds_ind,sep = ""))))
    seed_sim_mix<-seed_sim_mix+320
  }
  ds_ind<-ds_ind+1
}

### Simulate another 5 data sets as above with different seeds after Revisions


## DATA SET 6 ##
seed_sim_mix<-10001
for(i in seq(45,315,45)){
  assign(paste("sim_mix_75node_resc_",i,sep = ""),simulate_data_e(centr_list_75_resc,p,q,i,3,seed_sim_mix))
  assign(paste("prob_vec_mix_75node_resc_",i,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_75node_resc_",i,sep = ""))))
  seed_sim_mix<-seed_sim_mix+100
}

## DATA SET 7,8,9,10 ###
ds_ind<-2 #data set index
seed_sim_mix<-179409
for(k in 1:4){
  for(i in seq(45,315,45)){
    assign(paste("sim_mix_75node_resc_",i,"_dataset_",ds_ind,sep = ""),simulate_data_e_alter_seed(centr_list_75_resc,p,q,i,3,seed_sim_mix))
    assign(paste("prob_vec_mix_75node_resc_",i,"_dataset_",ds_ind,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_75node_resc_",i,"_dataset_",ds_ind,sep = ""))))
    seed_sim_mix<-seed_sim_mix+320
  }
  ds_ind<-ds_ind+1
}


#100 node representatives with rescaled affinity vec theta (to keep same expected degree, w kept same)
theta_100_resc<-c(22.5/100,1.75/100,22.5/100) #22.5 and 1.75 are the p_fixed resulting from multiplying 25(nodes)*0.9 and 25(nodes)*0.07 respectively
names(theta_100_resc)<-lab_theta

repr_100_resc<-gener_clust_repres(w,theta_100_resc,100,3,4000,40,29)

# Visualise 100-node network representatives Figure 13 in Section 5.2 of main paper, and Figure 4 in Supplement Section 2.2
a1_graph<-graph_reform(repr_100_resc[[2]][[1]],100)
a2_graph<-graph_reform(repr_100_resc[[2]][[2]],100)
a3_graph<-graph_reform(repr_100_resc[[2]][[3]],100)
plot(a1_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_100_resc[[1]][[1]]==1),which(repr_100_resc[[1]][[1]]==2)),vertex.label.cex=2)
plot(a2_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_100_resc[[1]][[2]]==1),which(repr_100_resc[[1]][[2]]==2)),vertex.label.cex=2)
plot(a3_graph,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,
     vertex.label.font=2, edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(repr_100_resc[[1]][[3]]==1),which(repr_100_resc[[1]][[3]]==2)),vertex.label.cex=2)

mtext("100-node cluster representatives",outer = TRUE,cex=1.5,line = -18)

#degree distribution
hist(degree(a2_graph))

#for the 3 100-NODE REPRESENTATIVES, Generate DIFFERENT SIZES of network populations(45,90,135,180,225,270,315)
centr_list_100_resc<-repr_100_resc[[2]]
c_tr_100_resc<-repr_100_resc[[1]]
p<-rep(0.08,3) #false positive
q<-rep(0.08,3) #false negative

## DATA SET 1 ##
seed_sim_mix<-10000
for(i in seq(45,315,45)){
  assign(paste("sim_mix_100node_resc_",i,sep = ""),simulate_data_e(centr_list_100_resc,p,q,i,3,seed_sim_mix))
  assign(paste("prob_vec_mix_100node_resc_",i,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_100node_resc_",i,sep = ""))))
  seed_sim_mix<-seed_sim_mix+100
}

## DATA SET 2,3,4,5 ###
ds_ind<-2 #data set index
seed_sim_mix<-26900
for(k in 1:4){
  for(i in seq(45,315,45)){
    assign(paste("sim_mix_100node_resc_",i,"_dataset_",ds_ind,sep = ""),simulate_data_e_alter_seed(centr_list_100_resc,p,q,i,3,seed_sim_mix))
    assign(paste("prob_vec_mix_100node_resc_",i,"_dataset_",ds_ind,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_100node_resc_",i,"_dataset_",ds_ind,sep = ""))))
    seed_sim_mix<-seed_sim_mix+320
  }
  ds_ind<-ds_ind+1
}

### Simulate another 5 data sets as above with different seeds after Revisions


## DATA SET 6 ##
seed_sim_mix<-100001
for(i in seq(45,315,45)){
  assign(paste("sim_mix_100node_resc_",i,sep = ""),simulate_data_e(centr_list_100_resc,p,q,i,3,seed_sim_mix))
  assign(paste("prob_vec_mix_100node_resc_",i,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_100node_resc_",i,sep = ""))))
  seed_sim_mix<-seed_sim_mix+100
}

## DATA SET 7,8,9,10 ###
ds_ind<-2 #data set index
seed_sim_mix<-269009
for(k in 1:4){
  for(i in seq(45,315,45)){
    assign(paste("sim_mix_100node_resc_",i,"_dataset_",ds_ind,sep = ""),simulate_data_e_alter_seed(centr_list_100_resc,p,q,i,3,seed_sim_mix))
    assign(paste("prob_vec_mix_100node_resc_",i,"_dataset_",ds_ind,sep = ""),prob_vec_fun(get(paste("z_tr_mix_",i,sep = "")),get(paste("sim_mix_100node_resc_",i,"_dataset_",ds_ind,sep = ""))))
    seed_sim_mix<-seed_sim_mix+320
  }
  ds_ind<-ds_ind+1
}

# Simulations were performed in parallel on cluster, below is the tuning of MCMC for simulated data above, as well as format for running on cluster

### For each seeded data set: set-up RDATA file with arguments for main MCMC for alternative simulation regimes 
N_data_jb<-c(rep(seq(45,315,45),4))
n_nodes_jb<-c(rep(25,7),rep(50,7),rep(75,7),rep(100,7))

G_data_jb<-list()
prob_vec_jb<-list()
ind<-1
for(k in seq(25,100,25)){
  if(k==25){
    for(i in seq(45,315,45)){
      # NOTE: check names used for assigning simulate data above to variables, for each deeferent seeded data set 1-10, 
      # to correspond with names below for generating RData with MCMC arguments for each data set
      G_data_jb[[ind]]<-get(paste("sim_mix_",k,"node_",i,"_dataset_2",sep = ""))
      prob_vec_jb[[ind]]<-get(paste("prob_vec_mix_",k,"node_",i,"_dataset_2",sep = ""))
      ind<-ind+1
    }
  }else{
    for(i in seq(45,315,45)){
      G_data_jb[[ind]]<-get(paste("sim_mix_",k,"node_resc_",i,"_dataset_2",sep = ""))
      prob_vec_jb[[ind]]<-get(paste("prob_vec_mix_",k,"node_resc_",i,"_dataset_2",sep = ""))
      ind<-ind+1
    }
  }
}

# hyperparameters for Beta priors on p,q,theta
a_0_jb<-b_0_jb<-c_0_jb<-d_0_jb<-e_0_jb<-f_0_jb<-lapply(1:28, function(x) rep(0.5,3))


z_init_jb<-rep(list(z_tr_mix_45,z_tr_mix_90,z_tr_mix_135,z_tr_mix_180,z_tr_mix_225,z_tr_mix_270,z_tr_mix_315),4)
c_init_jb<-c(lapply(1:7, function(x) c_tr_25),lapply(1:7, function(x) c_tr_50_resc),lapply(1:7, function(x) c_tr_75_resc),lapply(1:7, function(x) c_tr_100_resc))

p_init_jb<-q_init_jb<-rep(0.01,28) #initial values for p,q for each job 
theta_init_jb<-rep(0.3,28)

perturb_jb<-c(lapply(1:7, function(x) c(0.01,0.001,0.02,0.005,0.03)),lapply(1:7, function(x) c(0.001,0.007,0.003,0.005,0.01)),lapply(1:7, function(x) c(0.0005,0.002,0.001,0.0008,0.005)),lapply(1:7, function(x) c(0.0005,0.0001,0.002,0.001,0.0008)))
rw_steps_jb<-lapply(1:28, function(x) c(0.01,0.005,0.001,0.0001,0.0005,0.05))

save(N_data_jb,n_nodes_jb,G_data_jb,prob_vec_jb,a_0_jb,b_0_jb,c_0_jb,d_0_jb,e_0_jb,f_0_jb,
     z_init_jb,c_init_jb,p_init_jb,q_init_jb,theta_init_jb,perturb_jb,rw_steps_jb,file="/YOUR_DIRECTORY/Tuning_sim_1_2_dataset_1.RData")

# create job index for parallel computing on cluster, for each seeded dataset, index ranging from 1 to 28 for all simulation regimes
job_ind<-as.numeric(Sys.getenv("SOME_ARRAYID"))

d<-MCMC_mix_err_sbm(3,2,500000,0,n_nodes_jb[job_ind],N_data_jb[job_ind],G_data_jb[[job_ind]],prob_vec_jb[[job_ind]],a_0_jb[[job_ind]],
                    b_0_jb[[job_ind]],c_0_jb[[job_ind]],d_0_jb[[job_ind]],e_0_jb[[job_ind]],f_0_jb[[job_ind]],
                    z_init_jb[[job_ind]],c_init_jb[[job_ind]],p_init_jb[job_ind],q_init_jb[job_ind],theta_init_jb[job_ind],
                    perturb_jb[[job_ind]],rw_steps_jb[[job_ind]])

rep_1<-as.data.frame(d[[5]][[1]])
rep_2<-as.data.frame(d[[5]][[2]])
rep_3<-as.data.frame(d[[5]][[3]])

write.fst(rep_1,path = paste("/YOUR_DIRECTORY/rep_1_job_",job_ind,".fst",sep = ""),100)
write.fst(rep_2,path = paste("/YOUR_DIRECTORY/rep_2_job_",job_ind,".fst",sep = ""),100)
write.fst(rep_3,path = paste("/YOUR_DIRECTORY/rep_3_job_",job_ind,".fst",sep = ""),100)

saveRDS(d[[1]],file = paste("/YOUR_DIRECTORY/tau_job_",job_ind,".rds",sep = ""),version = 2)
saveRDS(d[[2]],file = paste("/YOUR_DIRECTORY/z_job_",job_ind,".rds",sep = ""),version = 2)
saveRDS(d[[3]],file = paste("/YOUR_DIRECTORY/p_job_",job_ind,".rds",sep = ""),version = 2)
saveRDS(d[[4]],file = paste("/YOUR_DIRECTORY/q_job_",job_ind,".rds",sep = ""),version = 2)
saveRDS(d[[6]],file = paste("/YOUR_DIRECTORY/c_job_",job_ind,".rds",sep = ""),version = 2)
saveRDS(d[[7]],file = paste("/YOUR_DIRECTORY/w_job_",job_ind,".rds",sep = ""),version = 2)
saveRDS(d[[8]],file = paste("/YOUR_DIRECTORY/theta_job_",job_ind,".rds",sep = ""),version = 2)



#####################!!!!!!!!!!!###############################!!!!!!!!!!!!!!!!!!!!!!!########################
### Visualisations of simulation results of Section 5.2, Figures 14 and 15, for different seeded data sets ###
#####################!!!!!!!!!!!###############################!!!!!!!!!!!!!!!!!!!!!!!########################

#function that identifies posterior mode for representative
library(mgcv)
post_mode<-function(mat){
  t1<-uniquecombs(mat)
  t1_ind<-which.max(table(attr(t1,"index")))
  if(length(table(attr(t1_ind,"index")))!=1){
    post_mode_cent<-t1[t1_ind,]
  }else{
    post_mode_cent<-t1
  }
  return(post_mode_cent)
}

library(ggplot2)
library(reshape2)


### 25-node networks ###

#FIRST get matrix of post means of very first repetition corresponding to jobs for 25-node networks

p_mat<-matrix(NA_real_,nrow =70,ncol = 3 ) # 70 (rows) = 10 (data sets from diff seeds) * 7 (sample sizes N for each n), 3 columns corresponding to number of clusters
q_mat<-matrix(NA_real_,nrow =70,ncol = 3 )
tau_mat<-matrix(NA_real_,nrow =70,ncol = 3 )
ind <- 1
for (i in 1:10){
  for (j in 1:7){
    # NOTE: depending on the names of files and directories used to save results from MCMC
    p<-readRDS(paste("/YOUR_DIRECTORY/DIRECRTORY_REPETITION_",i,"/p_job_",j,".rds",sep=""))
    q<-readRDS(paste("/YOUR_DIRECTORY/DIRECRTORY_REPETITION_",i,"/q_job_",j,".rds",sep=""))
    tau<-readRDS(paste("/YOUR_DIRECTORY/DIRECRTORY_REPETITION_",i,"/tau_job_",j,".rds",sep=""))
    p_mat[ind, ]<-apply(p[-c(1:150000),],2,mean)
    q_mat[ind, ]<-apply(q[-c(1:150000),],2,mean)
    tau_mat[ind, ]<-apply(tau[-c(1:150000),],2,mean)
    ind <- ind+1
  }
}


### ABSOLUTE ERROR ###
meas_p<-abs(0.08-p_mat) 
meas_q<-abs(0.08-q_mat)
meas_tau<-abs(0.3333333-tau_mat) 

p_col<-c(meas_p[,1],meas_p[,2],meas_p[,3])#make the matrix a single column
q_col<-c(meas_q[,1],meas_q[,2],meas_q[,3])
tau_col<-c(meas_tau[,1],meas_tau[,2],meas_tau[,3])

df_25_node<-data.frame("sample_size"=rep(rep(seq(45,315,45),10),3),p=p_col,q=q_col,'tau'=tau_col)
melt_df_25_node<-melt(df_25_node,id.var="sample_size")
melt_new<-data.frame(melt_df_25_node,cluster=rep(c(rep(1,70),rep(2,70),rep(3,70)),3)) #rep as many as variables 
melt_new$sample_size<-ordered(melt_new$sample_size,levels=seq(45,315,45)) 
melt_new$cluster<-as.factor(melt_new$cluster)

### Boxplots ###
# grey scale
ggplot(melt_new,aes(x=sample_size,y=value,fill=cluster))+geom_boxplot( )+facet_grid(variable~.,scales = "free_y",labeller=label_parsed)+scale_colour_manual(values = c("grey88", "grey73", "grey57"),"cluster")+scale_fill_manual(values = c("grey88", "grey73", "grey57"),"cluster")+
  xlab("sample size")+ylab("absolute error")+scale_x_discrete(expand = c(0.09,-0.3))+stat_summary(fun="mean", geom="line", aes(group=cluster),size=0.3)+theme(text = element_text(size = 10))+theme(text = element_text(size = 10),strip.text.y = element_text(angle =0),legend.position="bottom")


### NOTE: Repeat above for 50-node, 75-node and 100-node networks ###



########## Visualise Hamming distance results for representatives Figure 5 in Supplement Section 2.2 ##########


#function for Hamming distance
Hamm_dist<-function(v1,v2){ #takes arguments vectors of upper/lower triangle of adjacency 
  return(sum(abs(v1-v2))) 
}

# Get Hamming between posterior draws for representative and true representative
# 25-node case
rep_mat_repet_dh1<-matrix(NA_real_,nrow =70,ncol = 3 )
rep_mat_repet_dh5<-matrix(NA_real_,nrow =70,ncol = 3 )
rep_mat_repet_dh10<-matrix(NA_real_,nrow =70,ncol = 3 )
ind<-1
for(i in 1:10){ # corresponding to 10 repetitions (10 different seeded data sets)
  for (j in 1:7){ # corresponding to simulation regimes (different N sample sizes) for 25-node networks
    print(c("25node iter: ",ind))
    rep_1<-read.fst(paste("/YOUR_DIRECTORY/rep_1_job_",j,".fst",sep=""),from = 150001,to=500000)
    rep_2<-read.fst(paste("/YOUR_DIRECTORY/rep_2_job_",j,".fst",sep=""),from = 150001,to=500000)
    rep_3<-read.fst(paste("/YOUR_DIRECTORY/rep_3_job_",j,".fst",sep=""),from = 150001,to=500000)
    # Hamming from true representatives
    ham_1<-apply(rep_1,1,function(x)Hamm_dist(x,centr_list_25[[1]])) 
    ham_2<-apply(rep_2,1,function(x)Hamm_dist(x,centr_list_25[[2]]))
    ham_3<-apply(rep_3,1,function(x)Hamm_dist(x,centr_list_25[[3]]))
    # Proportion of times Hamming less than 1, 5 and 10
    rep_mat_repet_dh1[ind,]<-c(length(which(ham_1<=1))/350000,length(which(ham_2<=1))/350000,length(which(ham_3<=1))/350000)
    rep_mat_repet_dh5[ind,]<-c(length(which(ham_1<=5))/350000,length(which(ham_2<=5))/350000,length(which(ham_3<=5))/350000)
    rep_mat_repet_dh10[ind,]<-c(length(which(ham_1<=10))/350000,length(which(ham_2<=10))/350000,length(which(ham_3<=10))/350000)
    ind<-ind+1
  }
}

# Note: Repeat above for 50, 75 and 100 node 

# visual for Hamming<=1

# results for 25-node
rep_25_mat_col<-c(rep_mat_repet_dh1[,1],rep_mat_repet_dh1[,2],rep_mat_repet_dh1[,3])
# results for 50-node
rep_50_mat_col<-c(rep_mat_repet_dh1[,1],rep_mat_repet_dh1[,2],rep_mat_repet_dh1[,3])
# results for 75-node
rep_75_mat_col<-c(rep_mat_repet_dh1[,1],rep_mat_repet_dh1[,2],rep_mat_repet_dh1[,3])
# results for 100-node
rep_100_mat_col<-c(rep_mat_repet_dh1[,1],rep_mat_repet_dh1[,2],rep_mat_repet_dh1[,3])

df_rep<-data.frame("sample_size"=rep(rep(seq(45,315,45),10),3),'25'=rep_25_mat_col,'50'=rep_50_mat_col,'75'=rep_75_mat_col,'100'=rep_100_mat_col,check.names = FALSE)

melt_df_rep<-melt(df_rep,id.var="sample_size")
melt_new<-data.frame(melt_df_rep,cluster=rep(c(rep(1,70),rep(2,70),rep(3,70)),4)) 
melt_new$sample_size<-ordered(melt_new$sample_size,levels=seq(45,315,45)) 
melt_new$cluster<-as.factor(melt_new$cluster)

library(ggbeeswarm)
ggplot(melt_new,aes(x=sample_size,y=value,color=cluster))+geom_beeswarm(dodge.width=0.1,size =0.3)+facet_grid(variable~.,scales = "free_y",labeller=label_parsed)+scale_colour_manual(values = c("darkgoldenrod3", "darkmagenta", "mediumseagreen"),"cluster")+
  xlab("sample size")+ylab("Proportion of Hamming<=10")+scale_x_discrete(expand = c(0.09,-0.3))+stat_summary(fun="mean", geom="line", aes(group=cluster))+theme(text = element_text(size = 12),legend.position = "bottom")+ylim(0,1)

# Note: Repeat above for Hamming<=5 and Hamming<=10

##############################
# Sensitivity analysis for hyperparameters
##############################

# variance of distribution for different parameters
var(rbeta(1000,0.5,0.5))
apply(rdirichlet(1000,rep(0.5,2)),2,var)

var(rbeta(1000,1,1))
apply(rdirichlet(1000,rep(1,2)),2,var)
(0.12-0.08)/0.12

var(rbeta(1000,1.75,1.75))
apply(rdirichlet(1000,rep(1.75,2)),2,var)
(.082-0.056)/0.082

# Run parallel on cluster
# For one of the different seeded data sets, use job 28 index (corresponding to n=100 nodes, N=315 networks)
# and run MCMC for different hyperparameter settings

job_ind<-as.numeric(Sys.getenv("SOME_ARRAYID"))

# different hyperparameter settings
hyperpar_cases<-c(0.5,1,1.75)

d<-MCMC_mix_err_sbm(3,2,500000,0,n_nodes_jb[28],N_data_jb[28],G_data_jb[[28]],prob_vec_jb[[28]],rep(hyperpar_cases[job_ind],3),
                    rep(hyperpar_cases[job_ind],3),rep(hyperpar_cases[job_ind],3),rep(hyperpar_cases[job_ind],3),rep(hyperpar_cases[job_ind],3),
                    rep(hyperpar_cases[job_ind],3),
                    z_init_jb[[28]],c_init_jb[[28]],p_init_jb[28],q_init_jb[28],theta_init_jb[28],
                    perturb_jb[[28]],rw_steps_jb[[28]])

# Save results run on cluster
rep_1<-as.data.frame(d[[5]][[1]])
rep_2<-as.data.frame(d[[5]][[2]])
rep_3<-as.data.frame(d[[5]][[3]])

write.fst(rep_1,path = paste("/YOUR_DIRECTORY/rep_1_hyperpar_",hyperpar_cases[job_ind],".fst",sep = ""),100)
write.fst(rep_2,path = paste("/YOUR_DIRECTORY/rep_2_hyperpar_",hyperpar_cases[job_ind],".fst",sep = ""),100)
write.fst(rep_3,path = paste("/YOUR_DIRECTORY/rep_3_hyperpar_",hyperpar_cases[job_ind],".fst",sep = ""),100)

saveRDS(d[[1]],file = paste("/YOUR_DIRECTORY/tau_hyperpar_",hyperpar_cases[job_ind],".rds",sep = ""),version = 2)
saveRDS(d[[2]],file = paste("/YOUR_DIRECTORY/z_hyperpar_",hyperpar_cases[job_ind],".rds",sep = ""),version = 2)
saveRDS(d[[3]],file = paste("/YOUR_DIRECTORY/p_hyperpar_",hyperpar_cases[job_ind],".rds",sep = ""),version = 2)
saveRDS(d[[4]],file = paste("/YOUR_DIRECTORY/q_hyperpar_",hyperpar_cases[job_ind],".rds",sep = ""),version = 2)
saveRDS(d[[6]],file = paste("/YOUR_DIRECTORY/c_hyperpar_",hyperpar_cases[job_ind],".rds",sep = ""),version = 2)
saveRDS(d[[7]],file = paste("/YOUR_DIRECTORY/w_hyperpar_",hyperpar_cases[job_ind],".rds",sep = ""),version = 2)
saveRDS(d[[8]],file = paste("/YOUR_DIRECTORY/theta_hyperpar_",hyperpar_cases[job_ind],".rds",sep = ""),version = 2)