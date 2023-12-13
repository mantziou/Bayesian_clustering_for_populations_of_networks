# This script contains code for simulation results and visualisations presented in Section 5.1,
# as well as additional results presented in Supplementary material Section 2.1 corresponding to simulations of Section 5.1 
# Note: ordering of code in this script follows ordering of simulation results appearing in Section 5.1

# First, source R scripts to run functions needed for simulating the data and running the MCMC:
setwd("/YOUR_DIRECTORY")
source("simul_fun.R")
source("label_fun.R")
source("MCMC_mixture_e_sbm_rev.R")
source("vec_to_graph.R")
library(jsonlite)
library(LaplacesDemon)
library(xtable)  
library(latex2exp)

###########################################
 ### Simulation studies of Section 5.1 ### 
##########################################

# Simulate the network populations below:
# Data size: n=21 nodes, N=90 networks, C=3 clusters, B=2 blocks

# get block labels for pairs of nodes for undirected networks and thetas
toy_mat<-matrix(nrow = 21,ncol = 21)
lab<-which(upper.tri(toy_mat),arr.ind = TRUE)
lab_theta<-c("11","12","22")

# Generate representatives for each SBM structure 1 and 2

#function that creates representatives for each sim study of different SBM structures
create_repres<-function(memb,thet,labl,sd){#memb:list of c_1,c_2,c_3 thet:list of theta_1,... sd:seed for Bernoulli
  a<-lapply(1:3,function(x) vector(length = 210)) #create list of empty vectors corresponding to each representative
  for(i in 1:3){
    a[[i]]<-lab_adj(a[[i]],memb[[i]],labl) 
    a[[i]]<-gener_edges_sbm(a[[i]],thet[[i]],sd) 
    sd<-sd+5000
  }
  return(a)
}

# SBM 1 structure for representatives
c_1<-c(rep(1,5),rep(2,5),rep(1,5),rep(2,6))#each representative different block membership of their nodes
c_2<-c(rep(2,10),rep(1,11))
c_3<-c(rep(2,5),rep(1,10),rep(2,6))
c_study_1<-list(c_1,c_2,c_3)
w_study_1<-rep(list(c(0.5,0.5)),3) #corresponding to w_1,w_2,w_3
theta_study_1<-rep(list(c(0.8,0.2,0.8)),3) ##corresponding to theta_1,theta_2,theta_3
names(theta_study_1[[1]])<-lab_theta
names(theta_study_1[[2]])<-lab_theta
names(theta_study_1[[3]])<-lab_theta

repres_study_1<-create_repres(c_study_1,theta_study_1,lab,5000)

# SBM 2 structure for representatives
c_1<-c(rep(1,7),rep(2,3),rep(1,8),rep(2,3))#each representative different block membership of their nodes
c_2<-c(rep(2,10),rep(1,11))
c_3<-c(rep(2,7),rep(1,6),rep(2,8))
c_study_2<-list(c_1,c_2,c_3)
w_study_2<-c(list(c(0.7,0.3)),list(c(0.5,0.5)),list(c(0.3,0.7))) #corresponding to w_1,w_2,w_3
theta_study_2<-rep(list(c(0.7,0.05,0.8)),3) ##corresponding to theta_1,theta_2,theta_3
names(theta_study_2[[1]])<-lab_theta
names(theta_study_2[[2]])<-lab_theta
names(theta_study_2[[3]])<-lab_theta

repres_study_2<-create_repres(c_study_2,theta_study_2,lab,15000)


### Visualize SBM representatives for SBM 1 and SBM 2 structures, Figure 3 in paper, Section 5.1
for(i in 1:3){
  assign(paste("graph_study_1_c",i,sep = ""),graph_reform(repres_study_1[[i]],21))
  assign(paste("graph_study_2_c",i,sep = ""),graph_reform(repres_study_2[[i]],21))
}
par(mfrow=c(1,3))
plot(graph_study_1_c1,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,vertex.label.font=2,edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(c_study_1[[1]]==1),which(c_study_1[[1]]==2)),vertex.label.cex=1.8)
plot(graph_study_1_c2,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,vertex.label.font=2,edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(c_study_1[[2]]==1),which(c_study_1[[2]]==2)),vertex.label.cex=1.8)
plot(graph_study_1_c3,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,vertex.label.font=2,edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(c_study_1[[3]]==1),which(c_study_1[[3]]==2)),vertex.label.cex=1.8)
mtext("Representatives under SBM structure 1",outer = TRUE,cex=1.5,line = -18)
plot(graph_study_2_c1,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,vertex.label.font=2,edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(c_study_2[[1]]==1),which(c_study_2[[1]]==2)),vertex.label.cex=1.8)
plot(graph_study_2_c2,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,vertex.label.font=2,edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(c_study_2[[2]]==1),which(c_study_2[[2]]==2)),vertex.label.cex=1.8)
plot(graph_study_2_c3,vertex.size=15,edge.color="slategrey", vertex.color="lightblue", vertex.label.color="black",edge.curved=0.7,vertex.label.font=2,edge.arrow.size=1.5,vertex.label.degree=3.7, main=NULL,cex.sub=15,mark.groups = list(which(c_study_2[[3]]==1),which(c_study_2[[3]]==2)),vertex.label.cex=1.8)
mtext("Representatives under SBM structure 2",outer = TRUE,cex=1.5,line = -18)


# Tuning of MCMC for implementing it on above simulated data

repres<-list(repres_study_1,repres_study_2) 
c<-list(c_study_1,c_study_2)

#Tuning of hyperparameters of Beta priors for p_c,q_c
a_0=rep(0.5,3) #hyperparameter of Beta prior, same for all clusters
b_0=rep(0.5,3)

c_0=rep(0.5,3) #hyperparameter of Beta prior, same for all clusters
d_0=rep(0.5,3)

e_0=rep(0.5,3) #hyperparameter of Beta prior, same for all clusters
f_0=rep(0.5,3)

pert<-c(0.01,0.001,0.02,0.005,0.03)
rw_stps<-c(0.01,0.005,0.001,0.0001,0.0005,0.05)

# Generate the regimes and simulate the data
t<-seq(0.1,0.3,0.1)
#sim_1_1<-list()
data_seed<-1
regime_ind<-1
for(i in t){
  for(k in t[t!=i]){
    for (j in 1:2){
      p<-rep(i,3)
      q<-rep(k,3)
      # seed used to simulate data
      Sim_mix<-simulate_data_e(repres[[j]],p,q,180,3,data_seed)
      z_tr_mix<-c(rep(1,60),rep(2,60),rep(3,60))
      prob_vec_mix<-prob_vec_fun(z_tr_mix,Sim_mix)
      res<-MCMC_mix_err_sbm(3,2,500000,0,21,180,Sim_mix,prob_vec_mix,a_0,b_0,c_0,d_0,e_0,f_0,z_tr_mix,c[[j]],0.01,0.01,0.3,pert,rw_stps)
      #sim_1_1<-list.append(sim_1_1,res)
      rep_1<-as.data.frame(res[[5]][[1]])
      rep_2<-as.data.frame(res[[5]][[2]])
      rep_3<-as.data.frame(res[[5]][[3]])
      # Save results in your preferred directory
      # write.fst(rep_1,path = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/rep_1.fst",sep = ""),100)
      # write.fst(rep_2,path = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/rep_2.fst",sep = ""),100)
      # write.fst(rep_3,path = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/rep_3.fst",sep = ""),100)
      # 
      # saveRDS(res[[1]],file = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/tau.rds",sep = ""),version = 2)
      # saveRDS(res[[2]],file = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/z.rds",sep = ""),version = 2)
      # saveRDS(res[[3]],file = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/p.rds",sep = ""),version = 2)
      # saveRDS(res[[4]],file = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/q.rds",sep = ""),version = 2)
      # saveRDS(res[[6]],file = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/c.rds",sep = ""),version = 2)
      # saveRDS(res[[7]],file = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/w.rds",sep = ""),version = 2)
      # saveRDS(res[[8]],file = paste("/YOUR_DIRECTORY/regime_",regime_ind,"/theta.rds",sep = ""),version = 2)
      
      data_seed<-data_seed+50000
      regime_ind<-regime_ind+1
      print(c(i,k,j))
    }
  }
}


# Visualise violin plots for regimes 4 and 7, Figures 4 and 5, Section 5.1

# Regime 4
par(mfrow=c(1,2))
cluster_label<-rep(c("1","2","3"),each=350000)
post_p<-readRDS(paste("/YOUR_DIRECTORY/regime_",4,"/p.rds",sep=""))
distr<-c(post_p[-c(1:150000),1],post_p[-c(1:150000),2],post_p[-c(1:150000),3])
data=data.frame(cluster_label,distr)
ggplot(data, aes( y=distr, x=cluster_label)) + geom_violin() + geom_hline(yintercept=.1, linetype="dashed",color = "red") + 
  stat_summary(fun=mean,geom = "pointrange",fun.min = function(x) p.interval(x)[1,1], fun.max = function(x) p.interval(x)[1,2] ) +theme(axis.text.x = element_text(size=11),axis.text.y = element_text(size=11),
                                                                                                                                        axis.title.y=element_text(size=13),axis.title.x=element_text(size=13))+
  labs(x="c", y = "p")

post_q<-readRDS(paste("/YOUR_DIRECTORY/regime_",4,"/q.rds",sep=""))
distr<-c(post_q[-c(1:150000),1],post_q[-c(1:150000),2],post_q[-c(1:150000),3])
data=data.frame(cluster_label,distr)
ggplot(data, aes( y=distr, x=cluster_label)) + geom_violin() + geom_hline(yintercept=.3, linetype="dashed",color = "red") + 
  stat_summary(fun=mean,geom = "pointrange",fun.min = function(x) p.interval(x)[1,1], fun.max = function(x) p.interval(x)[1,2] ) +theme(axis.text.x = element_text(size=11),axis.text.y = element_text(size=11),
                                                                                                                                        axis.title.y=element_text(size=13),axis.title.x=element_text(size=13))+
  labs(x="c", y = "q")

# Regime 7

cluster_label<-rep(c("1","2","3"),each=350000)
post_p<-readRDS(paste("/YOUR_DIRECTORY/regime_",7,"/p.rds",sep=""))
distr<-c(post_p[-c(1:150000),1],post_p[-c(1:150000),2],post_p[-c(1:150000),3])
data=data.frame(cluster_label,distr)
ggplot(data, aes( y=distr, x=cluster_label)) + geom_violin() + geom_hline(yintercept=.2, linetype="dashed",color = "red") + 
  stat_summary(fun=mean,geom = "pointrange",fun.min = function(x) p.interval(x)[1,1], fun.max = function(x) p.interval(x)[1,2] ) +theme(axis.text.x = element_text(size=11),axis.text.y = element_text(size=11),
                                                                                                                                        axis.title.y=element_text(size=13),axis.title.x=element_text(size=13))+
  labs(x="c", y = "p")

post_q<-readRDS(paste("/YOUR_DIRECTORY/regime_",7,"/q.rds",sep=""))
distr<-c(post_q[-c(1:150000),1],post_q[-c(1:150000),2],post_q[-c(1:150000),3])
data=data.frame(cluster_label,distr)
ggplot(data, aes( y=distr, x=cluster_label)) + geom_violin() + geom_hline(yintercept=.3, linetype="dashed",color = "red") + 
  stat_summary(fun=mean,geom = "pointrange",fun.min = function(x) p.interval(x)[1,1], fun.max = function(x) p.interval(x)[1,2] ) +theme(axis.text.x = element_text(size=11),axis.text.y = element_text(size=11),
                                                                                                                                        axis.title.y=element_text(size=13),axis.title.x=element_text(size=13))+
  labs(x="c", y = "q")

# Create tables for posterior means and credible intervals for p,q,theta,w presented in Supplementary Material Section 2.1, Tables 1-15

pair <- function(v) sprintf("(%.2f,%.2f)", v[1],v[2]) #function that creates pairs to store in one position of table
triple<-function(v) sprintf("(%.2f,%.2f,%.2f)",v[1],v[2],v[3]) # %.xf floating point number with x digits after decimal point
post_means_p<-matrix(nrow = 12,ncol =1)
post_means_q<-matrix(nrow = 12,ncol = 1)
cred_int_p<-matrix(nrow = 12,ncol = 3)
cred_int_q<-matrix(nrow = 12,ncol = 3)
post_means_theta_c1<-matrix(nrow = 12,ncol =1)
post_means_theta_c2<-matrix(nrow = 12,ncol = 1)
post_means_theta_c3<-matrix(nrow = 12,ncol = 1)
post_means_w_c1<-matrix(nrow = 12,ncol =1)
post_means_w_c2<-matrix(nrow = 12,ncol = 1)
post_means_w_c3<-matrix(nrow = 12,ncol = 1)
cred_int_theta_c1<-matrix(nrow = 12,ncol = 3)
cred_int_theta_c2<-matrix(nrow = 12,ncol = 3)
cred_int_theta_c3<-matrix(nrow = 12,ncol = 3)
for (i in 1:12){
  post_means_p[i]<-triple(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/p.rds",sep=""))[-c(1:150000), ],2,mean))
  post_means_q[i]<-triple(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/q.rds",sep=""))[-c(1:150000), ],2,mean))
  post_means_theta_c1[i]<-triple(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/theta.rds",sep=""))[[1]][-c(1:150000), ],2,mean))
  post_means_theta_c2[i]<-triple(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/theta.rds",sep=""))[[2]][-c(1:150000), ],2,mean))
  post_means_theta_c3[i]<-triple(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/theta.rds",sep=""))[[3]][-c(1:150000), ],2,mean))
  post_means_w_c1[i]<-pair(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/w.rds",sep=""))[[1]][-c(1:150000), ],2,mean))
  post_means_w_c2[i]<-pair(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/w.rds",sep=""))[[2]][-c(1:150000), ],2,mean))
  post_means_w_c3[i]<-pair(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/w.rds",sep=""))[[3]][-c(1:150000), ],2,mean))
  for (j in 1:3){
    cred_int_p[i,j]<-pair(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/p.rds",sep=""))[-c(1:150000), ],2,p.interval)[ ,j])
    cred_int_q[i,j]<-pair(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/q.rds",sep=""))[-c(1:150000), ],2,p.interval)[ ,j])
    cred_int_theta_c1[i,j]<-pair(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/theta.rds",sep=""))[[1]][-c(1:150000), ],2,p.interval)[ ,j])
    cred_int_theta_c2[i,j]<-pair(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/theta.rds",sep=""))[[2]][-c(1:150000), ],2,p.interval)[ ,j])
    cred_int_theta_c3[i,j]<-pair(apply(readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/theta.rds",sep=""))[[3]][-c(1:150000), ],2,p.interval)[ ,j])
  }
}

#generate latex table code for matrix in R
xtable(post_means_p,type="latex")
xtable(dist_mat,type="latex")


### Tables for Hammming , entropy, purity presented in Supplementary Material Section 2.1, Tables 16, 17, 18

library(fst)
dist_mat<-matrix(nrow = 12,ncol = 3) #rows correspond to simulation regimes and columns to distance <=1 / <=5 / <=10
clust_mat<-matrix(nrow = 12,ncol = 2) #rows corr to sim regimes and col to clust entr/clust purity
for (i in 1:12){ #i corresponds to sim regime
  d_1<-rep(0,3)
  d_5<-rep(0,3)
  d_10<-rep(0,3)
  # obtain mean clustering entropy/purity for each sim regime
  lagg_z<-readRDS(paste("/YOUR_DIRECTORY/regime_",i,"/z.rds",sep=""))[-c(1:150000), ]
  lagg_z<-lagg_z[seq(1,350000,50), ]
  clust_mat[i,1]<-mean(Clustering_Entr_alt(lagg_z,z_tr_mix))
  clust_mat[i,2]<-mean(Clustering_Purity_alt(lagg_z,z_tr_mix))
  
  for (j in 1:3){ #j correpsonds to cluster
    lagg<-read.fst(paste("/YOUR_DIRECTORY/regime_",i,"/rep_",j,".fst",sep=""),from = 150001,to=500000)
    lagg<-lagg[seq(1,350000,50), ]
    if (i %in% seq(1,12,2)){
      tot_d<-apply(lagg,1,function(x)Hamm_dist(x,repres_study_1[[j]]))
      d_1[j]<-length(which(tot_d<=1))
      d_5[j]<-length(which(tot_d<=5))
      d_10[j]<-length(which(tot_d<=10))
    }else{
      tot_d<-apply(lagg,1,function(x)Hamm_dist(x,repres_study_2[[j]]))
      d_1[j]<-length(which(tot_d<=1))
      d_5[j]<-length(which(tot_d<=5))
      d_10[j]<-length(which(tot_d<=10))
    }
  }
  dist_mat[i,1]<-triple(d_1/7000)
  dist_mat[i,2]<-triple(d_5/7000)
  dist_mat[i,3]<-triple(d_10/7000)
}

#####################################################################################
#### Implement our model on Durante, Dunson and Vogelstein (2017) simulated data #### 
#####################################################################################

# load Durante, Dunson and Vogelstein (2017) data
load("/YOUR_DIRECTORY/simulated_data.RData")

#function that vectorizes adjacency matrices, excluding diagonal for UNdirected networks
mat_vect_undir<-function(mat){
  return(mat[which(upper.tri(mat),arr.ind = TRUE)])
}

# Make list of lists the array of adjacencies, and vectorise them to use for our model
adj_list<-list()
vec_adj_list<-list()
for (i in 1:dim(A)[3]){
  vec_adj_list[[i]]<-mat_vect_undir(A[,,i])
  adj_list[[i]]<-A[,,i]
}

# Tuning our MCMC to implement it on  Durante, Dunson and Vogelstein (2017) data
C<-4
a_0<-b_0<-c_0<-d_0<-chi_0<-e_0<-f_0<-rep(0.5,4)
pert<-c(0.01,0.001,0.02,0.005,0.03)
rw_steps<-c(0.01,0.005,0.001,0.0001,0.0005,0.05)

# Initialise z_init cluster membership of networks by taking Jaccard distance metric between networks
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

jacc_dur<-jaccard_mat(adj_list)
library(kmed) # use k means to get initial cluster membership for MCMC using Jaccard distance matrix
set.seed(12)
z_init<-fastkmed(jacc_dur,ncluster = 4,iterate = 1000)$cluster

prob_vec_mix_gen<-prob_vec_fun_gen(z_init,vec_adj_list,4)

##SBM estimation of network representative estimates
library("blockmodels")
library("ramify")
SBM_est<-function(prob_vec,C_max,seed,n_nodes,K){
  s<-seed
  c_init_list<-list()
  G_m_init<-list()
  for (i in 1:C_max){
    print(i)
    set.seed(s)
    centr_est<-rbinom((n_nodes*(n_nodes-1))/2,1,prob_vec[[i]])
    G_m_init[[i]]<-centr_est
    centr_est_adj<-adjac_reform_undir(centr_est,n_nodes)
    sbm_1<-BM_bernoulli("SBM",centr_est_adj)
    sbm_1$estimate()
    #K ist most probable number of classes with: K<-which.max(sbm_1$ICL) which is 1 for this representative however I set it to 2
    c_init_list[[i]]<-argmax(sbm_1$memberships[[K]]$Z)
    s<-s+1
  }
  return(list(c_init_list,G_m_init))
}

cc<-SBM_est(prob_vec_mix_gen,4,12,20,2)
c_init<-cc[[1]]
G_m_init<-cc[[2]]

save(vec_adj_list,prob_vec_mix_gen,a_0,b_0,c_0,d_0,e_0,f_0,z_init,c_init,pert,rw_steps,file = "/YOUR_DIRECTORY/Tuning_durante_data_mixerr.RData")
# Run MCMC on  Durante, Dunson and Vogelstein (2017) data
dur_data_mix_err<-MCMC_mix_err_sbm(4,2,500000,0,20,100,vec_adj_list,prob_vec_mix_gen,a_0,b_0,c_0,d_0,e_0,f_0,z_init,c_init,0.01,0.01,.3,pert,rw_steps)

dur_data_z<-readRDS("/YOUR_DIRECTORY/dur_data_z.rds")
memb_prop<-matrix(NA,nrow = 100,ncol = 4)
for(i in 1:100){
  te<-tabulate(dur_data_z[-c(1:150000),i]) #prop of time individ i in "1" and in "2"
  if(length(te)<4){te<-c(te,rep(0,4-length(te)))}
  #print(te)
  memb_prop[i,]<-te
}
memb_prop_df<-as.data.frame(memb_prop)
names(memb_prop_df)<-c("1","2","3","4")
library(ggplot2)
library(reshape2)
melted_df<-melt(memb_prop_df)
melted_df<-data.frame(melted_df,rep(seq(1,100,1),4))
names(melted_df)<-c("Cluster","value","Individual")
ggplot(melted_df, aes(fill=Cluster, y=value/350000, x=Individual)) + 
  geom_bar(position="stack", stat="identity")+ylab("proportion")+
  theme(axis.title.x = element_text(size=30),axis.title.y = element_text( size=30),axis.text.x = element_text(size = 25),axis.text.y = element_text(size = 25),legend.text=element_text(size=30),legend.title = element_text(size = 30),legend.position = "bottom")+
  scale_x_continuous(expand = c(0,0.05),breaks=seq(1, 100, by = 5))+scale_y_continuous(expand = c(0,0.01))


########################################################
###### High noise p=q=.4 and Smaller sample sizes ######
########################################################

# nested data sets of networks

# First generate population of 180 networks with high noise p=q=0.4
data_seed <- 151001
Sim_mix_new_highnoise_8 <- simulate_data_e(repres_study_1,rep(.4,3),rep(.4,3),180,3,data_seed)

# Then keep subpopulations of networks (nested data), tune and run MCMC: 
ss <- seq(36,180,36)
count <- 1
seed <- 100
for (i in 1:length(ss)){
  print(c("regime",i))
  no_data <- ss[count]%/%3
  Sim_temp <- c(Sim_mix_new_highnoise_8[1:no_data],Sim_mix_new_highnoise_8[61:(60+no_data)],Sim_mix_new_highnoise_8[121:(120+no_data)])
  print(c("sample_size ",length(Sim_temp)))
  # jaccard initialisation
  adj_list <- list()
  for (j in 1:length(Sim_temp)){
    adj_list[[j]] <- adjac_reform_undir(Sim_temp[[j]],21)
  }
  jacc_res<-jaccard_mat(adj_list)
  # initialise using either Jaccard or informed init
  set.seed(12) 
  z_init<-fastkmed(jacc_res,ncluster = 3,iterate = 1000)$cluster
  set.seed(seed)
  #num <- length(Sim_temp)/3
  #z_init <- c(rep(1,num),rep(2,num),rep(3,num))
  prob_vec_mix<-prob_vec_fun(z_init,Sim_temp)
  cc<-SBM_est(prob_vec_mix,3,12,21,2)
  c_init<-cc[[1]]
  G_m_init<-cc[[2]]
  a_0<-b_0<-c_0<-d_0<-chi_0<-e_0<-f_0<-rep(0.5,3)
  pert<-c(0.01,0.001,0.02,0.005,0.03)
  rw_steps<-c(0.01,0.005,0.001,0.0001,0.0005,0.05)
  res <- MCMC_mix_err_sbm(3,2,500000,0,21,length(Sim_temp),Sim_temp,prob_vec_mix,a_0,b_0,c_0,d_0,e_0,f_0,z_init,c_init,0.01,0.01,.3,pert,rw_steps)
  saveRDS(res,file=paste("/YOUR_DIRECTORY/res_sim_mix_new_highnoise_",i,"_jaccinit.rds",sep = ""),version=2)
  seed <- seed+10
  count <- count+1
}



### Visualise results from our MCMC with boxplots for entropy/purity, Figures 6,7, Section 5.1

library(NMF)
ent_highnoise_list <- vector("list",5)  
pur_highnoise_list <- vector("list",5)
count <- 1
ss <- seq(36,180,36)
for (reg in c(1,2,3,4,5)){
  # read saved results from initialisation with jaccard or informed
  cutoff <- readRDS(file=paste("/YOUR_DIRECTORY/res_sim_mix_new_highnoise_",reg,"_truez.rds",sep = ""))
  #cutoff <- readRDS(file=paste("/YOUR_DIRECTORY/res_sim_mix_new_highnoise_",reg,"_jaccinit.rds",sep = ""))
  print(reg)
  num <- ss[count]/3
  z_true <- c(rep(1,num),rep(2,num),rep(3,num))
  for (i in seq(150001,500000,87)){# step 50
    print(i)
    ent_highnoise_list[[count]] <- c(ent_highnoise_list[[count]],entropy(as.factor(cutoff[[2]][i,]),as.factor(z_true)))
    pur_highnoise_list[[count]] <- c(pur_highnoise_list[[count]],purity(as.factor(cutoff[[2]][i,]),as.factor(z_true)))
  }
  count <- count+1
}

ent_highnoise_df4 <- do.call(cbind.data.frame, ent_highnoise_list)
pur_highnoise_df4 <- do.call(cbind.data.frame, pur_highnoise_list)
colnames(ent_highnoise_df4) <- c("36","72","108","144","180")
colnames(pur_highnoise_df4) <- c("36","72","108","144","180")
library(reshape2)
library(ggplot2)
ent_highnoise_df_melt4 <- melt(ent_highnoise_df4)
pur_highnoise_df_melt4 <- melt(pur_highnoise_df4)
ggplot(ent_highnoise_df_melt4, aes(x = variable, y = value)) +           
  geom_boxplot()+ggtitle("noise p=q=0.4")+ylab("entropy")
ggplot(pur_highnoise_df_melt4, aes(x = variable, y = value)) +           
  geom_boxplot()+ggtitle("noise p=q=0.4")+ylab("purity")

########################################################################
############ Different noise levels, fixed sample size N=180 ############
########################################################################

# Tuning of MCMC
#Tuning of hyperparameters of Beta priors for p_c,q_c
a_0=rep(0.5,3) #hyperparameter of Beta prior, same for all clusters
b_0=rep(0.5,3)

c_0=rep(0.5,3) #hyperparameter of Beta prior, same for all clusters
d_0=rep(0.5,3)

e_0=rep(0.5,3) #hyperparameter of Beta prior, same for all clusters
f_0=rep(0.5,3)

pert<-c(0.01,0.001,0.02,0.005,0.03)
rw_stps<-c(0.01,0.005,0.001,0.0001,0.0005,0.05)

# Run MCMC:
# Fix q, range p
p_regimes<-c(0.01,seq(0.05,0.45,0.05))
data_seed<-1

for(i in p_regimes){
  q<-rep(0.1,3)
  p<-rep(i,3)
  Sim_mix<-simulate_data_e(repres_study_1,p,q,180,3,data_seed)
  z_tr_mix<-c(rep(1,60),rep(2,60),rep(3,60))
  prob_vec_mix<-prob_vec_fun(z_tr_mix,Sim_mix)
  res<-MCMC_mix_err_sbm(3,2,500000,0,21,180,Sim_mix,prob_vec_mix,a_0,b_0,c_0,d_0,e_0,f_0,z_tr_mix,c_study_1,0.01,0.01,0.3,pert,rw_stps)
  
  saveRDS(res[[1]],file = paste("/YOUR_DIRECTORY/tau_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[2]],file = paste("/YOUR_DIRECTORY/z_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[3]],file = paste("/YOUR_DIRECTORY/p_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[4]],file = paste("/YOUR_DIRECTORY/q_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[6]],file = paste("/YOUR_DIRECTORY/c_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[7]],file = paste("/YOUR_DIRECTORY/w_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[8]],file = paste("/YOUR_DIRECTORY/theta_",i,".rds",sep = ""),version = 2)
  
  data_seed<-data_seed+50000
  print(c("varying p ",i))
}

# Fix p, range q
q_regimes<-c(0.01,seq(0.05,0.45,0.05))
data_seed<-2

for(i in q_regimes){
  p<-rep(0.1,3)
  q<-rep(i,3)
  Sim_mix<-simulate_data_e(repres_study_1,p,q,180,3,data_seed)
  z_tr_mix<-c(rep(1,60),rep(2,60),rep(3,60))
  prob_vec_mix<-prob_vec_fun(z_tr_mix,Sim_mix)
  res<-MCMC_mix_err_sbm(3,2,500000,0,21,180,Sim_mix,prob_vec_mix,a_0,b_0,c_0,d_0,e_0,f_0,z_tr_mix,c_study_1,0.01,0.01,0.3,pert,rw_stps)
  
  saveRDS(res[[1]],file = paste("/YOUR_DIRECTORY/tau_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[2]],file = paste("/YOUR_DIRECTORY/z_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[3]],file = paste("/YOUR_DIRECTORY/p_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[4]],file = paste("/YOUR_DIRECTORY/q_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[6]],file = paste("/YOUR_DIRECTORY/c_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[7]],file = paste("/YOUR_DIRECTORY/w_",i,".rds",sep = ""),version = 2)
  saveRDS(res[[8]],file = paste("/YOUR_DIRECTORY/theta_",i,".rds",sep = ""),version = 2)
  
  data_seed<-data_seed+50000
  print(c("varying q ",i))
}

# Visualise error bars for posterior means and credible intervals for p,q Figures 8 and 9 in Section 5.1

for(i in p_regimes){
  p<-readRDS(paste("/YOUR_DIRECTORY/p_",i,".rds",sep = ""))
  if (i==0.01){
    p_mean<-c(mean(p[-c(1:150000),1]),mean(p[-c(1:150000),2]),mean(p[-c(1:150000),3]))
    p_ci1<-p.interval(p[-c(1:150000),1])
    p_ci2<-p.interval(p[-c(1:150000),2])
    p_ci3<-p.interval(p[-c(1:150000),3])
    p_min<-c(p_ci1[1,1],p_ci2[1,1],p_ci3[1,1])
    p_max<-c(p_ci1[1,2],p_ci2[1,2],p_ci3[1,2])
  }else{
    p_mean<-c(p_mean,mean(p[-c(1:150000),1]),mean(p[-c(1:150000),2]),mean(p[-c(1:150000),3]))
    p_ci1<-p.interval(p[-c(1:150000),1])
    p_ci2<-p.interval(p[-c(1:150000),2])
    p_ci3<-p.interval(p[-c(1:150000),3])
    p_min<-c(p_min,p_ci1[1,1],p_ci2[1,1],p_ci3[1,1])
    p_max<-c(p_max,p_ci1[1,2],p_ci2[1,2],p_ci3[1,2])
  }
}
p_tr<-rep(c(0.01,seq(0.05,0.45,0.05)),each=3)
assign(paste("Cluster","c",sep=" "),rep(c("1","2","3"),10))
data=data.frame(p_tr, `Cluster c` , p_mean, p_min, p_max)
p_gg<-ggplot(data, aes(x=p_tr, y=p_mean, colour=`Cluster c`),cex=5) + 
  geom_errorbar(aes(ymin=p_min, ymax=p_max), width=.035, position=position_dodge(0.05))

# color scale that appears also in black-white version
p_gg+theme(axis.text.x = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x=element_text(size=12), axis.text.y = element_text(size=10)) + 
  theme(legend.text=element_text(size=10),legend.title = element_text(size=12),legend.position = "bottom",  panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour="black"))+
  scale_x_continuous(name=TeX("true $p_{c}$"),breaks=seq(0, 0.45, by = 0.05),expand = c(0,0))+scale_y_continuous(name=TeX("posterior for $p_{c}$"),breaks=seq(0, 0.45, by = 0.05),expand = c(0,0))+
  labs(fill = "c") +geom_point(position=position_dodge(0.05),size=.7)+scale_color_manual(values=c("sienna2", "sienna4", "lightblue3"))+
  geom_abline(color="red", linetype="dashed",size=.2)

for(i in q_regimes){
  q<-readRDS(paste("/YOUR_DIRECTORY/q_",i,".rds",sep = ""))
  if (i==0.01){
    q_mean<-c(mean(q[-c(1:150000),1]),mean(q[-c(1:150000),2]),mean(q[-c(1:150000),3]))
    q_ci1<-p.interval(q[-c(1:150000),1])
    q_ci2<-p.interval(q[-c(1:150000),2])
    q_ci3<-p.interval(q[-c(1:150000),3])
    q_min<-c(q_ci1[1,1],q_ci2[1,1],q_ci3[1,1])
    q_max<-c(q_ci1[1,2],q_ci2[1,2],q_ci3[1,2])
  }else{
    q_mean<-c(q_mean,mean(q[-c(1:150000),1]),mean(q[-c(1:150000),2]),mean(q[-c(1:150000),3]))
    q_ci1<-p.interval(q[-c(1:150000),1])
    q_ci2<-p.interval(q[-c(1:150000),2])
    q_ci3<-p.interval(q[-c(1:150000),3])
    q_min<-c(q_min,q_ci1[1,1],q_ci2[1,1],q_ci3[1,1])
    q_max<-c(q_max,q_ci1[1,2],q_ci2[1,2],q_ci3[1,2])
  }
}
q_tr<-rep(c(0.01,seq(0.05,0.45,0.05)),each=3)
assign(paste("Cluster","c",sep=" "),rep(c("1","2","3"),10))
data=data.frame(q_tr, `Cluster c` , q_mean, q_min, q_max)
q_gg<-ggplot(data, aes(x=q_tr, y=q_mean, colour=`Cluster c`),cex=5) + 
  geom_errorbar(aes(ymin=q_min, ymax=q_max), width=.035, position=position_dodge(0.05))
q_gg+theme(axis.text.x = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x=element_text(size=12), axis.text.y = element_text(size=10)) + 
  theme(legend.text=element_text(size=10),legend.title = element_text(size=12),legend.position = "bottom")+scale_x_continuous(name=TeX("true $q_{c}$"),breaks=seq(0, 0.45, by = 0.05))+scale_y_continuous(name=TeX("posterior for $q_{c}$"),breaks=seq(0, 0.45, by = 0.05))+
  labs(fill = "c") +geom_point(position=position_dodge(0.05),size=.7)#+xlim(-.5,0.5)
# color scale that appears also in black-white version
q_gg+theme(axis.text.x = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x=element_text(size=12), axis.text.y = element_text(size=10)) + 
  theme(legend.text=element_text(size=10),legend.title = element_text(size=12),legend.position = "bottom",  panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour="black"))+
  scale_x_continuous(name=TeX("true $q_{c}$"),breaks=seq(0, 0.45, by = 0.05),expand = c(0,0))+scale_y_continuous(name=TeX("posterior for $q_{c}$"),breaks=seq(0, 0.45, by = 0.05),expand = c(0.01,0))+
  labs(fill = "c") +geom_point(position=position_dodge(0.05),size=.7)+scale_color_manual(values=c("sienna2", "sienna4", "lightblue3"))+
  geom_abline(color="red", linetype="dashed",size=.2)