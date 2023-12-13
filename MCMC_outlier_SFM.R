#############################################################################################################################
########### Main code for MCMC algorithm for Sparse Finite Mixture (SFM) extension of our outlier cluster model ###########
############################################ as presented in Sections 4.4 and 4.5  ############################################
#############################################################################################################################

library(stats)
library(boot)
library(igraph)
library(sand)
library(coda)
library(nnet)
library(Matrix)
library(MCMCprecision)
library(Hmisc)
library(rlist)
library(plyr)
library(stringr)
library(matrixStats)

# source R scripts containing functions called in MCMC 
setwd("/YOUR_DIRECTORY")
source("Log_post_e_sbm_rev.R")
source("Log_post_e0.R")
source("Membership_prob_function_e_sbm.R")
source("normalizing_log.R")
source("label_fun.R")
source("Trunc_beta.R")
source("Membership_c_prob_function_e_sbm.R")
source("MH_updates.R")

MCMC_mix_err_sbm_sfm<-function(C,B,iter,sub_iter,burn_in,n_nodes,N_data,G_data,prob_vec,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,z_init,c_init,p_init,q_init,theta_init,e0_init,G_m_init,perturb,rw_steps,rw_steps_e0){ 
  
  ### START ARGUMENTS DESCRIPTION ###
  # C: number of clusters, B: number of blocks, iter: number of iterations of MCMC, burn_in: number of iterations to burn,
  # sub_iter: scalar to save subset of iterations, set sub_iter equal to iter-burn_in in case not wanting to save subiterations
  # N_data: number of network observations, n_nodes: number of nodes of network data,
  # G_data: list of vectors/upper triangle matrices for symmetric networks, 
  # prob_vec: list of vector probabilities to produce proposal representative network
  # a_0, b_0,c_0,d_0,e_0,f_0 each a vector of C length, denoting hyperparameters for Beta prior, for different clusters for p,q respectively,
  # chi_0: vector of length C for hyperparameters of Dirichlet
  # a_e,b_e: scalars, hyperparameters for Gamma prior of e0
  # z_init: vector of initial cluster membership of networks, c_init: vector of n_nodes length for initial block membership of the nodes for representative
  # p_init, q_init: vectors with initial values of p and q for each cluster, theta_init: initial val for theta,
  # e0_init: scalar initialise e0, G_m_init: initialise single representative network
  # perturb: vector with elements the perturbations of the representative's edges,
  # rw_steps: vector with elements step sizes of RW proposals for p and q
  # rw_steps_e0: step sizes of RW proposals for e0
  ### END OF ARGUMENTS DESCRIPTION ###
  
  #chains for posterior samples                                          
  chain_G_m<-vector("list",sub_iter) 
  chain_p<-lapply(1:C, function(x) vector("list",sub_iter))
  chain_q<-lapply(1:C, function(x) vector("list",sub_iter))
  chain_taus<-vector("list",sub_iter)
  chain_z<-vector("list",sub_iter)
  chain_c<-vector("list",sub_iter)
  chain_w<-vector("list",sub_iter) 
  chain_theta<-vector("list",sub_iter) 
  chain_C_plus<-vector("list",sub_iter) 
  chain_e0<-vector("list",sub_iter) 
  
  #initialization of current values
  G_m_current<-list(G_m_init) # put vector in list format
  p_current<-rep(p_init,C) 
  q_current<-rep(q_init,C) 
  z_current<-z_init
  c_current<-c_init #vector of n_node length for the block membership of the representative
  theta_current<-rep(theta_init,B*(B+1)/2) # B(B+1)/2 corresponds to the number of elements in upper triangle and the diagonal
  e0_current<-e0_init # scalar
  taus_current<-rep(1,C)
  sam<-1 #counter for samples after burn_in
  h2<-tabulate(c_current) 
  if(length(h2)<B){h2<-c(h2,rep(0,B-length(h2)))} #if a block has remained with no allocation, fill with 0s missing due to tabulate fun
  w_current<-rdirichlet(1,chi_0+h2) 
  
  #create names for the represetative entries and for theta entries
  toy_mat<-matrix(nrow = n_nodes,ncol = n_nodes)
  lab<-which(upper.tri(toy_mat),arr.ind = TRUE) #gives the positions of the upper triang elements in matrix of size n_nodes, i.e. the pairs of nodes(possible edges) for undirected networks
  k<-1
  lab_theta<-array()
  for (i in 1:B){
    for (j in i:B){
      lab_theta[k]<-paste(i,j,sep = "") #labels for theta that are combinations of the blocks
      k<-k+1
    }
  }
  names(theta_current)<-lab_theta # give names to vector of theta probabilities 
  G_m_current[[1]]<-lab_adj(G_m_current[[1]],c_current,lab) 
  
  for (i in 1:iter){

    # Calculate current partition
    N_part<-tabulate(z_current) # Partition of N networks to clusters 1:C (counts number of 1s,2s etc in z_current)
    if(length(N_part)<C){N_part<-c(N_part,rep(0,C-length(N_part)))} #if a cluster has remained with no allocation, fill with 0s missing due to tabulate fun
    
    # Update C+ active clusters of current iteration given current partition N_part
    C_plus<-C-length(which(N_part==0))  
    
    # Update τ1,...,τc weights for clusters using Stick-breaking representation
    # Sample sticks from posterior Beta distr, and obtain weights taus from the stick breaking representation
    v_sticks <- rbeta(C-1, e0_current+N_part[1:length(N_part)-1], e0_current*(C-seq(1,C-1,1))+ rev(cumsum(rev(N_part)))[-1])
    taus_current[1]<-log(v_sticks[1])
    taus_current[2:(C-1)]<-log(v_sticks[2:(C-1)])+cumsum(log(1-v_sticks))[1:(C-2)]
    taus_current[C]<-cumsum(log(1-v_sticks))[(C-1)]
    taus_current<-exp(taus_current-logSumExp(taus_current))
    
    # MH Step: Update Gm NON cluster specific
    v<-rmultinom(1,1,c(0.5,0.3,0.4,0.5,0.3,0.3)) #to see which update of Gm to choose
    if (1 %in% v[1:5]){ 
      #obtain vector of p and q with respect to corresponding network membership to cluster 1 or outlier cluster 2
      p_vec_curr<-p_current[z_current[c(1:N_data)]]
      q_vec_curr<-q_current[z_current[c(1:N_data)]]
      G_m_current[[1]] <- MH_update_rep_out(G_m_current[[1]],G_data,perturb[which(v==1)],c_current,lab,
                                            p_vec_curr,q_vec_curr,theta_current,n_nodes)
    }
    
    if (v[6]==1){
      G_m_new<-rbinom(length(G_data[[1]]),1,prob_vec) # generate edges of proposal with diff probs according to prob_vec
      G_m_new<-lab_adj(G_m_new,c_current,lab) #give membership labels to proposed centroid
      d1<-sum(dbinom(G_m_current[[1]],1,prob_vec,log = TRUE)) #proposal distr density for current state
      d2<-sum(dbinom(G_m_new,1,prob_vec,log = TRUE)) #proposal distr density for proposal vector
      p_vec_curr<-p_current[z_current[c(1:N_data)]]
      q_vec_curr<-q_current[z_current[c(1:N_data)]]
      mhr=exp(Logpost_outlier_clust(G_m_new,p_vec_curr,q_vec_curr,G_data,theta_current,c_current,n_nodes,lab)+d1-
                Logpost_outlier_clust(G_m_current[[1]],p_vec_curr,q_vec_curr,G_data,theta_current,c_current,n_nodes,lab)-d2)
      prob<-min(1, mhr)
      gen<-rbinom(1, 1, prob)
      if (gen==1){
        G_m_current[[1]]<-G_m_new
      }
    }
    
    #MH Step: Update  p_c or q_c 
    for (c in 1:C){
      if (N_part[c]!=0){
        v<-rmultinom(1,1,rep(c(0.2,0.2,0.3,0.4,0.4,0.2),2)) #to see if update p or q for c class and in which way to update
        
        # MH Step: Update p_c
        if(1 %in% v[1:6]){ 
          ind_cl <- which(z_current==c) #find which graphs belong to c class, gives their index to identify them in list of G_data
          p_current[c] <- MH_update_p(p_current[c],rw_steps[which(v==1)],G_m_current[[1]],q_current[c],G_data,ind_cl,theta_current,
                                      c_current,n_nodes,lab,a_0[[c]],b_0[[c]],c_0[[c]],d_0[[c]])
        }
        

        # MH Step: Update q_c
        if(1 %in% v[7:12]){ 
          ind_cl <- which(z_current==c) #find which graphs belong to c class, gives their index to identify them in list of G_data
          q_current[c] <- MH_update_q(p_current[c],rw_steps[which(v==1)-6],G_m_current[[1]],q_current[c],G_data,ind_cl,theta_current,
                                      c_current,n_nodes,lab,a_0[[c]],b_0[[c]],c_0[[c]],d_0[[c]])
        }

      }else{
        # Sample from priors of param
        p_current[c]<-rbetat(1,c(0,0.5),a_0[c],b_0[c]) # Sample from truncated beta in range (0,0.5)
        q_current[c]<-rbetat(1,c(0,0.5),c_0[c],d_0[c])
      }
    }
    
    #Gibbs step: Update w-weights for block membership of nodes
    h2<-tabulate(c_current)
    if(length(h2)<B){h2<-c(h2,rep(0,B-length(h2)))} #if a block has remained with no allocation, fill with 0s missing due to tabulate fun
    w_current<-rdirichlet(1,chi_0+h2) 
    
    #Gibbs step: Update theta (probabil of edge for given block membership of nodes)
    for (b in 1:length(theta_current)){ #update each theta value for possible pairs of blocks
      ind<-names(theta_current[b])
      adj_ind<-G_m_current[[1]][which(names(G_m_current[[1]]) %in% ind)]
      n_ind<-length(adj_ind)
      adj_ind<-sum(adj_ind)
      theta_current[b]<-rbeta(1,adj_ind+e_0,f_0+n_ind-adj_ind) 
    }
    
    #Gibbs step: Update c (membership of nodes of representative)
    pr<-array() #array of probs for node i to belong to group b
    for(n in 1:n_nodes){
      adj_local<-G_m_current[[1]][c(which(match(lab[,"col"],as.character(n))==1),which(match(lab[,"row"],as.character(n))==1))] #identify entries of adjacency that involve node i
      j_nodes<-setdiff(seq(1,n_nodes,1),n) #define the j nodes with i!=j
      c_local<-c_current[j_nodes] #find the membership of the j nodes
      for(b in 1:B){
        nam<-ifelse(c_local<b,paste(c_local,b,sep = ""),paste(b,c_local,sep = "")) #paste the membership of the j nodes with b of current iteration
        theta_local<-theta_current[match(nam,names(theta_current))] #create a theta local that corresponds to the entries of the adjacency basis the memberships of the j nodes and b
        pr[b]<-mem_c_prob(w_current[b],theta_local,adj_local)
      }
      l_2<-logxpy_2(pr[1],pr[2]) #normalize to sum to 1 for 2 Blocks (B=2)
      pr<-c(exp(pr[1]-l_2),exp(pr[2]-l_2))
      c_current[n]<-which(rmultinom(1,1,pr)!=0,arr.ind=TRUE)[,"row"]
    }
    G_m_current[[1]]<-lab_adj(G_m_current[[1]],c_current,lab)
    
    # Update e0 using MH step
    e0_cand<-e0_current+runif(1,min = -rw_steps_e0[1],max = rw_steps_e0[1]) # candidate proposal for e0
    if (e0_cand<0){
      e0_new<-(-e0_cand)
    }else{
      e0_new<-e0_cand
    }
    mhr=exp(Logpost_e0(e0_new,N_data,N_part,a_e,b_e)-Logpost_e0(e0_current,N_data,N_part,a_e,b_e))
    prob<-min(1, mhr)
    gen<-rbinom(1, 1, prob)
    if (gen==1){
      e0_current<-e0_new
    }
    
    #Gibbs step: Update z1,...,zc
    p<-matrix(NA,nrow=N_data,ncol=C) #matrix with each row prob
    for (l in 1:N_data){
      for (m in 1:C){
        p[l,m]<-Membership_prob_e(taus_current[m],p_current[m],q_current[m],G_m_current[[1]],G_data[[l]])
      }
      p_temp<-exp(p[l,]-logxpy_C(p[l,]))    
      p[l, ]<-p_temp
    }
    p<-t(apply(p,1,function(x) {x/sum(x)})) #normalize the probability matrix in order for the rows to sum to 1
    z_current<-t(rMultinom(p,1))
    
    #Keep updates in chain if iteration>burn_in
    if(i>burn_in){
      chain_taus[[sam]]<-taus_current
      chain_C_plus[[sam]]<-C_plus
      chain_e0[[sam]]<-e0_current
      chain_z[[sam]]<-z_current
      chain_G_m[[sam]]<-as(G_m_current[[1]],"sparseVector")
      chain_c[[sam]]<-c_current
      chain_w[[sam]]<-w_current
      chain_theta[[sam]]<-theta_current
      for (k in 1:C){
        chain_p[[k]][[sam]]<-p_current[k]
        chain_q[[k]][[sam]]<-q_current[k]
      }
      sam<-sam+1
    }
    # UNCOMMENT IF WANTED TO SAVE SUBITER (FOR MANY MCMC ITER AND BIG NETWORKS)
    # if((i-burn_in)%%sub_iter==0){ 
    #   for(k in 1:C){
    #     saveRDS(do.call(rbind,chain_p[[k]]),file = paste("/YOUR_DIRECTORY/p_",k,"_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #     saveRDS(do.call(rbind,chain_q[[k]]),file = paste("/YOUR_DIRECTORY/q_",k,"_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #   }
    #   saveRDS(chain_G_m,file = paste("/YOUR_DIRECTORY/represent_outl_sfm_",i-burn_in,".rds",sep=""),version=2) #NO do.call(cBind,) very slow for that big vectors, the corresponding rds will result in a 0.8GB object.size() loaded in R (if not sparseVector it would be 3.7GB)
    #   saveRDS(do.call(rbind,chain_taus),file = paste("/YOUR_DIRECTORY/tau_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #   saveRDS(do.call(rbind,chain_z),file = paste("/YOUR_DIRECTORY/z_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #   saveRDS(do.call(rbind,chain_c),file = paste("/YOUR_DIRECTORY/c_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #   saveRDS(do.call(rbind,chain_w),file = paste("/YOUR_DIRECTORY/w_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #   saveRDS(do.call(rbind,chain_theta),file = paste("/YOUR_DIRECTORY/theta_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #   saveRDS(do.call(rbind,chain_C_plus),file = paste("/YOUR_DIRECTORY/C_plus_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #   saveRDS(do.call(rbind,chain_e0),file = paste("/YOUR_DIRECTORY/e0_outl_sfm_",i-burn_in,".rds",sep=""),version=2)
    #   sam<-1 #!!!!!! make sam=1 to refill already created list in-place (because each entry refilled is list)
    # }
    
  }# end of MCMC iterations i
  return(list('taus'=chain_taus,'z'=chain_z,'p'=chain_p,'q'=chain_q,'G_m'=chain_G_m,'c'=chain_c,'w'=chain_w,'theta'=chain_theta,'C_plus'=chain_C_plus,'e0'=chain_e0)) 
}

