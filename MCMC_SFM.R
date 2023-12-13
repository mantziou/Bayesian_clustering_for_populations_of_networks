#############################################################################################################################
########### Main code for MCMC algorithm for Sparse Finite Mixture (SFM) extension of our model ###########
############################################ as presented in Section 4.5  ############################################
#############################################################################################################################

library(stats)
library(boot)
library(igraph)
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

MCMC_mix_err_sbm_SFM<-function(C,B,iter,burn_in,n_nodes,N_data,G_data,prob_vec,a_0,b_0,c_0,d_0,e_0,f_0,chi_0,a_e,b_e,z_init,c_init,p_init,q_init,theta_init,e0_init,G_m_init,perturb,rw_steps,rw_steps_e0){ 
  
  ### START ARGUMENTS DESCRIPTION ###
  # C: number of clusters, B: number of blocks, iter: number of iterations of MCMC, burn_in: number of iterations to burn,
  # N_data: number of network observations, n_nodes: number of nodes of network data,
  # G_data: list of vectors/upper triangle matrices for symmetric networks, 
  # prob_vec: list of vector probabilities to produce proposal representative network in each cluster
  # a_0, b_0,c_0,d_0,e_0,f_0 each a vector of C length, denoting hyperparameters for Beta prior, for different clusters for p,q respectively,
  # chi_0: vector of length C for hyperparameters of Dirichlet
  # a_e,b_e: scalars, hyperparameters for Gamma prior of e0
  # z_init: vector of initial cluster membership of networks, c_init: list of vectors, one for each cluster, denoting initial block membership of the nodes for representative of each cluster
  # p_init, q_init: vectors with initial values of p and q for each cluster, theta_init: initial val for theta,
  # e0_init: scalar initialise e0, G_m_init: list of length C for initial representatives of each cluster
  # perturb: vector with elements the perturbations of the representatives'edges,
  # rw_steps: vector with elements step sizes of RW proposals for p and q
  # rw_steps_e0: step sizes of RW proposals for e0
  ### END OF ARGUMENTS DESCRIPTION ###
  
  # chains for posterior samples                                          
  chain_G_m<-lapply(1:C, function(x) matrix(NA, nrow = iter-burn_in, ncol = length(G_data[[1]]))) #list of chains for Gm for different clusters C
  chain_p<-matrix(NA, nrow = iter-burn_in, ncol = C)
  chain_q<-matrix(NA, nrow = iter-burn_in, ncol = C)
  chain_taus<-matrix(NA, nrow = iter-burn_in, ncol = C)
  chain_z<-matrix(NA, nrow = iter-burn_in, ncol = N_data)
  chain_c<-lapply(1:C, function(x) matrix(NA, nrow = iter-burn_in, ncol = n_nodes)) #list of chains for c membership of centroids' nodes, for each cluster C
  chain_w<-lapply(1:C, function(x) matrix(NA, nrow = iter-burn_in, ncol = B)) #list of chains for w prob of membership of centroids' nodes for each cluster C
  chain_theta<-lapply(1:C, function(x) matrix(NA, nrow = iter-burn_in, ncol = (B*(B+1))/2)) #list of chains for theta block probabilities for each cluster C
  chain_C_plus<-matrix(NA, nrow = iter-burn_in, ncol = 1) # chain for posterior draws for C+ keeping the non-empty/active clusters 
  chain_e0<-matrix(NA, nrow = iter-burn_in, ncol = 1) # chain for e0 which is the parameter of sparse Dirichlet
  
  
  #initialization of current values
  G_m_current<-G_m_init
  p_current<-rep(p_init,C) 
  q_current<-rep(q_init,C) 
  z_current<-z_init
  e0_current<-e0_init 
  taus_current<-rep(1,C)
  c_current<-c_init #list with C vectors of the block membership of the representatives' nodes in each of the C clusters
  theta_current<-lapply(1:C, function(x) rep(theta_init,B*(B+1)/2)) #  B(B+1)/2 corresponds to the number of elements in upper triangl and the diagonal
  sam<-1 #counter for samples after burn_in
  w_current<-list()
  for (c in 1:C){#Gibbs step Initialization of w_c for each cluster c using initialized c_current (weights for block membership of nodes)
    h2<-tabulate(c_current[[c]]) #counts number of 1s,2s etc in c_current[[c]]
    if(length(h2)<B){h2<-c(h2,rep(0,B-length(h2)))} #if a block has remained with no allocation, fill with 0s missing due to tabulate fun
    w_current[[c]]<-rdirichlet(1,chi_0[c]+h2) 
  }
  
  
  #create names for the representatives entries and for theta entries, corresponding to the block membership of the nodes, for each cluster
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
  for (k in 1:C){# give names to elements of vectors in list (theta vector for each cluster) where these names remain the same throughout
    names(theta_current[[k]])<-lab_theta # give names to each vector of theta probabilities for each cluster C
  }
  G_m_current<-lab_adj_list(G_m_current,c_current,lab) # call function lab_adj_list that returns a list with vectorized ajdacencies with names for each element, given the current membership for each cluster
  
  for (i in 1:iter){
    
    # Calculate current partition
    N_part<-tabulate(z_current) # Partition of N networks to clusters 1:C (counts number of 1s,2s etc in z_current)
    if(length(N_part)<C){N_part<-c(N_part,rep(0,C-length(N_part)))} #if a cluster has remained with no allocation, fullfil with 0s missing due to tabulate fun

    # Update C+ active clusters of current iteration given current partition N_part
    C_plus<-C-length(which(N_part==0))  
    
    # Update τ1,...,τc weights for clusters using Stick-breaking representation
    # Sample sticks from posterior Beta distr, and obtain weights taus from the stick breaking representation
    v_sticks <- rbeta(C-1, e0_current+N_part[1:length(N_part)-1], e0_current*(C-seq(1,C-1,1))+ rev(cumsum(rev(N_part)))[-1])
    taus_current[1]<-log(v_sticks[1])
    taus_current[2:(C-1)]<-log(v_sticks[2:(C-1)])+cumsum(log(1-v_sticks))[1:(C-2)]
    taus_current[C]<-cumsum(log(1-v_sticks))[(C-1)]
    taus_current<-exp(taus_current-logSumExp(taus_current))
    
    #MH Steps and Gibbs steps
    for (c in 1:C){

      if (N_part[c]!=0){
        v<-rmultinom(1,1,c(0.5,0.3,0.4,0.5,0.3,0.3,rep(c(0.2,0.2,0.3,0.4,0.4,0.2),2))) #to see if update Gm or w or q or q for c class and in which way to update
        
        # MH Step: Update Gm,c
        if (1 %in% v[1:6]){ # Update representative with one of the 6 possible proposals
          ind_cl <- which(z_current==c) # find which graphs belong to c class, gives their index to identify them in list of G_data
          G_m_current[[c]] <- MH_update_rep(ind_cl,G_m_current[[c]],G_data,perturb[which(v==1)],c_current[[c]],lab,
                                            p_current[c],q_current[c],theta_current[[c]],n_nodes,a_0[c],b_0[c],c_0[c],d_0[c])
          
        }
        
        # MH Step: Update p_c
        if(1 %in% v[7:12]){ 
          ind_cl <- which(z_current==c) #find which graphs belong to c class, gives their index to identify them in list of G_data
          p_current[c] <- MH_update_p(p_current[c],rw_steps[which(v==1)-6],G_m_current[[c]],q_current[c],G_data,ind_cl,theta_current[[c]],
                                      c_current[[c]],n_nodes,lab,a_0[[c]],b_0[[c]],c_0[[c]],d_0[[c]])
        }

        # MH Step: Update q_c
        if(1 %in% v[13:18]){ 
          ind_cl <- which(z_current==c) #find which graphs belong to c class, gives their index to identify them in list of G_data
          q_current[c] <- MH_update_q(p_current[c],rw_steps[which(v==1)-12],G_m_current[[c]],q_current[c],G_data,ind_cl,theta_current[[c]],
                                      c_current[[c]],n_nodes,lab,a_0[[c]],b_0[[c]],c_0[[c]],d_0[[c]])
        }
        
        #Gibbs step: Update w_c for each cluster c weights for block membership of nodes
        h2<-tabulate(c_current[[c]]) 
        if(length(h2)<B){h2<-c(h2,rep(0,B-length(h2)))} 
        w_current[[c]]<-rdirichlet(1,chi_0[c]+h2) 
        
        #Gibbs step: Update theta_c (probabil of edge for given block membership of nodes) for each cluster c
        for (b in 1:length(theta_current[[1]])){ #update each theta value for possible pairs of blocks
          ind<-names(theta_current[[c]][b])
          adj_ind<-G_m_current[[c]][which(names(G_m_current[[c]]) %in% ind)]
          n_ind<-length(adj_ind)
          adj_ind<-sum(adj_ind)
          theta_current[[c]][b]<-rbeta(1,adj_ind+e_0[c],f_0[c]+n_ind-adj_ind) 
        }
        
        #Gibbs step: Update c_c (membership of nodes of representative in each cluster c)
        pr<-array() #array of probs for node i to belong to group b
        for(n in 1:n_nodes){
          adj_local<-G_m_current[[c]][c(which(match(lab[,"col"],as.character(n))==1),which(match(lab[,"row"],as.character(n))==1))] #identify entries of adjacency that involve node i
          j_nodes<-setdiff(seq(1,n_nodes,1),n) #define the j nodes with i!=j
          c_local<-c_current[[c]][j_nodes] #find the membership of the j nodes
          for(b in 1:B){
            nam<-ifelse(c_local<b,paste(c_local,b,sep = ""),paste(b,c_local,sep = "")) #paste the membership of the j nodes with b of current iteration
            theta_local<-theta_current[[c]][match(nam,names(theta_current[[c]]))] #create a theta local that corresponds to the entries of the adjacency basis the memberships of the j nodes and b
            pr[b]<-mem_c_prob(w_current[[c]][b],theta_local,adj_local)
          }
          l_2<-logxpy_2(pr[1],pr[2]) #normalize to sum to 1 for 2 Blocks (B=2)
          pr<-c(exp(pr[1]-l_2),exp(pr[2]-l_2))
          c_current[[c]][n]<-which(rmultinom(1,1,pr)!=0,arr.ind=TRUE)[,"row"]
        }
        G_m_current[[c]]<-lab_adj(G_m_current[[c]],c_current[[c]],lab)
      }else{
        # Sample from priors of param
        p_current[c]<-rbetat(1,c(0,0.5),a_0[c],b_0[c]) # Sample from truncated beta in range (0,0.5)
        q_current[c]<-rbetat(1,c(0,0.5),c_0[c],d_0[c])
        theta_current[[c]][1:(B*(B+1)/2)]<-rbeta(B*(B+1)/2,e_0[c],f_0[c])
        w_current[[c]]<-rdirichlet(1,rep(chi_0[c],B)) 
        c_current[[c]]<-which(rmultinom(n_nodes,1,w_current[[c]])!=0,arr.ind=TRUE)[,"row"] # block membership of n_nodes
        G_m_current[[c]]<-lab_adj(G_m_current[[c]],c_current[[c]],lab) # give labels to centroid of c according to sampled labels above c_current[[c]]

        theta_temp<-theta_current[[c]][match(names(G_m_current[[c]]),names(theta_current[[c]]))]
        G_m_current[[c]][1:length(G_m_current[[c]])]<-rbinom(length(G_m_current[[c]]),1,theta_temp)
      }
    }# end of iterations over c
    
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
        p[l,m]<-Membership_prob_e(taus_current[m],p_current[m],q_current[m],G_m_current[[m]],G_data[[l]])
      }
      p_temp<-exp(p[l,]-logxpy_C(p[l,]))    
      p[l, ]<-p_temp
    }
    p<-t(apply(p,1,function(x) {x/sum(x)})) #normalize the probability matrix in order for the rows to sum to 1
    z_current<-t(rMultinom(p,1))

    
    #Keep updates in chain if iteration>burn_in
    if(i>burn_in){
      chain_taus[sam, ]<-taus_current
      chain_C_plus[sam,1]<-C_plus
      chain_e0[sam,1]<-e0_current
      for (k in 1:C){
        chain_G_m[[k]][sam, ]<-G_m_current[[k]] 
        chain_c[[k]][sam, ]<-c_current[[k]]
        chain_w[[k]][sam, ]<-w_current[[k]]
        chain_theta[[k]][sam, ]<-theta_current[[k]]
      }
      chain_p[sam, ]<-p_current
      chain_q[sam, ]<-q_current
      chain_z[sam, ]<-z_current
      sam<-sam+1
    }

  }# end of MCMC iterations i
  return(list('taus'=chain_taus,'z'=chain_z,'p'=chain_p,'q'=chain_q,'G_m'=chain_G_m,'c'=chain_c,'w'=chain_w,'theta'=chain_theta,'C_plus'=chain_C_plus,'e0'=chain_e0)) 
}

