
#log-likelihood of partition
Loglikel_partition<-function(e0,N_data,N_partition){ # e0:scalar, N_data: # of network observations, N_partition: number of network data in each cluster
  c_list<-which(N_partition!=0)
  term_1<-0
  for (c in c_list){
    term_1<-term_1+sum(log(e0+seq(0,N_partition[c]-1,1)))
  }
  term_2<-sum(log((length(N_partition)*e0)+seq(0,N_data-1,1)))
  return(term_1-term_2)
}

#log-prior for e0
Logprior_e0<-function(e0,a_e,b_e){
  # Gamma prior for e0
  e0_prior<-((a_e-1)*log(e0))-(b_e*e0)
  return(e0_prior)
}

#Log-posterior for e0
Logpost_e0<-function(e0,N_data,N_partition,a_e,b_e){ 
  return(Loglikel_partition(e0,N_data,N_partition)+Logprior_e0(e0,a_e,b_e))
}
