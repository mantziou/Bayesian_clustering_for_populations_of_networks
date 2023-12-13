### Description of code and data files for “Bayesian model-based clustering for populations of network data” 

The path names will need to be specified by the user depending on where they store the files on their computer.

The following files are included in the main branch of this repository:

1. MCMC_mixture_e_sbm_rev.R which contains main code for MCMC algorithm for mixture of measurement error models presented in Section 4.3.
2. Log_post_e_sbm_rev.R which contains functions for calculating log posteriors for Metropolis-Hastings ratio calculated in MCMC, both for clustering with cluster specific network representative and outlier cluster version of the model (Sections 4.3 and 4.4)
3. Membership_prob_function_e_sbm.R which contains function for calculating the probabilities for cluster membership of networks
4. normalizing_log.R which contains code for normalising probabilities
5. label_fun.R which contains functions for assigning names to elements of vectorised adjacency matrix, according to the block membership of the nodes involved in corresponding entry
6. Membership_c_prob_function_e_sbm.R which contains code for calculating the probabilities for block membership of nodes
7. MH_updates.R which contains code for Metropolis-Hastings steps in MCMC for false positive/false negative probabilities and network representatives (both for clustering model and outlier cluster model)
8. MCMC_SFM.R which contains main code for MCMC algorithm for Sparse Finite Mixture (SFM) extension of our model as presented in Section 4.5
9. Log_post_e0.R which contains code for calculating log posterior for Metropolis-Hastings ratio for e_0 parameter of symmetric Dirichlet in SFM extension of our model
10. Trunc_beta.R which contains function for drawing from Truncated Beta distribution
11. MCMC_outlier_SFM.R which contains main code for MCMC algorithm for Sparse Finite Mixture (SFM) extension of our outlier cluster model as presented in Sections 4.4 and 4.5
12. MCMC_outlier_cluster_sbm.R which contains main code for MCMC algorithm for outlier version of our mixture of measurement error models presented in Section 4.4
13. Simulations_Section_5.1.R which contains code for simulating data of Section 5.1, as well as tuning and running MCMC for simulated data, and visualising results as presented in Section 5.1 and Supplement Section 2.1
14. Simulations_Section_5.2.R which contains code for simulating data of Section 5.2, as well as tuning and running MCMC for simulated data, and visualising results as presented in Section 5.2 and Supplement Section 2.2
15. Simulations_Section_5.3.R which contains code for simulations of Section 5.3
16. simul_fun.R which contains functions for simulating data of Section 5
17. vec_to_graph.R which contains code for converting vectorised adjacency of undirected networks into graph and/or adjacency matrix, as well as vectorising adjacency
18. Clustering_Entropy.R which contains code for calculating clustering entropy index for the case of C=3 clusters (for general case of C clusters readily available functions are used from R package “NMF” in corresponding scripts)
19. Clustering_Purity.R which contains code for calculating clustering purity index for the case of C=3 clusters (for general case of C clusters readily available functions are used from R package “NMF” in corresponding scripts)
20. Brain_data_application.R which contains code for tuning and running MCMC for outlier cluster model on brain data, as well as code for visualisations presented in Section 6.2
21. neurodata_edge_lists which contains brain data of Section 6.2 (brain edge lists for 30 individuals, with 10 measurements taken for each individual)

Please note:
- Data and code of Section 6.1 (Tacita data) are not provided herein due to data provider restriction.
- No guarantees can be given about the stability of the code and all code is run at the user’s risk.