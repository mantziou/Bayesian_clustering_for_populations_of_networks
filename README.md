# Bayesian_clustering_for_populations_of_networks

This repository contains code related to the paper Mantziou, A., Lunagomez, S. and Mitra, R., 2021. Bayesian model-based clustering for populations of network data. Ann. Appl. Stat. 18(1): 266-302 (March 2024) [(link to paper)](https://doi.org/10.1214/23-AOAS1789).

This repository contains the data and code to reproduce the simulations presented in  Mantziou, A., Lunagomez, S. and Mitra, R., 2021.

Clustering of the network observations as well as inferences for the model parameters are made using an MCMC scheme using the following scripts:

1. *MCMC_mixture_e_sbm_rev.R* is the main script for performing inferences for finite mixture model, with predefined number of clusters, presented in Section 4.3 of the paper. 
2. *MCMC_outlier_cluster_sbm.R* which contains main code for MCMC algorithm for performing inferernces for outlier version of our mixture of measurement error models, presented in Section 4.4 of the paper
3. *MCMC_SFM.R* contains main code for MCMC algorithm for Sparse Finite Mixture (SFM) extension of our model which allows inferences for clustering of networks with non-predefined number of clusters, as presented in Section 4.5 of the paper. 

Please note that each of the above R scripts are sourcing other R scripts containing functions for the execution of the MCMC. 

For more details about source R scripts, please see [README.md](docs/README.md) file in *docs* .

In addition, in [README.md](docs/README.md) file in *docs*, description of the rest of the scripts for reproducing the simulation experiments of the paper and the neuroscience data application are given.

Neuroscience data are given in [neurodata_edge_lists](neurodata_edge_lists).
