metropolis_hastings

This folder contains code for implementing the Metropolis-Hastings MCMC algorithm
for the Lotka-Volterra ODE model.

This code was adapted from the R code by C.A.L. Bailer-Jones ('Astrostats: Bayesian 
parameter estimation and model comparison' translated to MATLAB Maya Watanabe &
Christina Catlett. Original translated Astrostats tutorial to be found in the astrostats folder.
==========================================================================
HaresLynxData.mat
Predator-prey data (from 1845-1935)

loglikeLVmodel.m
Function to compute the log-likelihood of the model

logpostLVmodel.m
Function to compute the log-posterior of the model

logprior.m
Function to compute the (uniform) log-priors of the model

lotka_volterra-MH.m
Code to implement MH algorithm
- load data
- initialize parameters, burn-in
- compute sensitivity and covariance matrices
- run metroMCMC
- plot results of chain
- compute density functions
- plot density functions

Lotka_Volterra_Model.m
System of equations functions

(makecovmatrix.m: unused in this implementation)

metroMCMC.m
Function to implement the Metropolis-Hastings MCMC algorithm

modelBasic.m
Alternative function to compute system of equations for model

senseq.m
Function to create the sensitivity matrix (based on the data) needed to compute the 
covariance matrix
