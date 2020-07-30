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

-----------------
1) lotka_volterra_mh.m:
Main function to fit model with parameters, perform plotting.

The following functions produce components of the posterior probability from Baye's Theorem:

2) logpostLVmodel.m:
Returns the log10 of the posterior probability (via Bayes Theorem) of the inputted parameters and observed data

3) loglikeLVmodel.m:
Returns the log10 of the likelihood function evaluated at the inputted set of parameters and observed data.

4) logprior.m:
Returns the log10 of the normalized prior computed from the inputted set of parameters.

5) makecovmatrix.m:
Returns a covariance matrix given a vector of the standard deviations and the global (scalar) correlation coefficient

6) metroMCMC.m:
Runs the Metropolis-Hastings MCMC algorithm on an input of starting parameters, sample size, a covariance matrix, burn-in, observed data and number of times to run the algorithm. It outputs an N-by-6 matrix containing the posteriors and parameter values at each step in the Markov chain that was accepted (for more details see code comments). This function implements the previous 4 functions.

7) Lotka_Volterra_Model.m:
Structure of Lotka Volterra system of ODEs to be called in ODE solver. Requires set of y values, as well as k: a vector of the parameters.

8) HaresLynxData.mat:
Data set for system where column 1 is time, column 2 is x, column 3 is y
