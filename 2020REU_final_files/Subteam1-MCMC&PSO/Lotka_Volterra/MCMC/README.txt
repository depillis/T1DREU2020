MCMC

This folder contains directories for 3 implementations of different MCMC methods on the Lotka-Volterra predator-prey model given the dataset HaresLynxData.mat

The following MCMC methods are implemented
1) Metropolis-Hastings: adapted from Astrostats R code and using MATLAB mhsample function
2) Delayed Rejection Adaptive Metropolis: adaptive from M. J. Laine algae example


PARAMETER ESTIMATION DESCRIPTION:
The goal is to fit 4 parameters to a linear model: a, the intrisic growth rate of the prey population, b, the rate of predation, c, the mortality rate of the predator population in the absence of prey, and d, the growth rate of the predator population per prey, to produce a posterior density function for each of the parameters.
Prior distributions are assumed to be uniform, and we use a Gaussian likelihood function. 


1) astrostats_tutorial
Adapted from code used to fit a linear regression model: astrostats_tutorial.m, which was created based on R code from 'Astrostats: Bayesian parameter estimation and model comparison'  by C.A.L Bailer Jones.


2) DRAM_tutorial **This is the folder in which code and figures for the results in the final write-up can be found**

The tutorial in this folder is geared toward the DRAM parameterization, however, the mcmcstat library does have an option for Metropolis-Hastings. The code here is adapted
From M.J.Laine code found in the folder MJLaine_Algae_example.


3) mhsample_code
This folder contains a tutorial for using the built-in MATLAB function mhsample that performs Metropolis MCMC. The mhsample function is implemented in Reuel_Smith_mhsample.m and based on code from: https://crr.umd.edu/bayesian-parameter-estimation-single-data-set-example-problem-52-matlab


