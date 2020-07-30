mhsample_code

This folder contains files needed to run the MATLAB mhsample function 
for the Lotka-Volterra predator-prey model.


The mhsample function is implemented in Reuel_Smith_mhsample.m and based on code
from:
https://crr.umd.edu/bayesian-parameter-estimation-single-data-set-example-problem-52-matlab

LL_mhsample.m
Computes the loglikelihoood of the data given an inputed theta value
Returns the scaled sum of the lognormal pdf

PSOPC_LVModel.m | PSOPC.m | least_squares.m
Re-indexes the data and computes a better initial parameter guess via least squares
code from
https://github.com/shtyllab/2020-HMC-REU-Codes/tree/master/GA_Parameter_Identification
