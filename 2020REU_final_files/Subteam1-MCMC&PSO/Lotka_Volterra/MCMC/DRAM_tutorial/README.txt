lotkaVolterra_mh_dram README.txt

This folder contains data and functions to parameterize the Lokta-Volterra predator-prey 
model. We use a script adapted from an algae example by M.J. Laine and rely on the 
functions from the mcmcstat library.

***Note that code from this folder was used to generate results and figures in the final REU write up.***
-------------------------------------------------------------------------------------------------
lotkaVolterraex_2.m
 Authors:      Edited by C. Catlett & M. Watanabe
               Original code from https://mjlaine.github.io/mcmcstat/ex/algaeex.html

 Date:         June 2020

 Descr:        Script to estimate parameters of the Lotka-Volterra system.
               This script can be used to perform the Metropolis-Hastings 
               and DRAM MCMC methods. Specify 'mh' or 'dram' in lines 87 
               and 99.
               Uses functions lotkaVolterrasys.m, lotkaVolterrafun.m, 
               lotkaVolterrass.m, and mcmcstat library

 Directions:   This script is meant to be run in full, ~15 min. User 
               specifications exist for DRAM vs Metropolis-Hastings and 
               naming the workspace. 


	
lotkaVolterraex_2.m CALLS functions:

	lotkaVolterrafun.m
		Function to solve the ODE
	
	lotkaVolterrass.m
		Calculates a sum of squares for use in the acceptance criteria
	
	lotkaVolterrasys.m
		Function that describes the system of equations.
	
	LVmse.m
		Calculates mean squared error of the model prediction
	
	[error_est.m
		Calculates standard error of the estimate (another measure of algorithm 
		performance)] * not used

	
	fitInitialParams.m 
		Using fmincon finds minimum value of parameter values for which system is 
		viable. 
		These values are used as initial guesses for parameter values. Also 
		calculates an estimate for the variance of the prior function (model.sigma2).
	
	logprior.m
		Function for a user-specified prior function. Default in this implementation is a
		log-uniform prior.
-------------------------------------------------------------------------------------------------
DATA functions & files
makeSimData.m
	Creates simulated data by solving the LV ODE and adding 10% Gaussian noise
	Creates the file makesimData.mat

HaresLynxData.mat
	Raw hare and lynx population data from The Hudson Bay Company study (reference in
	References/MCMC directory)
-------------------------------------------------------------------------------------------------
COMPARISON functions: these scripts are used to compare the results of a Metropolis-Hastings parameterization with a DRAM parameterization. Each is a self-contained script and must be run individually to produce desired figures.

compare_algs.m
	Function to plot clustered bar charts of the means of parameters with error bars created with 	1 standard deviation

PDFoverlay.m
	Script to plot the individual parameter PDFs from the DRAM and MH parameterizations. Note 	that this script explicitly uses the data and variables from the workspaces "final_mh.mat" 	and "final_dram.mat". If comparison of subsequent parameterizations are desired, they must be
	specified.
-------------------------------------------------------------------------------------------------
SUB-DIRECTORIES
finalResults&Workspaces
	contains the workspaces and resulting figures that are found in the final REU write up

figures&results
	contains figures and workspaces from past runs of the MCMC parameterization implementation

mcmcstat
	MATLAB library that our parameterization relies on

