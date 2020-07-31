helperFuncs
This folder contains functions called directly by the parameterization routine script
T1D_execute.m
------------------------------------------------------------------------------------------		
------------------------------------------------------------------------------------------		
1) MCMC & ODE Model functions

T1Dfun_r.m 
	Uses ode15s to evaluate the T1D ODE system (only for the time frame of the data). Uses
	sickEvents.m.
	
T1Dode_r.m
	Mathematical ODE model of T1D system
	
T1Dss_r.m
	Calculates the sum of squares between the original data and the evaluated ODE for use 
	in the acceptance criteria of the MCMC algorithm
		
T1Dfun_meanpred.m
	Uses ode15s to evaluate T1D ODE system for the time frame 0:1:350 using the post-DRAM
	mean parameter values. This is to visualize the overall behavior of the glucose
	prediction.

logprior.m
	User-specified function for an informative prior function. Currently the prior is 
	implemented as a LOG-NORMAL function
------------------------------------------------------------------------------------------		
2) Data functions

getData.m
	creates data set to be used in parameterization routine: choose from Li et al., Mathews
	et al., or simulated data.

eFastparambounds.m
	Finds +/-user-specified percentage above and below the eFast baseline values (An Do) 
	for parameter range. Uses the data from eFAST.mat

initParams.m 
	Uses fmincon to find best/minimal initial parameter values (not used here)
	
makeICbounds.m
	Finds +/-user-specified percentage above and below the non-zero initial conditions 
	for parameter range.
	
makepar_musig2.m
	Function to extract parameter means and standard deviations (from UKF results).
	For use in user-specified prior function.
	
makeparambounds.m
	Finds +/-user-specified percentage above and below the UKF baseline values (D. Shenker,
	R. Wander) for parameter range.
	
makeSimData.m
	Creates simulated data using the evaluated ODE and adding 5% Gaussian noise.	
------------------------------------------------------------------------------------------		
3) Error functions

error_est.m
	Calculates the standard error of the estimate.
	
error_mse.m
	Calculates the mean squared error.
------------------------------------------------------------------------------------------		
4) Misc functions

filename.m
	creates a unique string file name for saving post-routine workspaces
------------------------------------------------------------------------------------------		
5) Files

eFAST.mat
	contains baseline and upper/lower bounds of 43 parameters (An Do) * not used
	
Parameter_info.xlxs
	contains parameter baseline, upper/lower bounds (eFAST) values, and "critical parameters"
	as determined by the UKF algorithm tuning

Parameter_Distributions_updated.csv
	contains means and standard deviations from fitting normal distributions to 
	post-UKF parameterizations of the T1D model for use in creating an informative 
	prior function
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
