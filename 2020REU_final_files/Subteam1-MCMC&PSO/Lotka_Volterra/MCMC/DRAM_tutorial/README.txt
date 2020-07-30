lotkaVolterra_mh_dram README.txt

This folder contains data and functions to parameterize the Lokta-Volterra predator-prey 
model. We use a script adapted from an algae example by M.J. Laine and rely on the 
functions from the mcmcstat library.

lotkaVolterraex_2.m
	Run the parameterization routine from this function.

lotkaVolterrafun.m
	Function to solve the ODE
	
lotkaVolterrass.m
	Calculates a sum of squares for use in the acceptance criteria
	
lotkaVolterrasys.m
	Function that describes the system of equations.
	
LVmse.m
	Calculates mean squared error of the model prediction
	
makeSimData.m
	Creates simulated data by solving the LV ODE and adding 10% Gaussian noise
	
error_est.m
	Calculates standard error of the estimate (another measure of algorithm performance)
	
fitInitialParams.m 
	Using fmincon finds minimum value of parameter values for which system is viable. 
	These values are used as initial guesses for parameter values. Also calculates an
	estimate for the variance of the prior function (model.sigma2).
	
logprior.m
	Function for a user-specified prior function. Default in this implementation is a
	log-uniform prior.

PDFoverlay.m
	Function to plot the posterior distributions of each parameter from both a 		Metropolis-Hastings and DRAM parameterization on top of each other for comparison

compare_algs.m
	Function to plot clustered bar charts of the means of parameters with error bars 	created with 1 standard deviation

***
Note that within this tutorial there is an option to perform Metropolis-Hastings MCMC