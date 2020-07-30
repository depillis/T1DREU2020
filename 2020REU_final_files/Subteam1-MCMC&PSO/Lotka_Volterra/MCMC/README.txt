This folder contains folders for implementation of different MCMC methods on the Lotka-Volterra predator-prey model given the dataset HaresLynxData.mat

The following MCMC methods are implemented
1) Metropolis-Hastings: adapted from Astrostats R code
2) Delayed Rejection Adaptive Metropolis: adaptive from M. J. Laine algae example

***********************************
1) Metropolis-Hastings
***********************************
1a) astrostats_tutorial
Adapted from code used to fit a linear regression model: astrostats_tutorial.m, which was created based on R code from 'Astrostats: Bayesian parameter estimation and model comparison'  by C.A.L Bailer Jones.

The goal is to fit 4 parameters to a linear model: a, the intrisic growth rate of the prey population, b, the rate of predation, c, the mortality rate of the predator population in the absence of prey, and d, the growth rate of the predator population per prey, to produce a posterior density function for each of the parameters.
Prior distributions are assumed to be uniform, and we use a Gaussian likelihood function. Running lotka_volterra_mh.m will run the MCMC algorithm and produce plots of the Markov chains and parameter density functions, pausing between plots. Functions 2-6 are called within Function 1.

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

10) Subfolder figures_and_data:
Contains visuals and data produced by running lotka_volterra_mh.m

10.1) LV_MH_mcmc.fig:
MATLAB figure file plotting resultant Markov Chains for parameters

10.2) LV_MH_density.fig:
MATLAB figure file plotting posterior PDFs of parameters

10.3) postSamp.csv:
Final value of variable postSamp from lotka_volterra_mh.m, giving the log of the posterior probability in columns 1,2 and the corresponding candidate parameters: a, b, c, d in columns 3-6

10.4) thetaMap.csv:
Final value of variable thetaMAP from astrostats_tutorial.m, giving the Maximum A Posteriori set of parameters (i.e. parameters producing maximum posterior likelihood) found with MCMC

10.5) thetaMean.csv:
Final value of variable thetaMean from astrostats_tutorial.m, giving the mean-valued set of parameters (i.e. average of all accepted candidate parameters) found with MCMC

1b) DRAM_tutorial
The tutorial in this folder is geared toward the DRAM parameterization, however, the mcmcstat library does have an option for Metropolis-Hastings

1c) mhsample_code
This folder contains a tutorial for using the built-in MATLAB function mhsample that performs Metropolis MCMC.


***********************************
2) DRAM
***********************************
This folder contains data and functions to parameterize the Lokta-Volterra predator-prey 
model. We use a script adapted from an algae example by M.J. Laine and rely on the 
functions from the mcmcstat library.

lotkaVolterraex_2.m
	Run the parameterization routine from this folder.

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