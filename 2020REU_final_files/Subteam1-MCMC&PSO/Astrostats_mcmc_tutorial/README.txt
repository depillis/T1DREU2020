MATLAB code to parameterize a simple linear model with unknown noise using Markov Chain Monte Carlo Methods. Translated from R code in 'Astrostats: Bayesian parameter estimation and model comparison'  by C.A.L Bailer Jones.

The goal of this tutorial is to fit 3 parameters to a linear model: a0, the y-intercept, a1, the slope, and ysig, unknown Gaussian noise, to produce a posterior density function for each of the parameters. The linear model itself is fit to 10 simulated datapoints. Prior distributions are assumed to be uniform, and we use a Gaussian likelihood function. Running astrostats.m will simulate sample data, run the MCMC algorithm, and produce plots of the Markov chains and parameter density functions, pausing between plots. Functions 2-6 are called within Function 1. Functions 8 and 9 are used to sample from the final parameter PDFs produced by astrostats_tutorial.m and plot some proposed solutions to the linear model.

*Note that function sampleprior.m is not used within our tutorial but can be implemented if the user wishes to have a prior function with a gamma distribution

1) astrostats_tutorial.m:
A tutorial to fit parameters to a linear model with unknown noise using Bayesian techniques a Metropolis-Hastings MCMC algorithm. We simulate data and compute the parameter PDFs using the following internally implemented functions:\

The following functions produce components of the posterior probability from Baye's Theorem

2) logpostlinearmodel.m:
Returns the log10 of the posterior probability (via Bayes Theorem) of the inputted parameters and observed data

3) loglikelinear.m:
Returns the log10 of the likelihood function evaluated at the inputted set of parameters and observed data.

4) logprior.m:
Returns the log10 of the normalized prior computed from the inputted set of parameters.

5) makecovmatrix.m:
Returns a covariance matrix given a vector of the standard deviations and the global (scalar) correlation coefficient

6) metroMCMC.m:
This function runs the Metropolis-Hastings MCMC algorithm on an input of starting parameters, sample size, a covariance matrix, burn-in, observed data and number of times to run the algorithm. It outputs an N-by-5 matrix containing the posteriors and parameter values at each step in the Markov chain that was accepted (for more details see code comments). This function implements the previous 4 functions.

7) sampleprior.m: 
Returns a specified number of samples from a proper gamma prior. (This function is used in a later part of the tutorial, we do not use it with our MCMC algorithm)

The following functions are used to produce sample solutions to the linear model from the parameter distribution computed in functions 1-6:

8) randVals.m:
Produces an N-by-N matrix of proportionately sampled values from PDF estimated using multivariate kernel density estimation (mvksdensity); used to form potential fittings

9) runN.m:
"Runs" the astrostats tutorial N times by selecting and plotting N possible fittings of the data created through sampling the posterior PDFs

10) Subfolder figures_and_data:
Contains visuals and data produced by running astrostats.m (10.1-10.5), runN.m (10.6)

10.1) params.fig:
MATLAB figure file plotting resultant Markov Chains for parameters

10.2) density_plots_true.fig:
MATLAB figure file plotting posterior PDFs of parameters with true parameters marked

10.3) postSamp.csv:
Final value of variable postSamp from astrostats_tutorial.m, giving the log of the posterior probability in columns 1,2 and the corresponding candidate parameters: a0, a1, ysig in columns 3-5

10.4) thetaMap.csv:
Final value of variable thetaMAP from astrostats_tutorial.m, giving the Maximum A Posteriori set of parameters (i.e. parameters producing maximum posterior likelihood) found with MCMC

10.5) thetaMean.csv:
Final value of variable thetaMean from astrostats_tutorial.m, giving the mean-valued set of parameters (i.e. average of all accepted candidate parameters) found with MCMC

10.6) sample_lines_from_paramDen.fig:
MATLAB figure file overplotting 10 possible linear fittings sampled according to the calculated posterior distributions against simulated data