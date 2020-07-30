Subteam1: MCMC & PSO
------------------------------------------------------------------------------------------
	Student Collaborators: Christina Catlett & Maya Watanabe
	HMC REU 2020	
	Advisors: Profs. de Pillis (HMC), Shtylla (Pomona), Edholm (Scripps), and An Do (CGU)
	Funded by the NSF


This subteam worked on two parameter estimation techniques: Markov chain Monte Carlo 
methods and Particle Swarm Optimization.
------------------------------------------------------------------------------------------
*************************
Astrostats_mcmc_tutorial
*************************
We began the summer by researching the MCMC mathematics and theory, using code found in 
the Astrostats_mcmc_tutorial folder as a training exercise. We translated R code from 
'Astrostats: Bayesian parameter estimation and model comparison' by C.A.L Bailer Jones 
into MATLAB. The tutorial estimated parameters for a simple linear model.

***************
Lotka_Volterra
***************
This folder contains codes for MCMC and PSO implementations of parameter estimation of
the Lokta-Volterra ODE model. 
The MCMC folder contains 3 different techniques
	1) We tuned the astrostats tutorial to use on the L-V system
	2) Using the mcmcstat library, we were able to implement a DRAM parameterization
	3) We parameterized using the built-in MATLAB function "mhsample"
	
****************
MCMC_algorithms
****************
This folder contains subfolders of code (written by various authors) for other MCMC
techniques
	1) Adaptive MCMC
	2) Adaptive Tempered Parallel MCMC
	3) Another example for parameterizing a biological ODE
	
**********************
MJLaine_Algae_example
**********************
This folder contains codes written by M. J. Laine (author of the mcmcstat library) that
uses the mcmcstat library to parameterize a biological algae ODE model. This is the
tutorial which we used to help us learn how to implement the DRAM MCMC method.

****
T1D
****
This folder contains codes for MCMC and PSO implementations of parameter estimation of
the Type 1 Diabetes ODE model. It also contains data folders which hold mouse glucose .csv
files and a folder with methods and functions to compare UKF, PSO, and MCMC algorithm
results (mainly visual)