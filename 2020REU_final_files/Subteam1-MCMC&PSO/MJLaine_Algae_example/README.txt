MJLaine_Algae_example

This folder contains code for an implementation of delayed rejection adaptive metropolis MCMC (DRAM) by M. J. Laine (https://mjlaine.github.io/mcmcstat/ex/algaeex.html)

This script was used as a practice tool to prepare to parameterize the T1D ODE model using the DRAM method and the mcmcstat library

Parameterization routine is run from
	algaeex.m

and calls the functions
	algaefun.m
	algaess.m
	algaesys.m
	 

------------------------
Calibrate_DRAMfittin.m: This function will calibrate your model to the given data via
parameter estimation using the DRAM algorithm. You will need to have the
mcmcstat folder added to your current path for this to run.


mcmcstat
Contains Laine's code for DRAM and accompanying functions