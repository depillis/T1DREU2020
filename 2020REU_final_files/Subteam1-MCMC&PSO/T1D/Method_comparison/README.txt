compare_avg_T1Dresults

This folder contains the functions, files, and results for comparing the results of 
running 4 algorithms
- Joint UKF		- PSO
- Dual UKF		- DRAM MCMC
to parameterize the T1D ODE model. Here we compare the results that were run on Li et al
averaged data.

Files marked with *
------------------------------------------------------------------------------------------
Figures: 

final_avgComp.fig/.png
	figure comparing glucose predictions of all 4 algorithms. PSO/MCMC paramterized using
	the acute averaged data set and UKFs ran algorithms multiple times then averaged
	final parameter values
	
averagethenfit_pso_dram_comp.fig/.png
	figure comparing glucose predictions of PSO/MCMC paramterized using the acute averaged
	data set


Functions:

error_mse.m
	function to compute the mean squared error values for glucose predictions generated
	from each of the 4 algorithms
	
plot_avgComp.m
	Function to plot the glucose predictions from PSO, MCMC, UKFs parameterizations.
	PSO/MCMC run on averaged data. UKFs perform parameterization multiple times and 
	average final parameter values

Folders:
comparison_data contains the individual raw data used in the function plot_avgComp.m
	- avg_lietal_2.csv: raw averaged data
	- final_mse_avgComp.csv: final MSE values for each algorithm
	- joint_dual_dram_pso_dram_avgPredictions.xlsx/.csv: glucose prediction data from 
	  each alg
	- joint_dual_pso_dram_avgMSEdata.csv: data used to compute MSE values
