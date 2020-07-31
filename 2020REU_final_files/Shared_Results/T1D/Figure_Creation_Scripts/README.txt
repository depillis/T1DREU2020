Figure_Creation_Scripts

This folder contains the functions, files, and results for comparing the results of 
running 4 algorithms
- Joint UKF		- PSO
- Dual UKF		- DRAM MCMC

------------------------------------------------------------------------------------------
Functions:

error_mse.m
	function to compute the mean squared error values for glucose predictions generated
	from each of the 4 algorithms
	
plot_avgthenfit_Comp.m
	Function to plot the glucose predictions from PSO, MCMC, UKFs parameterizations.
	PSO/MCMC run on averaged data. UKFs perform parameterization multiple times and 
	average final parameter values
	USES averagethenfit_comparison_data

Mouse6_ComparisonFigure_Creation.m
	Function to create the plot comparing glucose predictions parameterzates on Mouse 6
	data from each of the 4 algorithms.
	USES
	- jul10_mouse_run1(noIC)_acute_NOD_waveOn_lietal_errorestData.csv
	- jul10_mouse_run1(noIC)_acute_NOD_waveOn_lietal_predmodData.csv
	- mouse6_DualUKFSol.csv
	- mouse6_joint.csv
	- mouse6PSOfit.csv
	- mouse6PSOmse.csv

No_Wave_Glucose_Comparison 
	Script to produce figure of the glucose simulations under each algorithm's final parameters when wave turned off.
	
UKF_Avg_Data_Comparison_Figure 
	Script to compare performance of averaging techniques for Joint, Dual, and PSO and plot figure

	

Folders:
averagethenfit_comparison_data contains the individual raw data used in the function 
plot_avgComp.m
	- avg_lietal_2.csv: raw averaged data
	- final_mse_avgComp.csv: final MSE values for each algorithm
	- joint_dual_dram_pso_dram_avgPredictions.xlsx/.csv: glucose prediction data from 
	  each alg
	- joint_dual_pso_dram_avgMSEdata.csv: data used to compute MSE values
