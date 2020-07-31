This folder contains the figures comparing results
running 4 algorithms
- Joint UKF		- PSO
- Dual UKF		- DRAM MCMC

------------------------------------------------------------------------------------------
Figures: 

averagethenfit_pso_dram_comp.fig/.png
	figure comparing glucose predictions of PSO/MCMC paramterized using the acute averaged
	data set

final_avgComp.fig/.png
	figure comparing glucose predictions of all 4 algorithms. PSO/MCMC paramterized using
	the acute averaged data set and UKFs ran algorithms multiple times then averaged
	final parameter values
	
Mouse6_ComparisonFigure_AllAlgorithms.fig
	figure comparing glucose predictions of all 4 algs parameterized on Mouse 6 data

NoWave_Glucose_Comparison.fig
	Figures comparing glucose predictions using ODE model options NOD+no wave using 
	each algorithms best parameter values; a biological check figures
	
T1D_Averaging Techniques_PSOUKF_New.fig
	plots comparing the parameterization methods of fitting parameter estimates then 
	averaging parameter values vs. averaging data first then parameterizing
