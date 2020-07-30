MCMC
==========================================================================================
		*** RUN PARAMETER ESTIMATION FROM T1D_execute.m ***
==========================================================================================

Type 1 Diabetes Model - Delayed Rejection Adaptive Metropolis MCMC
	Author:       Edited to fit T1D by M. Watanabe & C. Catlett
	Date:          June 2020
	Desc:         Use delayed rejection adaptive Metropolis (DRAM) MCMC 
                      methods to parameterize a 12 equation type 1 diabetes ODE system. 
                      This code produces posterior distributions for chosen parameters
                      and the predicted glucose model. 
                      This code utilizes the mcmcstat library (M.J. Laine).

	Sources: 
       - Algorithm adapted from code by Marko J. Laine
       (https://mjlaine.github.io/mcmcstat/ex/algaeex.html). 

       - Glucose data from Mathews et al (2015) and Li et al (2016).

       - T1D ODE model from Shtylla et al (2019).

  Note: This script uses the mcmcstat library, but the main functions,
  mcmcrun and mcmcpred, have been altered to accomodate user preferences.
  These functions are now mcmcrun_custom (lines 431, 440) and
  mcmcpred_custom (line 480). The original functions still exist in the
  mcmcstat folder.

************
helperFuncs
************
T1Dfun_r.m 
	Uses ode15s to evaluate the T1D ODE system (only for the time frame of the data). Uses
	sickEvents.m.
	
T1Dode_r.m
	Models T1D system
	
T1Dss_r.m
	Calculates the sum of squares between the original data and the evaluated ODE for use 
	in the acceptance criteria of the MCMC algorithm
		
T1Dfun_meanpred.m
	Uses ode15s to evaluate T1D ODE system for the time frame 0:1:350 using the post-DRAM
	mean parameter values. This is to visualize the overall behavior of the glucose
	prediction.
------------------------------------------------------------------------------------------		
2) Data functions

eFastparambounds.m
	Finds +/-user-specified percentage above and below the eFast baseline values (An Do) 
	for parameter range.

initParams.m 
	Uses fmincon to find best/minimal initial parameter values (not used here)
	
makeICbounds.m
	Finds +/-user-specified percentage above and below the non-zero initial conditions 
	for parameter range.
	
makepar_musig.m
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

****************************
validation&comparison_funcs	
****************************
1) Individual validation

vis_validation_lietal.m
	Plots the predicted glucose model (using post-DRAM mean parameter values) with the 
	raw Li et al data for visual comparison
------------------------------------------------------------------------------------------
2) Comparison of Mice 3, 6, 11, averaged parameterizations

comp_3_6_11_avg.m
	Visual comparison of the glucose predictions from DRAM
        parameterizations using the Mice 3, 6, 11, and averaged 
        data sets

compare_priorresults.m
	Function to compare the predictions of DRAM 
        parameterizations with and without an informative prior

comparing_avg_mouse6_DRAMresults.m
	Function to compare means of parameters from results of the DRAM algorithm run on both
	averaged and mouse 6 data.

PDFoverlay.m
	Script to plot overlaid parameter PDFs from 
        parameterizations using from Mice 3, 6, 11 and 
        averaged data
------------------------------------------------------------------------------------------
3) Biological Checks

StateOverlay.m
	Script to plot T1D state estimates using parameter values 
	from Mice 3, 6, 11 and averaged data parameterizations in
	order to determine biological feasibility of parameter
	values.

*********************
writeUp_figs&results
*********************
Contains the figures and workspaces that were used to create the figures in the Summer 2020 final report

********************
past_run_workspaces
********************
Contains all workspaces from various runs of the DRAM parameterization tutorial

****************
Figures&Results
****************
Contains all figures and results from various runs of the DRAM parameterization tutorial

***************
T1D_data_files
***************
A local copy of the mouse glucose data

*********
mcmcstat
*********
MATLAB library (M.J. Laine) containing functions to run DRAM parameterization.
NOTE: This script uses the mcmcstat library, but the main functions,
  mcmcrun and mcmcpred, have been altered to accomodate user preferences.
  These functions are now mcmcrun_custom (lines 431, 440) and
  mcmcpred_custom (line 480). The original functions still exist in the
  mcmcstat folder.





	