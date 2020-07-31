validation&comparison_funcs
This folder contains functions and figures for visually checking parameterization results
and for comparing results from different runs of the parameterization routine.
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

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
	
4) Figures

allState_avg.fig/.png
	 plots of 12 T1D ODE states predictions using parameter means from averaged data 
	 parameterization

allState_comp_3_6_11_avg.fig/.png
	plots of 12 T1D ODE states predictions using parameter means from Mice 3, 6, 11, and 
	averaged data parameterizations

allState_mouse6.fig/.png
	plots of 12 T1D ODE states predictions using parameter means from Mouse 6 data
	parameterization
	
mouse6_importantStates.fig/mcmc_mouse6_importantStates.png
	plots of 4 T1D ODE states predictions using parameter means from Mouse 6 data
	parameterization; states are Apoptotic beta cells, insulin, effector T cells,
	and Regulatory T cells
	
PDFs_3_6_11_avg.fig/.png
	overlaid plots of each parameter's PDF from the parameterizations using Mice 3, 6, 11,
	and averaged data sets

PDFs_6_avg.fig/.png
	overlaid plots of each parameter's PDF from the parameterizations using Mouse 6 and
	averaged data sets