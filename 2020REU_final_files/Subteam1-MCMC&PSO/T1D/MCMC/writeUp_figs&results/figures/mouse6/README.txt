mouse 6

This folder contains the figures from running the DRAM parameterization routine using
the Mouse 6 data set

Conditions of the parameterization were:
- Mouse 6
- Model: NOD w/wave
- uninformative prior: uniform

Figures and .csv were generated by variables in the workspace 
"workspaces/jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal.mat"

Figures marked with * appear in the final write up
------------------------------------------------------------------------------------------
95_CI_mouse6_only.fig
	mcmcstat library generated plot of glucose predictions with 95% confidence interval;
	averaged data is plotted over; time span is 0 to 350 days
	
* 95_CI_mouse6.fig 
	mcmcstat library generated plot of glucose predictions with 95% confidence interval;
	Mouse 6 and remaining 8 acute mouse raw data is plotted over; time span is 0 to 350 days
	

* final_dram_mouse6_traceplot.fig
	traceplots of 10 parameters (excludes burn-in)
	
jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal_chainstats.txt
	doc with the means, standard deviations, geweke values for the 10 parameters

* jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal_den.fig
	plots of PDFs for the 10 parameters
	
* jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal_meanpred.fig
	plot of the glucose prediction using the mean parameter values post-parameterization
	with Mouse 6 data plotted over

jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal_modpred.fig
	mcmcstat library generated plot of glucose predictions with 95% confidence interval;
	Mouse 6 data is plotted over; time span is 60 to 220 days (the time span of the 
	Mouse 6 set in days)

jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal_modpred_states1.fig
	plots of predictions of the states Resting Macrophages, Activated Macrophages, 
	Apoptotic beta cells, Necrotic beta cells, Healthy beta cells, and glucose

jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal_modpred_states2.fig
	plot of prediction of the states Insulin, Immunogenic DCs, Tolerogenic DCs, Effector T
	cells, Regulatory T cells, and Memory T cells
	
jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal_modpred_visvalfull.fig
	plot of glucose prediction with 9 acute mice and averaged data set plotted over top
	
jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal_modpred_predmodData.csv
	.csv file with the data of the glucose prediction from day 0 to 350

------------------------------------------------------------------------------------------
pngVersions contains .png versions of each of these figures
