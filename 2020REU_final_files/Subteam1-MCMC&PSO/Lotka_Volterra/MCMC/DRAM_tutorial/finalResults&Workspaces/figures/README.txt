figures

1) comparison
	- dram_mh_meanpredComp.fig: predictions of Hare/Lynx populations using mean parameter
		values from the DRAM and Metropolis parameterizations
		
		Plotted using the data in dram_mcmc_predictedData_LV.csv and mh_mcmc_predictedData_LV.csv
		
	- mh_dram_paramComp1.fig: Bar charts with error bars comparing means of alpha and gamma
		from Met and DRAM parameterizations

	- mh_dram_paramComp2.fig: Bar charts with error bars comparing means of beta and delta
		from Met and DRAM parameterizations
		
	- overlaid_param_pdfs.fig/overlaid_params_pdfs.png: Plots of PDFs of each parameter 
		from both Metropolis and DRAM parameterizations
		
		
2) final_dram: final figures from the DRAM parameterization
	- final_dram_burninchain.fig: trace plot of of the burn-in chain for each parameter
	
	- final_dram_chain.fig: trace plot of the sampling chain for each parameter
	
	- final_dram_chainstats.txt: contains means, stds, and Geweke convergence values for 
		each parameter
		
	- final_dram_den.fig: PDFs of each parameter
	
	- final_dram_meanpred.fig: plots of hare and lynx population predictions using the 
		mean parameter values found in the chainstats file
		
	- final_dram_modpred.fig: plots of the hare and lynx population predictions with a 
		95% confidence interval.
		
	- final_dram_samples.fig/final_dram_samples.png: plots of correlations between 
		parameter values sampled during parameterization

3) final_mh: final figures from the Metropolis parameterization

	- final_mh_burninchain.fig: trace plot of of the burn-in chain for each parameter
	
	- final_mh_chain.fig: trace plot of the sampling chain for each parameter
	
	- final_mh_chainstats.txt: contains means, stds, and Geweke convergence values for 
		each parameter
		
	- final_mh_den.fig: PDFs of each parameter
		
	- final_mh_modpred.fig: plots of the hare and lynx population predictions with a 
		95% confidence interval.
		
	- final_mh_samples.fig/final_dram_samples.png: plots of correlations between 
		parameter values sampled during parameterization
4) mse_files
	- dram_mcmc_MSE.csv: contains the MSE score of the DRAM parameterization. MSE 
		tells us how close our predictions captures the data.
	
	- mh_mcmc_MSE.csv: contains the MSE score the Metropolis parameterization. MSE 
		tells us how close our predictions captures the data.

