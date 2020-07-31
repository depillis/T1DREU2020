comparison
This folder contains the figures comparing the DRAM and Metropolis parameterizations of
the Lokta-Volterra system

Conditions of the parameterization were:
- 'dram' & 'mh' MCMC methods
- uninformative prior: uniform


Figures marked with * appear in the final write up
------------------------------------------------------------------------------------------
* dram_mh_meanpredComp.fig
	predictions of Hare/Lynx populations using mean parameter values from the DRAM and 
	Metropolis parameterizations
		
	- Plotted using the data in dram_mcmc_predictedData_LV.csv and 
	  mh_mcmc_predictedData_LV.csv
		
* mh_dram_paramComp1.fig
	Bar charts with error bars comparing means of alpha and gamma from Met and DRAM 
	parameterizations

* mh_dram_paramComp2.fig
	Bar charts with error bars comparing means of beta and delta from Met and DRAM 
	parameterizations
		
* overlaid_param_pdfs.fig/overlaid_params_pdfs.png
	Plots of PDFs of each parameter from both Metropolis and DRAM parameterizations
		
dram_mcmc_predictedData_LV.csv
	.csv containing the glucose prediction data from running the DRAM parameterization

mh_mcmc_predictedData_LV.csv
	.csv containing the glucose prediction data from running the Metropolis 
	parameterization