Averaged_acute
This folder contains the .csv files for aligning and averaging the acute Li et al data sets. It also contains several figures of the averaged data.

------------------------------------------------------------------------------------------
FILES

acute_pivotPoints.fig/acute_pivotPoints.png: plot of each of the 9 acute mice with "pivot 	points" (chosen point of alignment) designated

avg_lietal_2.csv: raw averaged acute data with time in weeks

avgData.fig/avgData.png: plot of the averaged data (from avg_lieta_2.csv) time converted 	into days

consolidated_Lietal.xlsx: contains sheets that show the process of shifting and averaging 	the acute mouse data
	- absTime_alignment: aligning 11 mice into absolute time
	- pivot_alignment: determining pivot point, aligning according to pivot points, a	averaging glucose and time data

lietal_meandat.m: our first attempt to create a function that aligns and averages the Li 	et al data at the onset of diabetes point. However, for our tutorial we have done 	this alignment by hand.

Lietal_w_avgData.fig: plot of 11 Li mice and the averaged acute data

plot_pivotPoints.m: script to plot 9 acute mice raw data and their pivot points
------------------------------------------------------------------------------------------
PROCESS OF AVERAGING
How we arrived at the data set in avg_lietal_2.csv
for all mice except 1 and 5
1) Spaced out mice data according to an absolute time (see sheet1)
2) I aligned according to the last (largest time, but first 2 rows in sheet2) 2 data points because these are where the shape matches. 
- Essentially this takes the mouse out of absolute time so we can average the rest of the data too.
3) Since each mouse is no longer in absolute time, I averaged glucose readings across each row 
4) To get the mouse back into absolute time (for plotting and mcmc purposes), I averaged times across each row as well
5) Note that because of this averaging, each time and glucose reading is not strictly increasing, so must perform sortrows when using data

* A more detailed description of this process can be found in the final REU write up under the section Problem Set Up


