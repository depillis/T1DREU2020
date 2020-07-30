How I arrived at the data set in avg_lietal_1.csv
for all mice except 1 and 5
1) Spaced out mice data according to an absolute time (see sheet1)
2) I aligned according to the last (largest time, but first 2 rows in sheet2) 2 data points because these are where the shape matches. 
- Essentially this takes the mouse out of absolute time so we can average the rest of the data too.
3) Since each mouse is no longer in absolute time, I averaged glucose readings across each row 
4) To get the mouse back into absolute time (for plotting and mcmc purposes), I averaged times across each row as well
5) Note that because of this averaging, each time and glucose reading is not strictly increasing, so must perform sortrows when using data