

%============================
% ACUTE MICE ================
%============================

% open the csv file for acute mice
acute_df = textscan( regexprep( fileread('acute-age-onset.csv'), '\$', '0' ), '%f%f', 'delimiter', ',', 'HeaderLines', 1);
% retrieve only the ages from the second row
acute_ages = acute_df{2};

% get normal distribution from data
acute_dist = fitdist(acute_ages, 'Normal')



%============================
% PROGRESSIVE MICE ==========
%============================

% open the csv file for progressive mice
progressive_df = textscan( regexprep( fileread('progressive-age-onset.csv'), '\$', '0' ), '%f%f', 'delimiter', ',', 'HeaderLines', 1);
% retrieve only the ages from the second row
progressive_ages = progressive_df{2};

% get normal distribution from data
progressive_dist = fitdist(progressive_ages, 'Normal')

figure;
subplot(1,2,1)
% generate a histogram from age onset column, 15 bins, lognormal fit
acute_his = histfit(acute_ages, 15, 'lognormal');
% limit x-axis scale
set(gca, 'Xlim', [8, 30]);
% label x and y axises
ylabel('Frequencies');
xlabel('Day Of T1D Onset');
title('Acute Age Onset', 'fontsize', 20);
subplot(1,2,2)
% generate a histogram from age onset column, 30 bins, lognormal fit
progressive_his = histfit(progressive_ages, 30, 'lognormal');
% limit x-axis scale
set(gca, 'Xlim', [9, 30]);
% label x and y axises
ylabel('Frequencies');
xlabel('Day Of T1D Onset');
title('Progressive Age Onset', 'fontsize', 20);
