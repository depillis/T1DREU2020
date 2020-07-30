%% Compare effect of prior
%  Author:       M. Watanabe
%  Date:         July 2020
%  Desc:         Function to compare the predictions of DRAM 
%                parameterizations with and without an informative prior


addpath('../writeUp_figs&results/workspaces');

% load data for DRAM: uniformative prior
load('jul10_avg_run1(noIC)_acute_NOD_waveOn_lietal.mat')
meanParams = results.theta;
[t, y] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);

% load data for DRAM: informative (log-normal) prior
load('jul14_avg_wPriors_acute_NOD_waveOn_lietal.mat')
meanParams = results.theta;
[t2, y2] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);

% plot
figure(1); clf
scatter(data(:,1), data(:,2), 25, 'k', 'filled','o');
hold on
plot(t2,y2(:,6), 'color', 'r', 'LineWidth',2);
hold on
plot(t,y(:,6), 'color', 'b', 'LineWidth',2);


set(gca, 'FontSize', 15)
xlabel('day from birth');
ylabel('glucose (mg/dl)');
legend('Li data, averaged acute', 'Prediction with prior', 'Prediction no prior','Location', 'southeast');

    
