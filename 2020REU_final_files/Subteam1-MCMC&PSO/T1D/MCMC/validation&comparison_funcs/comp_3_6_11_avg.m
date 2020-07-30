%% Compare glucose predictions of DRAM parameterizations on Mice 3, 6, 11, and averaged data sets
%  Author:       M. Watanabe
%  Date:         July 2020
%  Desc:         Visual comparison of the glucose predictions from DRAM
%                parameterizations using the Mice 3, 6, 11, and averaged 
%                data sets

addpath('../writeUp_figs&results/workspaces');

% Load Mouse 3 data
load('jul13_mouse3_acute_NOD_waveOn_lietal.mat')
% Get prediction
meanParams = results.theta;
[t3, y3] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);
% Plot
plot(t3,y3(:,6), 'color', [0.31 0.58 0.83], 'LineWidth',2);
hold on
scatter(data(:,1), data(:,2), 40, [0.01 0.22 0.73], 'filled','o', 'LineWidth',2)
hold on

% Load Mouse 6 data
load('jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal.mat')
% Get prediction
meanParams = results.theta;
[t6, y6] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);
% Plot
plot(t6,y6(:,6), 'color', [0.87, 0.31, 0.31], 'LineWidth',2);
hold on
scatter(data(:,1), data(:,2), 40, [0.6 0.22 0.22], 'filled','o', 'LineWidth',2)
hold on

% Load Mouse 11 data
load('jul13_mouse11_acute_NOD_waveOn_lietal.mat')
% Get prediction
meanParams = results.theta;
[t11, y11] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);
% Plot
plot(t11,y11(:,6), 'color', [0.62 0.45 0.71], 'LineWidth',2);
hold on
scatter(data(:,1), data(:,2), 40, [0.41 0.3 0.47], 'filled','o', 'LineWidth',2)
hold on

% Load Averaged data
load('jul10_avg_run1(noIC)_acute_NOD_waveOn_lietal.mat')
% Get prediction
meanParams = results.theta;
[t, y] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);
plot(t,y(:,6), 'color', [0.46 0.72 0.58], 'LineWidth',2);
hold on
scatter(data(:,1), data(:,2), 40, [0.33 0.52 0.42], 'filled','o', 'LineWidth',2)
hold on

set(gca, 'FontSize', 15)
xlabel('day from birth');
ylabel('glucose (mg/dl)');
legend('Mouse 3 prediction', 'Mouse 3 raw', 'Mouse 6 prediction', 'Mouse 6 raw', 'Mouse 11 prediction', 'Mouse 11 raw','Average prediction', 'Average raw', 'Location', 'southeast')

