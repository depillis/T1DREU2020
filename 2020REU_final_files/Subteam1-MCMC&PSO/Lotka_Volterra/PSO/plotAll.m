%% Plotting routine for Lotka-Volterra
%  Author:      Christina Catlett
%  Date:        July 20, 2020
%  Desc:        Used to plot best fits of all techniques against obs data;
%               comment/uncomment, change legend to plot subsets of data

function plotAll
% Load data files
load('HaresLynxData.mat');
load('dual.mat');
load('joint.mat');
load('dram.mat');
load('mh.mat');
load('PSOdata_SS.mat');
load('PSOdata_LL.mat');

% Define solution data
time = Lotka_Volterra_Data(:,1);
ydata = Lotka_Volterra_Data(:,2:end);
MH_sol = table2array(mhmcmcLV(:,2:end));
DRAM_sol = table2array(drammcmcLV(:,2:end));
Joint_sol = table2array(JointUKFSolPredatorPrey(:,2:end));
Joint_t = shiftTime(JointUKFSolPredatorPrey(:,1));
Dual_sol = table2array(DualUKFSolPredatorPrey(:,2:end));
Dual_t = shiftTime(DualUKFSolPredatorPrey(:,1));

figure(1); clf
tiledlayout(2,1);

%  Plot hares
ax1 = nexttile;
p = plot(ax1, time, ydata(:,1),'ok');
p.MarkerSize = 5;
p.MarkerFaceColor = 'k';
hold on;
plot(ax1, time, bfdataLL(:,1), '-g', 'LineWidth',2);
plot(ax1, time, bfdataSS(:,1), '-c', 'LineWidth',2);
plot(ax1, time, MH_sol(:,1), '-', 'LineWidth',2);
plot(ax1, time, DRAM_sol(:,1), '-m', 'LineWidth',2);
plot(ax1, Joint_t, Joint_sol(:,1), '-b', 'LineWidth',2);
plot(ax1, Dual_t, Dual_sol(:,1), '-r', 'LineWidth',2);

legend(ax1, {'Raw','PSO LL', 'PSO SS', 'MH MCMC', 'DRAM MCMC', 'Joint UKF', 'Dual UKF'}, 'Location', 'northwest');
set(gca,'FontSize',18);
xlabel(ax1, 'year');
ylabel(ax1, 'population (in thousands)');
title(ax1, 'Hare Population');

% Plot lynx
ax2 = nexttile;
q = plot(ax2, time, ydata(:,2),'ok');
q.MarkerSize = 5;
q.MarkerFaceColor = 'k';
hold on;
plot(ax2, time, bfdataLL(:,2), '-g', 'LineWidth',2);
plot(ax2, time, bfdataSS(:,2), '-c', 'LineWidth',2);
plot(ax2, time, MH_sol(:,2), '-', 'LineWidth',2);
plot(ax2, time, DRAM_sol(:,2), '-m', 'LineWidth',2);
plot(ax2, Joint_t, Joint_sol(:,2), '-b', 'LineWidth',2);
plot(ax2, Dual_t, Dual_sol(:,2), '-r', 'LineWidth',2);

legend(ax2, {'Raw','PSO SS', 'DRAM MCMC', 'Dual UKF'}, 'Location', 'northwest');
set(gca,'FontSize', 18)
xlabel(ax2, 'year');
ylabel(ax2, 'population (in thousands)');
title(ax2, 'Lynx Population');

