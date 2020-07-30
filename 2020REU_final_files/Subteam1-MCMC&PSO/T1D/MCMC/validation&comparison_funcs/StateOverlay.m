%% Biological check
%  Author:       M. Watanabe
%  Date:         July 2020
%  Desc:         Script to plot T1D state estimates using parameter values 
%                from Mice 3, 6, 11 and averaged data parameterizations in
%                order to determine biological feasibility of parameter
%                values.

% load past work spaces
clear all
addpath('../writeUp_figs&results/workspaces');
addpath('../helperFuncs');
addpath('../T1D_data_files/Lietal_data');

% Average
load('jul10_avg_run1(noIC)_acute_NOD_waveOn_lietal'); % red
params_avg = results.theta;
[tavg, yavg] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);

% Mouse 3
load('jul13_mouse3_acute_NOD_waveOn_lietal'); % green
params3 = results.theta;
[t3, y3] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);

% Mouse 6
load('jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal'); % magenta
params6 = results.theta;
[t6, y6] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);

% Mouse 11
load('jul13_mouse11_acute_NOD_waveOn_lietal'); % blue
params11 = results.theta;
[t11, y11] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);

% Baseline 
[tbase, ybase] = baselineSolver;

% Titles and plotting data
plotTitles{1}='Resting Macrophages';
units{1} = 'cells/ml';
plotTitles{2}='Activated Macrophages';
units{2} = 'cells/ml';
plotTitles{3}='Apoptotic beta cells';
units{3} = 'cells/ml';
plotTitles{4}='Necrotic beta cells';
units{4} = 'cells/ml';
plotTitles{5}='Healthy beta cells';
units{5} = 'mg';
plotTitles{6}='Glucose';
units{6} = 'mg/dl';
plotTitles{7}='Insulin';
units{7} = '\muU';
plotTitles{8}='Immunogenic DCs';
units{8} = 'cells/ml';
plotTitles{9}='Tolerogenic DCs';
units{9} = 'cells/ml';
plotTitles{10}='Effector T cells';
units{10} = 'cells/ml';
plotTitles{11}='Regulatory T cells';
units{11} = 'cells/ml';
plotTitles{12}='Memory T cells (Em)';
units{12} = 'cells/ml';

figure(1)
for i = 1:12
h=subplot(4,3,i);
    % Baseline
    plot(tbase, ybase(:,i), '--', 'LineWidth', 2, 'color', 'k');
    hold on 
    
    % Average
    plot(tavg, yavg(:,i), 'LineWidth',2, 'color', [0.46 0.72 0.58])
    hold on
    
    % Mouse 3
    plot(t3, y3(:,i), 'LineWidth', 2, 'color', [0.31 0.58 0.83])
    hold on
    
    % Mouse 6
    plot(t6, y6(:,i), 'LineWidth',2, 'color', [0.87, 0.31, 0.31])
    hold on

    
    % Mouse 11
    plot(t, y(:,i), 'LineWidth',2, 'color', [0.62 0.45 0.71])
   
    set(gca, 'YScale', 'log')
    set(gca, 'FontSize', 15);
    ylabel(units{i})
    title(plotTitles{i});
end
subplot(4,3,12)
legend('baseline', 'Average acute', 'Mouse 3', 'Mouse 6', 'Mouse 11')
han=axes(figure(1),'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('days from birth')
set(gca, 'FontSize', 20)
sgtitle('DRAM MCMC', 'FontSize', 20)
hold off

% plot states of interest for Mouse 6 for comparison with UKF
j=1;
figure(2)
for i = [3,7,10,11]
h=subplot(2,2,j);
    % Baseline
    plot(tbase, ybase(:,i), '--', 'LineWidth', 2, 'color', 'k');
    hold on
    % Mouse 6
    plot(t6, y6(:,i), 'LineWidth',2, 'color', 'm')
    hold on
     
    set(gca, 'YScale', 'log')
    set(gca, 'FontSize', 15);
    ylabel(units{i})
    title(plotTitles{i});
    
    j=j+1;
end

han2=axes(figure(2),'visible','off');
han2.XLabel.Visible='on';
han2.YLabel.Visible='on';
xlabel('days from birth')
set(gca, 'FontSize', 20)
subplot(2,2,4)
legend('baseline', 'fitted')
sgtitle('DRAM MCMC', 'FontSize', 20)
hold off
