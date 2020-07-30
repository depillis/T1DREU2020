%% Example: fitting Lotka-Volterra with PSO
%  Author:       Christina Catlett
%  Date:         June 9, 2020
%  Desc:         Main function to fit coefficients a, b, c, d, sigmaH, sigmaL to Lotka Volterra
%                system for 'HareLynxData.mat'. Assumes uniform priors, Gaussian likelihood,
%                independent variance for H & L residuals if log-likelihood obj func is chosen.
%                Requires Global Optimization Toolbox.

% Reset at beginning of run
clear all
rng default 

% Load datasets
load('HaresLynxData.mat')    % Loads as 'Lotka_Volterra_Data'
load('LVsimData.csv')        % True Params: [.7 .1 .7 .1 0 0]
load('LVsimData_wnoise.csv') % True Params: [.7 .1 .7 .1 .05 .05]

% Select dataset: specify data from above files, retype as string
data = Lotka_Volterra_Data;
data_str = 'Lotka_Volterra_Data'; % Formatting for switch-statement

% Select objective function: ss -> sum of squares, ll -> normal log-likelihood
objfunchoice = 'ss';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Format particleswarm input, options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch objfunchoice
    case 'll'
        npars = 6; % Number of parameters being optimized

        % Set lower, upper bounds for pars: a, b, c, d, sigmaH, sigmaL
        switch data_str
            case 'LVsimData'
                lb = [.1; 5e-3; .1; 5e-3; 0; 0];
                ub = [1; .2; 1; .2; 1; 1];
            case 'LVsimData_wnoise'
                lb = [.1; 5e-3; .1; 5e-3; 0; 0];
                ub = [1; .2; 1; .2; 1; 1];
            case 'Lotka_Volterra_Data'
                lb = [.1; 5e-3; .1; 5e-3; 0; 0];
                ub = [1; .2; 1; .2; 5000; 1000];
        end
        
        % Objective func -> Maximize log likelihood by minimizing |loglikelihood|
        fn = @(params) LVabsllike(params, data);
    case 'ss'
        npars = 4; % Number of parameters being optimized

        % Set lower, upper bounds for pars: a, b, c, d
        lb = [.1; 5e-3; .1; 5e-3];
        ub = [1; .2; 1; .2];
        
        % Objective func -> Minimize sum of squares difference
        fn = @(params) LVss(params, data);  
end

options = optimoptions('particleswarm', 'SwarmSize', 100, 'MaxIterations', 500, ...
    'UseParallel', true, 'Display', 'iter', 'DisplayInterval', 10, 'PlotFcn', @pswplotbestf, 'MinNeighborsFraction', 1);
options = optimoptions(options,'PlotFcn',@pswplotbestf); % Plot iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run particleswarm, plot iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out, fval, exitflag, output] = particleswarm(fn, npars, lb, ub, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overplot best fit on data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate ODE solved for best fit params
bfdata = LVfun([0:length(data(:,1)) - 1], out, data(:,2:end));

figure(1); clf
tiledlayout(2,1)
xaxis = data(:,1)';

ax1 = nexttile;
plot(ax1, xaxis, data(:,2),'sb');
hold on;
title(ax1, 'PSO best fit for observed hare population');
plot(ax1, xaxis, bfdata(:, 1), '-k');
xlabel(ax1, 'year');
ylabel(ax1, 'population');

ax2 = nexttile;
plot(ax2, xaxis, data(:,3),'sb');
hold on;
title(ax2, 'PSO best fit for observed lynx population');
plot(ax2, xaxis, bfdata(:, 2), '-k');
xlabel(ax2, 'year');
ylabel(ax2, 'population');
