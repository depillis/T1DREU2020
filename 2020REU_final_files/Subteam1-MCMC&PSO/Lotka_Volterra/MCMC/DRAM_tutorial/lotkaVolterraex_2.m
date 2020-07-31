%% Lotka-Volterra Model
% Authors:      Edited by C. Catlett & M. Watanabe
%               Original code from https://mjlaine.github.io/mcmcstat/ex/algaeex.html
% Date:         June 2020
%
% Descr:        Script to estimate parameters of the Lotka-Volterra system.
%               This script can be used to perform the Metropolis-Hastings
%               and DRAM MCMC methods. Specify 'mh' or 'dram' in lines 87
%               and 99.
%               Uses functions lotkaVolterrasys.m, lotkaVolterrafun.m,
%               lotkaVolterrass.m, and mcmcstat library
%
% Directions:   This script is meant to be run in full, ~15 min. User
%               specifications exist for DRAM vs Metropolis-Hastings and
%               naming the workspace.

%% Loading and plotting observed data

clear model data params options
addpath('mcmcstat'); % add mcmcstat library

% Specify a unique file name to save work space at conclusion of the
% parameter estimation
filename = 'dram_LV_nopriors';

% SELECT DATA
% observed predator-prey data: d = 0
% simulated data:              d = 1
d = 0;

if d==0
    load HaresLynxData.mat
    data = Lotka_Volterra_Data;
else
    load makeSimData.mat
    data = simData;
end

% Plot original data
figure(1); clf
plot(data(:,1),data(:,2:end),'o-');
title('Observed Populations Over Time');
legend({'Hares', 'Lynx'},'Location','best');
xlabel('years');


%% Calculating model sum of squares
% This is the likelihood function

model.ssfun = @(theta, data) lotkaVolterrass(theta, data);

%% Defining intial parameters
% All parameters are constrained to be positive, uniformly distibuted.

% {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
params = {
    {'alpha', 0.7, 0, 1}
    {'beta',  0.1, 0, 1}
    {'gamma', 0.7, 0, 1}
    {'delta', 0.1, 0, 1}
    };
%% Defining initial variance
% We assume having at least some prior information on the
% repeatability of the observation and assign rather non informational
% prior for the residual variances of the observed states. The default
% prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (see for example Gelman et al.). The predator and
% prey components have separate variances.

% Minimize L-V ODE to find best initial parameter values
[mse, minparams] = fitInitialParams(params, data);
model.sigma2 = mse;
params = minparams;

% Option to specify an informative prior function
% model.priorfun = @(theta, thetamu, thetasig) logprior(theta, thetamu, thetasig);

% Specify mean and var of prior distribution
model.S20 = [1];
model.N0 = [1];


%% Burn-in iterations
% First generate an burn-in initial chain.

options.nsimu = 1e4;
options.method = 'dram';
[results, chain, s2chain, ss2chain]= mcmcrun(model,data,params,options);

% plot burn-in chain
plot(chain(1:options.nsimu,:));
set(gca, 'FontSize', 20)
title('Burn-in Parameter Chain Values');
legend('alpha', 'beta', 'gamma', 'delta');
%% Running MCMC chain
% Then re-run starting from the results of the previous run,

options.nsimu = 5e5;
options.method = 'dram';
[results, chain, s2chain, ss2chain] = mcmcrun(model,data,params,options,results);

% plot chain
plot(chain(1:options.nsimu,:));
set(gca, 'FontSize', 20)
title('Parameter Chain Values');
legend('alpha', 'beta', 'gamma', 'delta');
%% Plotting chain results, densities
% Chain plots should reveal that the chain has converged and we can
% use the results for estimation and predictive inference.

figure(2); clf
mcmcplot(chain,[],results,'pairs');
figure(3); clf
mcmcplot(chain,[],results,'denspanel',2);

%% Compute statistics for resultant parameters
% Function |chainstats| calculates mean ans std from the chain and
% estimates the Monte Carlo error of the estimates. Number |tau| is
% the integrated autocorrelation time and |geweke| is a simple test
% for a null hypothesis that the chain has converged.

chainstats
cs = chainstats(chain,results);
pName = {'alpha'; 'beta'; 'gamma'; 'delta'};
mean = cs(:,1);
std = cs(:,2);
MC_err = cs(:,3);
tau = cs(:,4);
geweke = cs(:,5);
T = table(pName, mean, std, MC_err, tau, geweke);
writetable(T, strcat(filename, 'chainstats.txt'))

% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));
ind = ind(1);
fitParams = chain(ind,:); %Fitted parameter values


%% Plotting prediction vs. data

% In order to use the |mcmcpred| function we need
% function |modelfun| with input arguments given as
% |modelfun(xdata,theta)|. We construct this as an anonymous function.
modelfun = @(d,th) lotkaVolterrafun(d(:,1),th,d(:,2:3));

% We sample 500 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.

nsample = 200;
out = mcmcpred(results,chain,s2chain,data,modelfun,nsample);

% plot predator-prey predictions w/95% confidence interval
figure(4); clf
mcmcpredplot(out);
set(gca, 'FontSize', 20)
% add the 'y' observations to the plot
hold on
for i=1:2
  subplot(2,1,i)
  hold on
  scatter(data(:,1),data(:,i+1), 30, 'b', 'filled', 'o');
  set(gca, 'FontSize', 20)
  ylabel('');
  hold off
end
xlabel('year');
hold on
subplot(2,1,1)
title('Predicted Hare Population')
subplot(2,1,2)
title('Predicted Lynx Population')
%% Evaluation of Parameterization using mean squared error
meanparams = results.theta;
[mse_prey mse_pred] = LVmse(meanparams, data);

predmod = lotkaVolterrafun(data(:,1), meanparams, data(:,2:end));

for i=1:2
  subplot(2,1,i)
  hold on
  plot(data(:,1),data(:,i+1),'s');
  hold on
  plot(data(:,1), predmod(:,i));
  set(gca, 'FontSize', 20)
  xlabel('year');
  ylabel('population'); %title(titles{i});
end

hold on
subplot(2,1,1)
title('Mean prediction for Hare population')
subplot(2,1,2)
title('Mean prediction for Lynx population')

%% Save workspace variables
save(filename);











