%% Lotka-Volterra Model
% Example code from https://mjlaine.github.io/mcmcstat/ex/algaeex.html
% Edited to fit Lotka-Volterra model by Christina Catlett & Maya Watanabe (June 2020)
% Uses functions lotkaVolterrasys.m, lotkaVolterrafun.m, lotkaVolterrass.m,
% and mcmcstat library

%% Loading and plotting observed data

clear model data params options
load HaresLynxData.mat

figure(1); clf
plot(Lotka_Volterra_Data(:,1),Lotka_Volterra_Data(:,2:end),'o-');
title('Observed Populations Over Time');
legend({'Hares', 'Lynx'},'Location','best');
xlabel('years');


%% Calculating model sum of squares

model.ssfun = @(theta, data) lotkaVolterrass(theta, Lotka_Volterra_Data);

%% Defining intial parameters
% All parameters are constrained to be positive, uniformly distibuted.
% Initial parameter values from PSO
% {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
params = {
    {'alpha', 0.525771872777156, 0, 1}
    {'beta',  0.0202508516159993, 0, 1}
    {'gamma', 0.888556080190351, 0, 1}
    {'delta', 0.0223647386536971, 0, 1}
    };

%% Defining initial variance 
% We assume having at least some prior information on the
% repeatability of the observation and assign rather non informational
% prior for the residual variances of the observed states. The default
% prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (see for example Gelman et al.). The predator and
% prey components have separate variances.

[mse, minparams] = fitInitialParams(params, Lotka_Volterra_Data);
model.sigma2 = mse;
%params = minparams;


%% Burn-in iterations
% First generate an initial chain.

options.nsimu = 1e4;
options.method = 'dram';
[results, chain, s2chain, ss2chain]= mcmcrun(model,Lotka_Volterra_Data,params,options);
%% Running MCMC chain
% Then re-run starting from the results of the previous run,

options.nsimu = 5e4;
options.method = 'dram';
[results, chain, s2chain, ss2chain] = mcmcrun(model,Lotka_Volterra_Data,params,options,results);

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

chainstats(chain,results)

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
out = mcmcpred(results,chain,s2chain,Lotka_Volterra_Data,modelfun,nsample);
figure(4); clf
mcmcpredplot(out);

% add the 'y' observations to the plot
hold on
for i=1:2
  subplot(2,1,i)
  hold on
  plot(Lotka_Volterra_Data(:,1),Lotka_Volterra_Data(:,i+1),'s');
  titles = ['Hares', 'Lynx'];
  ylabel(''); title(titles(i));
  hold off
end
xlabel('year');
