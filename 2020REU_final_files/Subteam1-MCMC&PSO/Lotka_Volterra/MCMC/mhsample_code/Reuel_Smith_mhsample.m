%% Matlab Routine to Perform Bayesian Parameter Estimation of Lotka-Volterra ODE model
% Authors: Maya Watanabe & Christina Catlett 
% June 2020
% Based on code by Reuel Smith 2015-2017
% v. MATLAB R2015b through 2017a
% Description: This file implements the MATLAB mhsample function using an
% external loglikelihood function. It plots the resulting MCMC chain,
% calculates and plots the parameter densities
% ========================================================================
% Initialize data and constants
% ========================================================================
% Definitions of different variables used in this code
% Using the definition of the Lotka-Volterra predator-prey model, the
% parameters are:

% theta(1): A - intrinsic rate of prey population increase
% theta(2): B - predation rate coefficient
% theta(3): C - reproduction rate of predators per 1 prey eaten
% theta(4): D - predator mortality rate

%
% load the data
load('HaresLynxData.mat');
obsdata=Lotka_Volterra_Data; 

% Establish an initial parameter guess
% this thetaInit was found using D. Shenker & R. Wander's PSOPC_LVModel.m 
%  code found at: 
%  https://github.com/shtyllab/2020-HMC-REU-Codes/tree/master/GA_Parameter_Identification
%thetaInit=[0.6954 0.3865 0.5766 0.0150];

% ========================================================================
% Initialize constants
% ========================================================================
% initial parameter estimate
thetaInit=[.7 .1 .7 .1];

% number of parameters in likelihood equation)
n=4;

% sigma for the proposed logPDF, an nxn identity matrix
prop_sig=eye(n);

% number of randomly generated samples
nsamples = 20000;

% Establish burn-in period
K = 1000;

% This parameter controls the size of the new markov chain which omits M-1
% out of M values.  This will curb the effect of autocorrelation.
% M-1 out of M values omitted in the generated sequence
M = 10;



% =========================================================================
% Implement proposal and target distributions
% =========================================================================
% DEFINE THE PROPOSAL DISTRIBTUION (logproppdf) and random number generator
% (proprnd).  The standard deviation for each parameter in the proprnd is
% to be chosen such that the updated guess covers a reasonable range of
% estimates.  These should be small if the initial estimate is a close one.
logproppdf = @(x,y) log10(mvnpdf(x,y,prop_sig)); % lognormal proposal dist.
proprnd = @(y) [normrnd(y(1),0.015),normrnd(y(2),0.04),...
    normrnd(y(3),0.0211),normrnd(y(4),0.04)]; % random sample from proppdf

% DEFINE THE LIKELIHOOD X PRIORS (posteriors)
% The prior distributions for the parameters are assumed to be
% non-informative, therefore they are treated as uniform distributions
% between 1e-5 and 1. These priors multiplied by the product
% of the lognormal PDF of each data point form the LikelihoodxPriors term
% that undergoes Bayesian updating.
% LL_mhsample is the function in which we compute the loglikelihood
% function
logpdf = @(x) log10(unifpdf(x(1),1e-5, 1))+log10(unifpdf(x(2),1e-5, 1))+...
    log10(unifpdf(x(3),1e-5, 1))+log10(unifpdf(x(4),1e-5, 1))+...
    (LL_mhsample(thetaInit, obsdata));

% =========================================================================
% MH-sampling routine
% =========================================================================
[result,accept] = mhsample(thetaInit,nsamples,'logpdf',logpdf,...
    'logproppdf',logproppdf,'proprnd',proprnd,'burnin',K,'thin',M);



% Alternatively the BURN-IN, and THIN value may be left out for a faster,
% but less precise MH-sampling procedure
% [result,accept] = mhsample(init_param,nsamples,'pdf',pdf,'proppdf',proppdf,'proprnd',proprnd);


% =========================================================================
% Plot the results
% =========================================================================

% plot the MCMC chain
figure(1)
plot(result)

% plot histograms of parameter values
figure(2)
subplot(1,4,1)
hist(result(:,1),100)
xlabel('Parameter A')
subplot(1,4,2)
hist(result(:,2),100)
xlabel('Parameter B')
subplot(1,4,3)
hist(result(:,3),100)
xlabel('Parameter C')
subplot(1, 4, 4)
hist(result(:,4),100)
xlabel('Parameter D')
set(gcf,'color','w');

%calculate the densities
[a, aden]=mvksdensity(result(:,1));
[b, bden]=mvksdensity(result(:,2));
[c, cden]=mvksdensity(result(:,3));
[d, dden]=mvksdensity(result(:,4));

% plotting the densities
figure(3)
x1axis = aden;
x2axis = bden;
x3axis = cden;
x4axis = dden;
adenplot = a;
bdenplot = b;
cdenplot = c;
ddenplot = d;
tiledlayout(4,1)

% a plot
ax1 = nexttile;
plot(ax1,x1axis,adenplot)
title(ax1,'a density')
xlabel(ax1, 'a')
ylabel(ax1,'density')


% b plot
ax2 = nexttile;
plot(ax2,x2axis,bdenplot)
title(ax2,'b density')
xlabel(ax2, 'b')
ylabel(ax2,'density')


% c plot
ax3 = nexttile;
plot(ax3,x3axis,cdenplot)
title(ax3,'c density')
xlabel(ax3, 'c')
ylabel(ax3,'density')


% d plot
ax4 = nexttile;
plot(ax4,x4axis,ddenplot)
title(ax4,'d density')
xlabel(ax4, 'd')
ylabel(ax4,'density')
pause

% Find MAP solution and mean solution
% MAP = Maximum A Posteriori - peak of posterior (not peak in 1D PDF, but
% peak of the 3D PDF)
max1=max(result(:,1));
max2=max(result(:,2));
max3=max(result(:,3));
max4=max(result(:,4));
thetaMAP=[max1 max2 max3 max4];
 % extract row#posMAP and columns 3-5 of postSamp:
%   extracting the parameters at the maximum value of the posterior 

thetaMean=[mean(result(:,1)) mean(result(:,2)) mean(result(:,3)) mean(result(:,4))]; %﻿Returns a vector or array or list of values obtained by applying a function to margins