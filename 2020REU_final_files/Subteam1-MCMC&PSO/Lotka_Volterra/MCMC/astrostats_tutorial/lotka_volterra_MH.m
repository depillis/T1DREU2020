%%Metropolis-Hastings MCMC for Lotka Volterra
% Authors: Maya Watanabe & Christina Catlett
%
% Date: 5/29/2020
% Description: R Code by C.A.L. Bailer-Jones from Sections 3-5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB
% Original R code intended for linear model. Here we adapt the
% Metropolis-Hastings MCMC algorithm to the Lotka-Volterra ODE model


%% Initialize the Lotka Volterra ODE Model

% Initialize data
load('HaresLynxData.mat');
obsdata=Lotka_Volterra_Data; 
idx = 60:91;
mytime = (1:length(idx))';
mydata(:,1) = obsdata(idx,2);
mydata(:,2) = obsdata(idx,3);
   
init=[mydata(1,1) mydata(1,2)];
y0(1) = mydata(1,1); y0(2) = mydata(1,2);

% Initialize a set of parameters (starting guess)
% thetaInit=[a b c d e]
% a: intrinsic rate of prey pop. increase
% b: predation rate coefficient
% c: reproduction rate of predators
% d: predator mortality rate

thetaInit=[0.6954 0.3865 0.5766 0.0150];

% Initialize covariance matrix
[sens, sensrel, flagg, y, sol] = senseq(thetaInit,idx, mydata);
sampleCov=inv(sens'.*sens)*(.01);

% metropolis-hastings mcmc algorithm
postSamp = metroMCMC(thetaInit, 0, 1e4, 1e3, sampleCov, obsdata);

% Create data and 2-by-1 tiled chart layout
xaxis = 1:1e4;
aplot = postSamp(:,3);
bplot = postSamp(:,4);
cplot = postSamp(:,5);
dplot = postSamp(:,6);

tiledlayout(4,1)

% a plot
ax1 = nexttile;
plot(ax1,xaxis,aplot)
title(ax1,'a density')
xlabel(ax1, 'Iteration')
ylabel(ax1,'a')

% b plot
ax2 = nexttile;
plot(ax2,xaxis,bplot)
title(ax2,'b per Iteration')
xlabel(ax2, 'Iteration')
ylabel(ax2,'b')

% c plot
ax3 = nexttile;
plot(ax3,xaxis,cplot)
title(ax3,'c per Iteration')
xlabel(ax3, 'Iteration')
ylabel(ax3,'c')

% d plot
ax4 = nexttile;
plot(ax4,xaxis,dplot)
title(ax4,'d per Iteration')
xlabel(ax4, 'Iteration')
ylabel(ax4,'d')
pause

% Compute a probability density estimate of the parameter data
[a, aden]=mvksdensity(postSamp(:,3));
[b, bden]=mvksdensity(postSamp(:,4));
[c, cden]=mvksdensity(postSamp(:,5));
[d, dden]=mvksdensity(postSamp(:,6));

% Plot densities

% Format data and 4-by-1 tiled chart layout
figure
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
[~,posMAP]=max(postSamp(:, 1) + postSamp(:,2)); % returns the maximum value of the posterior
thetaMAP=postSamp(posMAP, 3:5); % extract row#posMAP and columns 3-5 of postSamp:
%   extracting the parameters at the maximum value of the posterior 

thetaMean=[mean(postSamp(:,3)) mean(postSamp(:,4)) mean(postSamp(:,5)) mean(postSamp(:,6))]; %﻿Returns a vector or array or list of values obtained by applying a function to margins