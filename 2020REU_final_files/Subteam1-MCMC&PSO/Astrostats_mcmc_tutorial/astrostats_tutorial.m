%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Sections 3-5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB


%% 5: Fitting a Linear Model with unknown noise

% 1) Define true model and simulate experimental data

rng(50);% set random number generator
Ndat=10;%  number of data points
x=sort(unifrnd(0,1:Ndat)); % sorted vector of 10 numbers between 0 and 10
sigTrue=1; % standard deviation of the noise of the simulated data y (113)
modMat=[0 1]; % row vector
% long way to create this 2d matrix, need a better way
ytemp=[ones(10,1) x']; % Ndat by 2 matrix with variable x data from line 108 in the second column
y=mtimes(ytemp, modMat') + normrnd(0, sigTrue, Ndat, 1); % final data (points) with noise added

% next they import some data that we don't have access to --> this is the
% TRUE MODEL

% True parameters transformed 
thetaTrue=[modMat(1) atan(modMat(2)) log10(sigTrue)]; % transformed parameters of the true data: 2 real parameters a_0, a_1, and sigTrue is the noise
% why are the parameters arctan transformed?
obsdata=table(x', y); % table with x: uniform randomly generated points, y: transformed uniform data with noise 
obsdataMat=obsdata{:,:}; % matrix version of obsdata

% 2) Define Model and infer posterior PDF over its parameters
% Model to infer: linear regression with Gaussian noise
% Parameters: intercept a_0, gradient a_1; Gaussian noise sigma: ysig

% ysig = Jeffreys prior used when you don't have a suitable prior
%   distribution; has little impact on the posterior

% theta is a 1xJ vector of all model parameters (J=3); the sampling is
%   performed on theta defined to make sampling from a Gaussian more convenient as
%   1)linear steps in alpha better than in tan(alpha) - What are these
%   Linear steps?
%   2)ysig cannot be neg and additive steps in log10(ysig) - the noise
%   cannot be negative

% define covariance matrix of MCMC sampleing PDF:
sampleCov=makecovmatrix([0.1 0.02 0.1],1); % matrix to estimate parameter distribution
% related to the prior? evalutaed prior?
% set starting point
thetaInit=[2 pi/8 log10(3)]; % initial guess for parameters

% run the MCMC to find postSamp: samples of the posterior PDF
rng(150);

% metropolis-hastings mcmc algorithm
postSamp = metroMCMC(thetaInit, 0, 1e4, 1e3, sampleCov,  obsdata, Ndat);
% metroMCMC(func - function, thetaInit-inital guess for params, Nburnin - 0, Nsamp - 1e4,sampleCov, obsdata  
% why is the burnin 0?
% returns a 1e3 by 5 matrix:
%   first 2 columns contain the posteriors for each sample in the algorithm
%   second 3 columns contain each of the 3 parameters for each sample

% plot MCMC chains & use density estimation to plot 1D posterior PDFs 
% plot postSamp

figure
% Create data and 2-by-1 tiled chart layout
xaxis = 1:1e4;
a0plot = postSamp(:,3);
aplot = postSamp(:,4);
sigplot = postSamp(:,5);
tiledlayout(3,1)

% a0 plot
ax1 = nexttile;
plot(ax1,xaxis,a0plot)
title(ax1,'a_0 density')
xlabel(ax1, 'Iteration')
ylabel(ax1,'a_0')

% a plot
ax2 = nexttile;
plot(ax2,xaxis,aplot)
title(ax2,'a per Iteration')
xlabel(ax2, 'Iteration')
ylabel(ax2,'arctan(a)/rad')

% sig plot
ax3 = nexttile;
plot(ax3,xaxis,sigplot)
title(ax3,'sigma per Iteration')
xlabel(ax3, 'Iteration')
ylabel(ax3,'log10(sigma)')
pause

% Compute a probability density estimate of the parameter data
[a0 a0den]=mvksdensity(postSamp(:,3));
[a1,a1den]=mvksdensity(postSamp(:,4));
[ysig,ysigden]=mvksdensity(postSamp(:,5));
% plot the densities
figure
% Create data and 2-by-1 tiled chart layout
x1axis = a0den;
x2axis=a1den;
x3axis=ysigden;
a0denplot = a0;
a1denplot = a1;
ysigdenplot = ysig;
tiledlayout(3,1)

% a0 plot
ax1 = nexttile;
plot(ax1,x1axis,a0denplot)
title(ax1,'a_0 density')
xlabel(ax1, 'a_0')
ylabel(ax1,'density')
xline(thetaTrue(1));

% a plot
ax2 = nexttile;
plot(ax2,x2axis,a1denplot)
title(ax2,'a_1 density')
xlabel(ax2, 'alpha = arctan(a_1)/radians')
ylabel(ax2,'density')
xline(thetaTrue(2));

% sig plot
ax3 = nexttile;
plot(ax3,x3axis,ysigdenplot)
title(ax3,'sigma density')
xlabel(ax3, 'log10(ysig)')
ylabel(ax3,'density')
xline(thetaTrue(3));
pause

% Find MAP solution and mean solution
% MAP = Maximum A Posteriori - peak of posterior (not peak in 1D PDF, but
% peak of the 3D PDF)
[~,posMAP]=max(postSamp(:, 1) + postSamp(:,2)); % returns the maximum value of the posterior
thetaMAP=postSamp(posMAP, 3:5); % extract row#posMAP and columns 3-5 of postSamp:
%   extracting the parameters at the maximum value of the posterior 

thetaMean=[mean(postSamp(:,3)) mean(postSamp(:,4)) mean(postSamp(:,5))]; %﻿Returns a vector or array or list of values obtained by applying a function to margins of an array or matrix.;

% Overplot these solutions with original data and true model
errorbar(obsdataMat(:,1), obsdataMat(:,2), sigTrue);
% add models
refline(modMat(1), modMat(2)); % true model
refline(thetaMAP(1), tan(thetaMAP(2))); % MAP model
refline(thetaMean(1), tan(thetaMean(2))); % mean model

% Compare this with the result from ML estimation from lm()
refline(fitlm(obsdataMat(:,1), obsdataMat(:, 2)));

%%% Make prediction: determine PDF(ycand | xnew, obsdata)

% Set up uniform grid of candidate values of y, ycand[], and
% at each of these, calculate the probability density function (a scalar)
% by integrating the likelihood over the posterior. We do this with the
% Monte Carlo approximation of the integration. Model and likelihood used
% here must be consistent with logpost.linearmodel() !

xnew = 6; % then repeat with xnew = 15
dy = 0.01;
ymid = thetaMAP(1) + xnew*tan(thetaMAP(2)); % for centering choice of ycand
ycand = ymid-5:dy:ymid+5; % uniform sampling of y with step size dy
ycandPDF = zeros(1,length(ycand),'double');

% predict y at all values of paramters drawn from post. by applying model.
% dimensions in matrix mult: [Nsamp x 1] = [Nsamp x P] * [P x 1]

modPred = [postSamp(:,3), tan(postSamp(:,3))] * transpose([1, xnew]);
for k = 1:length(ycand)
    like = normpdf(modPred-ycand(k), 0, 10^postSamp(:,5)); % [Nsamp x 1]
    ycandPDF(k) = mean(like);
end

% Note that ycandPDF[k] is normalized, i.e. sum(dy*ycandPDF) = 1
figure
plot(ycand, ycandPDF)
% Find peak and approximate confidence intervals at 1sigma
% REFACTOR: peak.ind, lower.ind, upper.ind
[~,peak] = max(ycandPDF);
[~,lower] = max(find(cumsum(dy*ycandPDF) < normpdf(1)));
[~,upper] = max(find(cumsum(dy*ycandPDF) < normpdf(-1)));

% Overplot this prediction with original data and the models

