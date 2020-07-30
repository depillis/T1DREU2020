%   Written by Khoa T. Tran, School of EE & CS
%                            University of Newcastle
%        		             Australia.
%
%   Copyright (C) Khoa T. Tran
%% Initialization
clear; close all;clc  
global font_size
font_size=14;
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultUIControlFontName', 'Times New Roman')
% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', font_size)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', font_size)

%Use Latex as default typeset
set(0,'DefaultTextInterpreter','Latex') 

%White figure background
set(0,'DefaultFigureColor','w');
%Render plot
set(0,'DefaultFigureRenderer', 'painters');
set(0,'DefaultLineLineWidth',1.1);
set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-|--|:|-.|.')
%% Provide the target log-density function, e.g. Multivariate Gaussian
var  = 1e-2;
logpdf=@(X)...                          % log-density function handle
    -0.5*(sum((X).^2)/var + length(X)*log(var*2*pi));              
%% Provide the initial state X and initial covariance matrix V
Dims = 10;X = zeros(Dims,1);            % User input
V = abs(var+1e-1*var*randn)*eye(Dims);  % User input
%% Check covariance matrix V for positive definiteness
L = tril(V,-1);V=diag(diag(V))+L+L';    % Fix symmetry
try
    V2=chol(V,'lower');
catch % Fix positive definiteness
    [V3,D] = eig(V);d=diag(D);
    d(d<=0)=10*realmin;V= V3*diag(d)*V3';
    L = tril(V,-1);V=diag(diag(V))+L+L'; % Fix symmetry again
    V2=chol(V,'lower');
end
%% Try to throw 1 million X's to initialize the chain inside the density support
if logpdf(X)==-Inf
    fprintf('Initial state is outside density support.\nThrowing random states for initilisation...\n')
    throws=1;X0=X;
    while logpdf(X)==-Inf && throws < 1e6
        X = X0+V2*randn(size(X));throws=throws+1;
    end
    if throws>=1e6
        error('Cannot initialize inside the support of the target distribution')
    end
end
%% Specify Options for MCMC Sampling
OPT.logpdf = logpdf;    % Handle to log-density function
OPT.V = V;              % Covariance matrix
OPT.Dims = Dims;        % Number of dimensions
OPT.dsp = 1;            % Display status of MCMC runs, default value = 1
OPT.Mmax=1e5;           % Number of MCMC samples, recommend 1e6 sample for slice sampling
OPT.div = ...           % Divide into smaller sestions to save RAM
    max(20,ceil(OPT.Mmax/1e6));
OPT.hgrd= 1e2;          % Number of grid point in histograms, must be multiple of 10
%% Perform MCMC sampling
sampler = -1;
while isempty(sampler)||~isa(sampler,'double')||sampler>2||sampler<0
    sampler = input(['Specify type of MCMC sampler by entering:\n' ...
                    ' 0 for Factor Slice Sampling\n' ...
                    ' 1 for Adaptive Metropolis\n' ...
                    ' 2 for Random Walk Metropolis\n']);
end
sampler = round(sampler);
switch sampler
    case 0 % Factor Slice Sampling
        OPT.filename='FSS_output';
        G = FSS(X,OPT);
    case 1 % Adaptive Metropolis 
        OPT.filename='AM_output';
        G = AM(X,OPT);
    case 2 % Random Walk Metropolis
        OPT.filename='RWM_output';
        G = RWM(X,OPT); 
end
%% Show raw TH data 
% Load back the last batch of MCMC samples to RAM
D = load([G.filename '_' num2str(OPT.div) '.mat']); % Load TH data
G.TH = D.TH; 
% Plotting histograms
for k=1:(Dims)
    figure
    subplot(1,2,1)
    plot(G.TH(k,:))
    line([0 size(G.TH,2)],G.best_TH(k)*[1 1],'color','g')
    ylabel(['$$\theta_{' num2str(k) '}$$'])
    xlabel('Iteration')
    axis tight
    subplot(1,2,2)
    counts = G.hist.p(k,:);
    bins = G.hist.x(k,:);
    barh(bins,counts)
    line([0 max(counts)],G.best_TH(k)*[1 1],'color','g')
    ylabel(['$$\theta_{' num2str(k) '}$$'])
    xlabel('Probability density')
    axis tight
end