% This function will calibrate your model to the given data via
% parameter estimation using the DRAM algorithm. You will need to have the
% mcmcstat folder added to your current path for this to run.

% Inputs:
%   ODEModel: specifies the ode's of the model being fit
%   params: parameter structure from paramStruct()
%   data: selected data from dataSelection()

% Outputs:
%   fitParams: optimized parameters for model fit
%   results: a structure containing info from DRAM
%   chain: a (# of simulations x # of params) matrix containing MCMC chains
%   s2chain: a (# of simulations x 1) vector containing error variance chain
%   ss2chain: a (# of simulations x 1) vector containing likelihood
%   posterior chain

%Call parameters
params=paramStruct;

%Load data structure

load algaedata.mat;

% The model sum of squares in file <algaess.html |algaess.m|> is
% given in the model structure.
model.ssfun = @algaess;

%%
% We assume having at least some prior information on the
% repeatability of the observation and assign rather non informational
% prior for the residual variances of the observed states. The default
% prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (see for example Gelman et al.). The 3
% components (_A_, _Z_, _P_) all have separate variances.
model.S20 = [1 1 2];
model.N0  = [4 4 4];

%%
% First generate an initial chain.
options.nsimu = 1000;
[results, ~, ~, ~]= mcmcrun(model,data,params,options);
%%
% Then re-run starting from the results of the previous run,
% this will take couple of minutes.
options.nsimu = 5000;
[results, chain, s2chain, ss2chain] = mcmcrun(model,data,params,options, results);

% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
fitParams = chain(ind,:); %These are the fitted parameter values


%%
% Chain plots should reveal that the chain has converged and we can
% use the results for estimation and predictive inference.
figure(2); clf
mcmcplot(chain,[],results,'pairs');
figure(3); clf
mcmcplot(chain,[],results,'denspanel',2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = algaess(theta,data)

%extract the data
time   = data.ydata(:,1);
ydata  = data.ydata(:,2:end);
xdata  = data.xdata;

% 3 last parameters are the initial states
y0 = theta(end-2:end);

ymodel = algaefun(time,theta,y0,xdata);

ss = sum((ymodel - ydata).^2);
end

function y=algaefun(time,theta,y0,xdata)
% algae model function

[t,y] = ode15s(@algaesys,time,y0,[],theta,xdata);

end

function ydot = algaesys(t,y,theta,xdata)
% ode system function for MCMC algae example

A = y(1);
Z = y(2);
P = y(3);

% control variables are assumed to be saved
% at each time unit interval
QpV = xdata(ceil(t),2);
T   = xdata(ceil(t),3);
Pin = xdata(ceil(t),4);

mumax = theta(1);
rhoa  = theta(2);
rhoz  = theta(3);
k     = theta(4);
alpha = theta(5);
th    = theta(6);

mu = mumax*th^(T-20)*P/(k+P);

dotA = (mu - rhoa - QpV -alpha*Z)*A;
dotZ = alpha*Z*A - rhoz*Z;
dotP = -QpV*(P-Pin) + (rhoa-mu)*A + rhoz*Z;

ydot=[dotA;dotZ;dotP];
end

%end