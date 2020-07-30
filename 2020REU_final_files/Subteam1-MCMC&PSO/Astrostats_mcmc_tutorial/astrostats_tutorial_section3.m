%% TRANSLATING MCMC TUTORIAL ASTROSTATS SECTION 3 ONLY
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Sections 3-5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB


%% 3) IS THIS COIN FAIR?
% In n coin tosses, we observe r heads

% 3.1 Uniform Prior

% Plot the posterior PDF
n = 20;
r = 7;
% Creating distribution of parameter p
p = linspace(0,1,201);
% creating binomal density function
pd = binopdf(r,n,p);
% Plot
plot(p,pd);

% normalizing parameter density
pdense = binopdf(r,n,p); % unnormalized in p
pmean = sum(p.*pdense)/sum(pdense);
% repeat the example: posterior PDFs for the probability p (uniform)
n2 = 20;
Nsamp = 201;
p2 = linspace(0,1,201);
deltap = 1/(Nsamp-1); % step size between samples of p
for r = 0:2:20
    pdense=binopdf(r,n2,p2);
    pdense=pdense/(deltap*sum(pdense)); %normalize posterior
    hold on; % plot all on one
    plot(p2, pdense);
    pmean=sum(p.*pdense)/sum(pdense);
end

% 3.2 Beta Prior
% same example with beta prior alpha=beta=10
p=linspace(0,1,201);
alphaprior=10;
betaprior=10;
% beta pdf
betapd=betapdf(p,alphaprior,betaprior);
plot(p,betapd);
% running the example again, reusing some variables from the previous
% example
n=20;
for r=0:2:20
    pdense=betapdf(p,alphaprior+r,betaprior+n-r);
    hold on;
    plot(p, pdense);
    p.mean=sum(p.*pdense)/sum(pdense);
    % Verify that this is the same as the direct calculation
    pdense2=binopdf(r,n,p)*betapdf(p,alphaprior,betaprior);
    pdense2=pdense2/(deltap*sum(pdense2)); % normalize posterior
end

% Plot the likelihood and posterior together
n=20;
alphaprior=20;
betaprior=20;
Nsamp=201; % no. of points to sample at
p=linspace(0,1,201);
deltap=1/(Nsamp-1); % step size between samples of p
prior=betapdf(p, alphaprior, betaprior);
for r=0:2:20
   like=binopdf(r,n,p);
   like=like/(deltap*sum(like)); % for plotting convenience in R, not sure how it translates in matlab
   post=betapdf(p,alphaprior,betaprior+n-r);
   plot(p,post);
   % add lines to represent like and prior
   hold on;
   plot(p, like);
   hold on;
   plot(p, prior);
end

% How the likelihood and posterior evolve as we get more data
alphaprior=10;
betaprior=10;
Nsamp=201; % num of points to sample at
p=linspace(0,1,Nsamp);
deltap=1/(Nsamp-1); % step size between samples of p
prior=dbeta(p,alphaprior,betaprior);
for i=1:9
   r=2^i;
   n=(3/2)*r;
   like=binopdf(r,n,p);
   like=like/(deltap*sum(like)); % for plotting convenience in R
   post=dbeta(p,alphaprior+r,betaprior+n-r);
   plot(p,post);
   % add lines to figure
   hold on;
   plot(p, like);
   hold on;
   plot(p, prior);
end % would still need to fix font size + add title to figures

