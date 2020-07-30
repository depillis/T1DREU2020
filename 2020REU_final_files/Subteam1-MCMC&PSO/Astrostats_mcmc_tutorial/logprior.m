%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Sections 5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB


% return log10(unnormalized prior)
function logPrior = logprior(theta)
a0Prior = normpdf(theta(1), 0, 2);
alphaPrior = 1;
logysigPrior = 1;
logPrior = log10(a0Prior) + log10(alphaPrior) +log10(logysigPrior);
end