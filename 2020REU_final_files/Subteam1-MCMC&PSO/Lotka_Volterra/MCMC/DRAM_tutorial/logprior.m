%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Sections 5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB


% return log10(unnormalized prior)
% return log10(unnormalized prior)
function logPrior = logprior(theta, thetamu, thetasig)
% Uninfomative priors for params: uniform between 1e5 and 1
aPrior = unifpdf(theta(1), thetamu(1), thetasig(1)); 
bPrior = unifpdf(theta(2), thetamu(2), thetasig(2));
cPrior = unifpdf(theta(3), thetamu(3), thetasig(3));
dPrior = unifpdf(theta(4), thetamu(4), thetasig(4));
% Adding log priors == Multiplying non-log priors
% Gives log likelihood of entire candidate parameter set
logPrior = log(aPrior) + log(bPrior) + log(cPrior) + log(dPrior);
end