%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Section 5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB

% return Nsamp samples from prior
% NOTE: This uses a proper gamma prior on ysig instead of an improper
% Jeffreys prior, which is inconsistent with logprior.linearmodel(). This
% is generally a BAD IDEA and really we should modify
% logprior.linearmodel() to use proper priors.

% REFACTOR: sampleprior.linearmodel
function sampleprior = sampleprior(Nsamp)
a0 = normrnd(0, 4, [1, Nsamp]);
a1 = unifrnd(-pi/2, pi/2, [1, Nsamp]);
logysig = gamrnd(1.5, 1, [1, Nsamp]);
sampleprior = [a0, a1, logysig];
end

