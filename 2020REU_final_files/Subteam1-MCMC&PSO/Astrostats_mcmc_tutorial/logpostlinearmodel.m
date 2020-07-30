%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Sections 5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB

function [lplm] = logpostlinearmodel(theta, obsdata, Ndat)
% logpost.linearmodel(theta, obsdata, ind=NaN)
%   Return log10(unnormalized posterior) of the linear model
    lp = logprior(theta); % returns the log10 of the unnormalized prior
    if(isfinite(lp))
       lplm=loglikelinearmodel(theta, obsdata, Ndat) + lp; % evaluate model only if parameters are sensible
        % sum of the log of the likelihood and the log of the prior =
        % numerator of Bayes Theorem
    else
        lplm=NaN;
    end
end

