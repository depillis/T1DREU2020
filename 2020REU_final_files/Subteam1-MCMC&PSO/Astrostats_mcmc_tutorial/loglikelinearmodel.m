%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Sections 5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB

function llik = loglikelinearmodel(theta, obsdata, Ndat)
% loglike.linearmodel(theta, obsdata, ind=NaN)
%   Return log10(likelihood) for parameters theta for rows ind in obsdata
        %obsdata needs to be in table format
        %ind=1:size(obsdata,1); I think ind is only used in the R code for
        %a different way of matrix indexing which we can do differently in
        %matlab
        % convert alpha to a_1 and log10(ysig) to ysig
        theta(2)=tan(theta(2)); % taking the tan of the parameter a_1 (slope of linear model)
        theta(3)=10^theta(3); % undoing the log10 of the ysig (noise)

        dat=[ones(10,1) obsdata.Var1]; % in order to do matrix multiplication, need to add an extra column
        modPred=mtimes(theta(1:2),dat'); % plugging in parameters into data
        norm=log10(normpdf(modPred'-obsdata.y,0,theta(3))); 

        % creates a normal distribution from the difference bt the
        %   predicted model data and the observed data, with mean 0 and
        %   standard dev ysig
        %   evaluating the likelihood function (?) given the data
        llik = sum(norm); % scaled sum of the normal pdf : probability of a distribution
end

