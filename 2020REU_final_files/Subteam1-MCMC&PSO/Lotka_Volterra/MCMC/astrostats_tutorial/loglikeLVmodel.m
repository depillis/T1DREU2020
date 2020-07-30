%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Sections 5 of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB

function llik = loglikeLVmodel(theta, obsdata)
%   Return log10(likelihood) for parameters theta over all of obsdata
        %time_range_cont = [obsdata(1,1) obsdata((length(obsdata(:,1))),1)]; % Continuous time range for data
        %time_range_disc = obsdata(1,1):1:obsdata((length(obsdata(:,1))),1); % Discrete time range for data
        obspops = [obsdata(:,2) obsdata(:,3)]; % Collected data for population size, 2x90 matrix
        sampleCov = cov(obsdata(:,2), obsdata(:,3)); %covariance for data
        pops_init = [obsdata(1,2) obsdata(1,3)]';

        [tmod,ymod] = ode45(@(t,y) Lotka_Volterra_Model(t,y,theta),obsdata(:,1), pops_init);
        %ymodthin = interp1(tmod, ymod, time_range_disc);
        norm = log10(mvnpdf(ymod-obspops,0,sampleCov));  
        
        % creates a normal distribution from the difference bt the
        %   predicted model data and the observed data, with mean 0 and
        %   standard dev ysig
        %   evaluating the likelihood function (?) given the data
        llik = sum(norm); % scaled sum of the normal pdf : probability of a distribution
end

