%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Appendix C of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB

% Metropolis (MCMC) algorithm
%   The first # The first argument of func must be a real vector of parameters, the initial values
% of which are provided by the real vector thetaInit.
% func() returns a two-element vector, the logPrior and logLike (log base 10), the sum
% of which is taken to be the log of the density function (i.e. unnormalized posterior).
% The MCMC sampling PDF is the multivariate Gaussian with fixed covariance, sampleCov.
% A total of Nburnin+Nsamp samples are drawn, of which the last Nsamp are kept.
% As the sampling PDF is symmetric, the Hasting factor cancels,
% leaving the basic Metropolis algorithm.
% Diagnostics are printed very verbose^th sample: sample number, acceptance rate so far.
% ... is used to pass data, prior parameters etc. to func()
% Return a Nsamp * (2+Ntheta) matrix (no names), where the columns are
% 1:  log10 prior PDF
% 2:  log10 likelihood
% 3+: Ntheta parameters
% (The order of the parameters in thetaInit and sampleCov must match, of course.)
function [funcSamp] = metroMCMC(thetaInit, Nburnin, Nsamp, verbose, sampleCov, obsdata)
    
    Ntheta=length(thetaInit); % number of parameters (in this case 3)
    thetaCur=thetaInit; % setting up a current theta (parameters)
    funcCur=logpostLVmodel(thetaCur,obsdata); % in this case calls logpostlinearmodel: 
    %   returns numerator of Bayes Theorem aka the posterior (unnormalized)
    
    funcSamp=zeros(Nsamp,2+Ntheta); % initialize an empty matrix to store posteriors later on
    % why are there so many columns?
    nAccept=0; % # of accepted parameters
    acceptRate=0; % percentage of accepted parameters
    nb=Nburnin+Nsamp; % number of times to sample
    % burnin is 0 in the tutorial, why is that?
    
    % For each sample:
    for n=1:nb
        % Metropolis alg. No Hastings factor for symmetric sampling
        %   distributions
        a0Prop = mvnrnd(thetaCur,sampleCov); % giving us a new proposed paramters (theta) using proposal distribution
        aProp = mvnrnd(thetaCur,sampleCov);
        sigProp = mvnrnd(thetaCur,sampleCov);
        sigProp2 = mvnrnd(thetaCur,sampleCov);
        %thetaProp = [a0Prop(1) aProp(2) sigProp(3), sigProp2(4)];
        thetaProp=[a0Prop(1) a0Prop(2) a0Prop(3) a0Prop(4)];
        %   using prior (covariance matrix)
        funcProp=logpostLVmodel(thetaProp, obsdata); % gives a new posterior distribution using the new proposed theta
        %   this is a scalar probability
        logMR=sum(funcProp)-sum(funcCur); % log10 of the Metropolis ratio
        % the dif betw the posterior probabilities of our new parameters and our old
        %   parameters
        
        % deciding whether accept or reject a candidate theta (proposed
        %   parameters: thetaProp)
        if logMR>=0 || logMR>log10(unifrnd(0,1,1)) % if posteriors probabilities of new is > old
            %new more likely than old
        % logMR>= 0: accept if the thetaProp is more likely/posterior probability is
        %   greater than the thetaCur (old theta)
        % OR 
        % logMR<0: draw from a uniform probability and accept with
        %   probability logMR
            thetaCur=thetaProp; % proposed theta becomes the new theta
            funcCur=funcProp; % proposed posterior becomes the new posterior
            nAccept=nAccept+1; % incrementing the number of accepted posteriors
            acceptRate=nAccept/n; % incrementing the acceptance rate
        end
        
        % filling in the empty matrix
        if n>Nburnin
           funcSamp(n-Nburnin, 1:2)=funcCur; % fills in the empty matrix with the posteriors from each sample 
           funcSamp(n-Nburnin, 3:(2+Ntheta))=thetaCur; % filling in the 'found' parameters for each sample
        end
        
        % diagnostics 
        %if(isfinite(verbose) && (mod(n,verbose)==0 || n==Nburnin+Nsamp))
            %fprintf(n,' of ' ,Nburnin, ' + ', Nsamp,' ',acceptRate);
        %end
        
    end %end for loop
end

