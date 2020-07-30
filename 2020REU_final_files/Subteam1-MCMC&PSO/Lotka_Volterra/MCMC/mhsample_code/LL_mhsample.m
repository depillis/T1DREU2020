%% Loglikelihood function for mhsample
% Authors: Maya Watanabe & Christina Catlett 
% June 2020
% Based on code by Reuel Smith 2015-2017
% v. MATLAB R2015b through 2017a
% Description: function to compute the loglikelihood function for the
% Lotka-Volterra model

function llik = LL_mhsample(theta,obsdata)
%   Return log10(likelihood) for parameters theta over all of obsdata
        % Solve ode model
        % Re-index Lotka-Volterra data
        % From PSOPC_LVModel
        %idx = 60:91;
        %mytime = (1:length(idx))';
        %mydata(:,1) = obsdata(idx,2);
        %mydata(:,2) = obsdata(idx,3);
   
        %init=[mydata(1,1) mydata(1,2)];
        %y0(1) = mydata(1,1); y0(2) = mydata(1,2);
        time=0:length(obsdata)-1;
        init=[obsdata(1,2) obsdata(1,3)];
        mydata=[obsdata(:,2) obsdata(:,3)];
        % ode solver
        [t,sol] = ode45(@(t,y) Lotka_Volterra_Model(t,y,theta),time, init);

            
        % creates a normal distribution from the difference bt the
        %   predicted model data and the observed data, with mean 0 and
        %   standard dev ysig
        %   evaluating the likelihood function (?) given the data
        sampleCov=cov(sol-mydata);
        norm = log10(mvnpdf(sol-mydata,0,sampleCov));  

        llik = sum(norm); % scaled sum of the normal pdf : probability of a distribution
end