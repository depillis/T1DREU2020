%% T1D ODE Evaluation for 350 day time span
%  Author:       Edited by M. Watanabe, from Shtylla et al 2019
%  Date:         June 2020
%  Desc:         Function to evaluate T1D ODE model for a large time period
%                to evaluate DRAM model. Typically params are from the
%                result of DRAM. Used in validation of model.

function [tmodel, ymodel] = T1Dfun_meanpred(modeltype, parVer, waveOn, params, data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get parameter values
initParams;
% Initial Values for DCs and tDCs
D0 = 0;
tD0 = 0;

%Initial values
%Tspan = data(1,1):350;
Tspan = 0:350;
%IC = [params(end-3) 0 0 0 params(end-2) params(end-1) params(end) D0 tD0 0 0 0];% start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve ODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waveOn == 0 % wave off
% Wave is off
wave = 0;
    if(modeltype==1) % SOLVE ODE FOR WILD TYPE - no wave
    [Tnowave, Ynowave] = ode15s(@(t,y)T1Dode_r(t,y,f1,f2,wave,params,parVer),Tspan,IC); % Solve ODE
    tmodel = Tnowave;
    ymodel = Ynowave;


    elseif(modeltype==0) % SOVLE ODE FOR NOD - no wave
    [Tnnowave, Ynnowave] = ode15s(@(t,y)T1Dode_r(t,y,f1n,f2n,wave,params,parVer),Tspan,IC); % Solve ODE
    tmodel = Tnnowave;
    ymodel = Ynnowave;
    
    end 

elseif waveOn == 1
% Wave is on
wave = 0.75;
    if(modeltype==1) % SOLVE ODE FOR WILD TYPE - wave
    [Twave, Ywave] = ode15s(@(t,y)T1Dode_r(t,y,f1,f2,wave,params,parVer),Tspan,IC); % Solve ODE
    tmodel = Twave;
    ymodel = Ywave;

    elseif(modeltype==0)% SOLVE ODE FOR NOD - wave
    options=odeset('Events',@sickEvents);
     [Tnwave, Ynwave,Te,Ye,ie] =...
         ode15s(@(t,y)T1Dode_r(t,y,f1n,f2n,wave,params,parVer),Tspan,IC,options); % Solve ODE
    tmodel = Tnwave;
    ymodel = Ynwave;
    end
 
else
   msg='Invalid model type: enter 1, 2, 3, or 4';
     error(msg) 
end
end