%% T1D ODE Model Solver
%  Author:       Edited by M. Watanabe, from Shtylla et al. 2019
%  Date:         June 2020
%  Desc:         Function to evaluate T1D ODE model in DRAM algorithm

function [tmodel, ymodel] = T1Dfun_r(modeltype, parVer, waveOn, params, data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get parameter values
initParams;
% Initial Values for DCs and tDCs
D0 = 0;
tD0 = 0;

%Initial values
%Tspan = data(:,1);
Tspan_pre = 0:14:data(1,1);
Tspan_post = data(:,1);
%IC = [params(end-3) 0 0 0 params(end-2) params(end-1) params(end) D0 tD0 0 0 0]; % initial conditions
IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve ODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if waveOn == 0 % wave off
wave = 0;
    if(modeltype==1) % SOLVE ODE FOR WILD TYPE - no wave
    [~, IC] = ode15s(@(t,y)T1Dode_r(t,y,f1,f2,wave,params,parVer),Tspan_pre,IC);% Solve ODE
    IC = IC(end,:);
    [Tpost, Ypost] = ode15s(@(t,y)T1Dode_r(t,y,f1,f2,wave,params,parVer),Tspan_post,IC);
    tmodel = Tpost;
    ymodel = Ypost(:,6);


    elseif(modeltype==0) % SOVLE ODE FOR NOD - no wave
    [~, IC] = ode15s(@(t,y)T1Dode_r(t,y,f1n,f2n,wave,params,parVer),Tspan_pre,IC); % Solve ODE
    IC = IC(end,:);
    [Tpost, Ypost] = ode15s(@(t,y)T1Dode_r(t,y,f1n,f2n,wave,params,parVer),Tspan_post,IC);
    tmodel = Tpost;
    ymodel = Ypost(:,6);
    
    end 

elseif waveOn == 1 % wave on
wave = 0.75;
    if (modeltype==1) % SOLVE ODE FOR WILD TYPE - wave
    [~, IC] = ode15s(@(t,y)T1Dode_r(t,y,f1,f2,wave,params,parVer),Tspan_pre,IC); % Solve ODE
    IC = IC(end,:);
    [Tpost, Ypost] = ode15s(@(t, y) T1Dsys(t, y, f1, f2, 0, params), Tspan_post, IC);
    tmodel = Tpost;
    ymodel = Ypost(:,6);

    elseif(modeltype==0) % SOLVE ODE FOR NOD - wave
    options=odeset('Events',@sickEvents);
%      [Tnwave, Ynwave,Te,Ye,ie] =...
%          ode15s(@(t,y)T1Dode_r(t,y,f1n,f2n,wave,params,parVer),Tspan,IC,options); % Solve ODE
%     tmodel = Tnwave;
%     ymodel = Ynwave(:,6);
     [~, IC] =...
         ode15s(@(t,y)T1Dode_r(t,y,f1n,f2n,wave,params,parVer),Tspan_pre,IC,options); % Solve ODE
     IC = IC(end,:);
     [Tpost, Ypost] = ode15s(@(t,y)T1Dode_r(t,y,f1n,f2n,wave,params,parVer),Tspan_post,IC,options);
    tmodel = Tpost;
    ymodel = Ypost(:,6);
    end

else
   msg='Invalid model type: enter 1, 2, 3, or 4';
     error(msg) 
end
end