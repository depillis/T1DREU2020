function initialConditions = Get_Initial_Conditions(T)
%This function gets the conditions at the specified time point

wave = 0;
T1D_ODE_Parameters;
truepar = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8); 
                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21]; %True parameter values - taken from original parameters file
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0; truepar];
%T = 280;
tspan = [0 T];
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, solutions] = ode15s(@(t, y) ODE_allparams(t, y, f1n, f2n, wave, truepar), tspan, InitialState,options); %Use ODE solver
initialConditions = solutions(end,:);
end