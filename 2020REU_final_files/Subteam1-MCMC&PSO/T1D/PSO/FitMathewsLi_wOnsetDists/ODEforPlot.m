%% T1D ODE Solver *for plotting convenience*
%  Author:       Altered by C. Catlett, unsure of original author
%  Date:         Altered June 11, 2020
%  Desc:         Solve T1D system given options (waveon, type, wICs) and
%                given parameters (params)

function [Tpre, Ypre, Tpost, Ypost] = ODEforPlot(params, wICs)
% Initial Values for DCs and tDCs
D0 = 0;
tD0 = 0;

%Macrophages
f1ns         = params(29); 
f1s          = params(30); 
f2ns         = params(31);
f2s          = params(32);
scale_factor = params(43); 

f1  = scale_factor*f1s*10^(-5); % ml cell^-1d^-1 %wt basal phagocytosis rate per M
f1n = scale_factor*f1ns*10^(-5);% ml cell^-1d^-1 %nod basal phagocytosis rate per M 
f2  = scale_factor*f2s*10^(-5); % ml cell^-1d^-1 %wt activated phagocytosis rate per Ma
f2n = scale_factor*f2ns*10^(-5);% ml cell^-1d^-1 %nod activated phagocytosis rate per Ma 
wave = 0.75;

if wICs == 1
                IC = params(47:end-2);  
else
                % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
                IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0]; 
end
            options = odeset('Events', @sickEvents);
            [Tpre, Ypre] = ode15s(@(t, y) T1Dsys(t, y, f1n, f2n, wave, params), 0:175, IC, options);
            [Tpost, Ypost] = ode15s(@(t, y) T1Dsys(t, y, f1n, f2n, wave, params), 175:400, Ypre(end, :), options);
end
