%% Plotting final PSO fits (Mathews & Li Data)
%  Author:       Christina Catlett
%  Date:         June 19, 2020
%  Desc:         Plot all equations of best fit parameters vs. original
%                system provided in Shtylla et al.; Used as a visual check
%                of the model's reasonableness

function plotAllPSO(params, waveon, type, wICs, bfTime_post)
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

% Select initial conditions
if wICs == 1
    IC = params(47:end-2);  
else
    % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
    IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0]; 
end

% Solve ODE
switch type
    % Solve for NOD mice
    case 'NOD'
        % Solve w/ apoptotic wave
        if waveon == 1
            options = odeset('Events', @sickEvents); 
            [T, Y] = ode15s(@(t, y) T1Dsys(t, y, f1n, f2n, wave, params), 0:bfTime_post(end), IC, options);
        % Solve w/o apoptotic wave
        else
            [T, Y] = ode15s(@(t, y) T1Dsys(t, y, f1n, f2n, 0, params), 0:bfTime_post(end), IC);
        end
    % Solve for wild type mice
    otherwise
        % Solve w/ apoptotic wave
        if waveon == 1
            [T, Y] = ode15s(@(t, y) T1Dsys(t, y, f1, f2, wave, params), 0:bfTime_post(end), IC);
        % Solve w/o apoptotic wave
        else
            [T, Y] = ode15s(@(t, y) T1Dsys(t, y, f1, f2, 0, params), 0:bfTime_post(end), IC);
        end
end

% Plotting 
tiledlayout(6,2);
for i = 1:12
    nexttile
    plot(T, Y(:,i));
    title(string(i))
end
end