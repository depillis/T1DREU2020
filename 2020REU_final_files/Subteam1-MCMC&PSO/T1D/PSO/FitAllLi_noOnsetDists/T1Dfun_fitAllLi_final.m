%% T1D ODE Solver (at avg times)
%  Author:       Altered by C. Catlett, unsure of original author
%  Date:         Altered July 11, 2020
%  Desc:         Provides final eval of T1D system given options (waveon, type, wICs) and
%                given parameters (params) at points in avg
%                data *used for MSE vals*

function answer = T1Dfun_fitAllLi_final(params, waveon, type, wICs)
% Initial Values for DCs and tDCs
D0 = 0;
tD0 = 0;

load('avglietal.mat')
time = avglietal2.time.*7

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
    
% Solve ODE
switch type
    % Solve for NOD mice
    case 'NOD'
        % Solve w/ apoptotic wave
        if waveon == 1
            if wICs == 1
                IC = params(47:end);  
            else
                % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
                IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0]; 
            end
            options = odeset('Events', @sickEvents);  
            % Run model
            modPred = ode15s(@(t, y) T1Dsys_noDist(t, y, f1n, f2n, wave, params), [0 350], IC, options);
            % Eval at avg times
            answer = deval(modPred, time);
        % Solve w/o apoptotic wave
        else
            if wICs == 1
               IC = params(47:end-2);  
            else
                % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
                IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0]; 
            end 
            % Run model
            modPred = ode15s(@(t, y) T1Dsys_noDist(t, y, f1n, f2n, 0, params), time, IC);
            % Eval at avg times
            answer = deval(modPred, time);
        end
    % Solve for wild type mice
    otherwise
        % Solve w/ apoptotic wave
        if waveon == 1
            if wICs == 1
                IC = params(47:end-2);  
            else
                % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
                IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0]; 
            end
            % Run model
            modPred = ode15s(@(t, y) T1Dsys_noDist(t, y, f1, f2, wave, params), time, IC);
            % Eval at avg times
            answer = deval(modPred, time);
        % Solve w/o apoptotic wave
        else
            if wICs == 1
                IC = params(47:end-2);  
            else
                % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
                IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0];
            end
            % Run model
            modPred = ode15s(@(t, y) T1Dsys_noDist(t, y, f1, f2, 0, params), time, IC);
            % Eval at avg times
            answer = deval(modPred, time);
        end
end
end