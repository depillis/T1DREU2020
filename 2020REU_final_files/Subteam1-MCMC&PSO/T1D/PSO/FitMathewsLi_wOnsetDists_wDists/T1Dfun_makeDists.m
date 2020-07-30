%% T1D ODE Solver
%  Author:       Altered by C. Catlett
%  Date:         Altered June 11, 2020
%  Desc:         Solve T1D system given options (waveon, type, wICs) and
%                given parameters (params)

function [Tpre, Ypre, Tpost, Ypost] = T1Dfun_makeDists(params, waveon, type, wICs, datalen)
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

% Initialize inputs for solver
mu = params(end-1);
sig = params(end);
onsettime = 7*ceil(exp(mu+sig^2/2)); % Convert lognormal params to normal (IN DAYS)

% Pre and pos-onset time spans (IN DAYS)
Tspan_pre = onsettime-23:onsettime+1;
Tspan_post = [onsettime+1 makePostTspan(onsettime, Tspan_pre, datalen)];
    
% Solve ODE
switch type
    % Solve for NOD mice
    case 'NOD'
        % Solve w/ apoptotic wave
        if waveon == 1
            if wICs == 1
                IC = params(47:end-2);  
            else
                % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
                IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0]; 
            end
            options = odeset('Events', @sickEvents); 
            [~, IC] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1n, f2n, wave, params), 0:Tspan_pre(1), IC, options);
            IC = IC(end, :);
            options = odeset('Events', @sickEvents);  
            [Tpre, Ypre] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1n, f2n, wave, params), Tspan_pre, IC, options);
            [Tpost, Ypost] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1n, f2n, wave, params), Tspan_post, Ypre(end, :), options);
        % Solve w/o apoptotic wave
        else
            if wICs == 1
               IC = params(47:end-2);  
            else
                % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
                IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0]; 
            end      
            [~, IC] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1n, f2n, 0, params), 0:Tspan_pre(1), IC);
            IC = IC(end, :);
            [Tpre, Ypre] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1n, f2n, 0, params), Tspan_pre, IC);
            [Tpost, Ypost] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1n, f2n, 0, params), Tspan_post, Ypre(end, :));
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
            [~, IC] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1, f2, wave, params), Tspan_pre, IC);
            IC = IC(end, :);
            [Tpre, Ypre] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1, f2, wave, params), Tspan_pre, IC);
            [Tpost, Ypost] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1, f2, wave, params), Tspan_post, Ypre(end, :));
        % Solve w/o apoptotic wave
        else
            if wICs == 1
                IC = params(47:end-2);  
            else
                % Start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin
                IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0];
            end
            [~, IC] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1, f2, 0, params), 0:Tspan_pre(1), IC);
            IC = IC(end, :);
            [Tpre, Ypre] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1, f2, 0, params), Tspan_pre, IC);
            [Tpost, Ypost] = ode15s(@(t, y) T1Dsys_makeDists(t, y, f1, f2, 0, params), Tspan_post, Ypre(end, :));
        end
end
end

% Helper function to build post-onset timespan
function Tspan_post = makePostTspan(onsettime, Tspan_pre, datalen)
Tspan_post(1) = onsettime+7;
while length(Tspan_post) + length(Tspan_pre) < datalen
    Tspan_post = [Tspan_post Tspan_post(end)+7];
end
end