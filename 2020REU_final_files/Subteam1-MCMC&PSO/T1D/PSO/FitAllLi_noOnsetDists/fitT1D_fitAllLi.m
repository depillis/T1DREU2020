%% Fitting T1D Model with (All 9 Li)
%  Author:       C. Catlett, using code by B. Shtylla as reference
%  Date:         June 10, 2020 (reference code: April 24, 2020)
%  Desc:         Main function to fit parameters to T1D
%                model for 'T1Ddata.mat' for 9 Li et al. acute mice. Requires Global
%                Optimization Toolbox, various user choices on lines 30-38

% Reset at beginning of run
clear all;
rng default;

% USER CHOICE: NOD vs. wild type mice, wave on/off, optimizing ICs or not,
% subset of params to optimize for, bounds choice
%
%       'NOD'   -> NOD
%       'WT'    -> Wild type
% ---------------------------------------------------
%       1       -> Apoptotic wave
%       0       -> No apoptotic wave
% ---------------------------------------------------
%       1       -> Treat ICs as parameters
%       0       -> Constant ICs
% ---------------------------------------------------
%       ' '     -> Optimize all possible params (does not supercede ICs choice)
%       'a'     -> Optimize only active params according to UKF
% ---------------------------------------------------
%       'eFAST' -> Use eFAST analysis to determine upper, lower bounds
%       'var'   -> Use variance from UKF to determine upper, lower bounds

type = 'NOD';
waveon = 1;
wICs = 1;
subset = 'a';
bounds = 'var';

% USER CHOICE: Define default percent variation for parameters/ICs w/o defined range (IN DECIMAL)
ICRange = .05;
paramRange = .03; % Superceded by vals in eFAST/variance table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Format particleswarm input, options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default lower, upper bounds for pars
switch bounds
    case 'eFAST' % Use eFAST
        [lb, ub, orig] = makebounds_eFAST(paramRange); % +/- paramRange
    case 'var' % Use variance from UKFs
        [lb, ub, orig] = makebounds_var(paramRange); % +/- paramRange
end

% Select parameter subset if applicable
switch subset
    case 'a' % active parameters
        % Allow all to vary by +/- paramRange
        lb_temp = orig.*(1-paramRange);
        ub_temp = orig.*(1+paramRange);
        
        % eta
        lb_temp(46) = lb(46);
        ub_temp(46) = ub(46);        
       % alpha_eta
        lb_temp(15) = lb(15);
        ub_temp(15) = ub(15);        
       % beta_eta
        lb_temp(21) = lb(21);
        ub_temp(21) = ub(21);  
    
        % GI
        lb_temp(3) = lb(3);
        ub_temp(3) = ub(3);
        % SI
        lb_temp(10) = lb(10);
        ub_temp(10) = ub(10);
        % mues_e
        lb_temp(38) = lb(38);
        ub_temp(38) = ub(38);
        % mues_r
        lb_temp(39) = lb(39);
        ub_temp(39) = ub(39);
        % e1
        lb_temp(27) = lb(27);
        ub_temp(27) = ub(27);
        % e2
        lb_temp(28) = lb(28);
        ub_temp(28) = ub(28);
        % delta_B
        lb_temp(26) = lb(26);
        ub_temp(26) = ub(26);
end    

% Adjust bounds for ICs
if wICs == 1
    % Load file
    ICs = load('ICs.mat');
    [lbICs, ubICs] = makeICbounds(ICs, ICRange); % +/- 10%
    lb = horzcat(lb,lbICs);
    ub = horzcat(ub, ubICs);
end

% Initialize array for params
out_all = [];
% Loop through each acute mouse
for n=[2,3,4,6,7,8,9,10,11]
    
    n_str = int2str(n);
    file = strcat('dat',n_str,'.csv'); %create file name
    X = readtable(file);
    X = X{:,:};
    X = flip(X);
    
    t = X(:,1).*7;
    y = X(:,2);
    
% Objective func -> Minimize sum of squares for y6 (glucose), determine
% shift w/ onset time distribution
objfn = @(params) T1Dss_noDist(params, y, waveon, type, wICs, t);
npars = length(ub);

options = optimoptions('particleswarm', 'SwarmSize', 20, 'MaxIterations', 200, ...
    'UseParallel', true, 'Display', 'iter', 'DisplayInterval', 10, 'PlotFcn', @pswplotbestf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run particleswarm, plot iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run PSO
[out, fval, exitflag, output] = particleswarm(objfn, npars, lb, ub, options);

%%%%%%%%%%%%
%% Store fit
%%%%%%%%%%%%

% Add optimal params to array
out_all = [out_all; out]
end