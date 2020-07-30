%% Fitting T1D Model with (Mathews & Li Data)
%  Author:       C. Catlett, using code by B. Shtylla as reference
%  Date:         June 10, 2020 (reference code: April 24, 2020)
%  Desc:         Main function to fit parameters to T1D
%                model. Requires Global
%                Optimization Toolbox, various user choices on lines 45-57

% Reset at beginning of run
clear all;
rng default;

% Load dataset
load('mathewsli_prog.mat');

% Select time series, data of interest
t = mathewsliprog.Day;
y = mathewsliprog.Glucose;
datalen = length(y);

% USER CHOICE: NOD vs. wild type mice, type of disease, wave on/off, optimizing ICs or not,
% subset of params to optimize for, bounds choice
%
%       'NOD'   -> NOD
%       'WT'    -> Wild type
% ---------------------------------------------------
%       'prog'  -> Progressive mice
%       'acute' -> Acute mice
% ---------------------------------------------------
%       1       -> Apoptotic wave
%       0       -> No apoptotic wave
% ---------------------------------------------------
%       1       -> Treat ICs as parameters
%       0       -> Constant ICs
% ---------------------------------------------------
%       'a'     -> Optimize all possible params (does not supercede ICs choice)
%       'g'     -> Optimize only params involved in glucose equation
%       's'     -> Optimize only params considered sensitive by eFAST
%       'gs'    -> Optimize both glucose, sensitive parameters
%       'e'     -> Optimize eta parameters
%       'v'     -> Optimize 'volatile' parameters accoring to UKF
% ---------------------------------------------------
%       'eFAST' -> Use eFAST analysis to determine upper, lower bounds
%       'var'   -> Use variance from UKF to determine upper, lower bounds

type = 'NOD';
disease = 'prog';
waveon = 1;
wICs = 1;
subset = 'a';
bounds = 'var';

% USER CHOICE: Define default percent variation for parameters/ICs w/o defined range (IN DECIMAL)
ICRange = .05;
paramRange = .03; % Superceded by vals in eFAST/variance table

% USER CHOICE: Anon function to weight time/shape score to determine goodness of fit
compositeScorefunc = @(t, shape) (t  + log10(shape));

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
    case 'g' % Only glucose params
        % Allow all to vary by +/- paramRange
        lb_temp = orig.*(1-paramRange);
        ub_temp = orig.*(1+paramRange);
        
        % EG0
        lb_temp(2) = lb(2);
        ub_temp(2) = ub(2);
        % R0
        lb_temp(9) = lb(9);
        ub_temp(9) = ub(9);
        % SI
        lb_temp(10) = lb(10);
        ub_temp(10) = ub(10);
        
        lb = lb_temp;
        ub = ub_temp;
    case 's' % Only sensitive params
        % Allow all to vary by +/- paramRange
        lb_temp = orig.*(1-paramRange);
        ub_temp = orig.*(1+paramRange);
        
        % D_ss
        lb_temp(1) = lb(1);
        ub_temp(1) = ub(1);
        % Qpanc
        lb_temp(7) = lb(7);
        ub_temp(7) = ub(7);
        % EG0
        lb_temp(2) = lb(2);
        ub_temp(2) = ub(2);
        % mues_e
        lb_temp(38) = lb(38);
        ub_temp(38) = ub(38);
        % mues_r
        lb_temp(39) = lb(39);
        ub_temp(39) = ub(39);
        % SI
        lb_temp(10) = lb(10);
        ub_temp(10) = ub(10);
        % alpha_B
        lb_temp(14) = lb(14);
        ub_temp(14) = ub(14);
        
        lb = lb_temp;
        ub = ub_temp;        
    case 'gs' % Both glucose, sensitive params
        % Allow all to vary by +/- paramRange
        lb_temp = orig.*(1-paramRange);
        ub_temp = orig.*(1+paramRange);
        
        % EG0
        lb_temp(2) = lb(2);
        ub_temp(2) = ub(2);
        % R0
        lb_temp(9) = lb(9);
        ub_temp(9) = ub(9);
        % SI
        lb_temp(10) = lb(10);
        ub_temp(10) = ub(10);
        % alpha_B
        lb_temp(14) = lb(14);
        ub_temp(14) = ub(14);
        % mues_e
        lb_temp(38) = lb(38);
        ub_temp(38) = ub(38);
        % mues_r
        lb_temp(39) = lb(39);
        ub_temp(39) = ub(39);
        % D_ss
        lb_temp(1) = lb(1);
        ub_temp(1) = ub(1);
        % Qpanc
        lb_temp(7) = lb(7);
        ub_temp(7) = ub(7);
        
        lb = lb_temp;
        ub = ub_temp;
    case 'e' % eta-related parameters
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
    case 'v' % Volatile parameters from UKF
        % Allow all to vary by +/- paramRange
        lb_temp = orig.*(1-paramRange);
        ub_temp = orig.*(1+paramRange); 
        
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

% Create time to onset distribution parameters
[mu, sig, mulb, muub, siglb, sigub] = onset_dist_mathews_christina(disease);

lb(end+1) = mulb;
lb(end+1) = siglb;
ub(end+1) = muub;
ub(end+1) = sigub;

% Objective func -> Minimize sum of squares for y6 (glucose), determine
% shift w/ onset time distribution
objfn = @(params) T1Dss(params, y, waveon, type, wICs, datalen);
%objfn = @(params) T1Dscore(params, t, y, waveon, type, wICs, datalen);
npars = length(ub);

options = optimoptions('particleswarm', 'SwarmSize', 20, 'MaxIterations', 200, ...
    'UseParallel', true, 'Display', 'iter', 'DisplayInterval', 10, 'PlotFcn', @pswplotbestf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run particleswarm, plot iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run PSO
[out, fval, exitflag, output] = particleswarm(objfn, npars, lb, ub, options);

%%%%%%%%%%%%
%% Score fit
%%%%%%%%%%%%

% Solve ODE for best fit (bf) parameters from 'out'
[bfTime_pre, bfData_pre, bfTime_post, bfData_post] = T1Dfun(out, waveon, type, wICs, datalen);

% Determine goodness of fit through comparison of onset time, shape
[tScore, shapeScore] = scoreFit(bfTime_pre, bfTime_post(2:end), bfData_pre(:,6), bfData_post(2:end,6), t , y);
compositeScore = compositeScorefunc(tScore, shapeScore); % Scale scores, weight components

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot best fit params, ODE system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Alert user if fitted params of interest match ub/lb (bounds may need to be changed)
indsofinterest = find(ub ~= ((1+paramRange)/(1-paramRange)).*lb); % Where default range not used, param is of interest
checkEq = @(x)isequal(length(x), length(unique(x)));
if ~checkEq([out(indsofinterest) lb(indsofinterest)]) || ~checkEq([out(indsofinterest) ub(indsofinterest)])
    warning('Best fit parameters exist at boundaries; ub/lb may need to be adjusted');
end    

% Used to visually compare results/bounds
compbounds = [out(indsofinterest); lb(indsofinterest); ub(indsofinterest)];

% Plot glucose data vs. ODE prediction
figure(1); clf
plot([bfTime_pre; bfTime_post(2:end)], y, 'sb')
hold on;
title('PSO best fit for observed glucose');
plot([bfTime_pre; bfTime_post(2:end)], [bfData_pre(:, 6); bfData_post(2:end,6)], '-k');
xlabel('days');
ylabel('glucose (mg/dl)');

% Plot all outcomes of ODE vs. orig model from Shtylla et al.
figure(2); clf
plotAllPSO(out, waveon, type, wICs, bfTime_post);