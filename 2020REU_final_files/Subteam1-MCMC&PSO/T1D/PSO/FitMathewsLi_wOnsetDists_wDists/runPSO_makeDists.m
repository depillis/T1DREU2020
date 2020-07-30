function out = runPSO_makeDists(y, datalen, type, disease, waveon, wICs, subset, bounds, ICRange, paramRange)
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
        lb = orig.*(1-paramRange);
        ub = orig.*(1+paramRange);
        
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
        lb = orig.*(1-paramRange);
        ub = orig.*(1+paramRange);
        
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
        lb = orig.*(1-paramRange);
        ub = orig.*(1+paramRange);
        
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
        lb = orig.*(1-paramRange);
        ub = orig.*(1+paramRange);
        
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
        lb = orig.*(1-paramRange);
        ub = orig.*(1+paramRange); 
        
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
objfn = @(params) T1Dss_makeDists(params, y, waveon, type, wICs, datalen);
npars = length(ub);

options = optimoptions('particleswarm', 'SwarmSize', 10, 'MaxIterations', 200, ...
    'UseParallel', true);


d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run particleswarm, plot iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run PSO
[out, ~, ~, ~] = particleswarm(objfn, npars, lb, ub, options);
end