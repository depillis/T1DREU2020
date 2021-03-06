%% Type 1 Diabetes Model - Delayed Rejection Adaptive Metropolis MCMC
%  Author:       Edited to fit T1D by M. Watanabe & C. Catlett
%  Date:         Summer 2020
%  Desc:         Use delayed rejection adaptive Metropolis (DRAM) MCMC 
%                methods to
%                parameterize a 12 equation type 1 diabetes ODE system. 
%                This code produces posterior distributions for chosen 
%                parameters and the predicted glucose model. 
%                This code utilizes the mcmcstat library (M.J. Laine).
%
%  Sources: 
%       Algorithm adapted from code by Marko J. Laine
%       (https://mjlaine.github.io/mcmcstat/ex/algaeex.html). 
%
%       Glucose data from Mathews et al (2015) and Li et al (2009).
%
%       T1D ODE model from Shtylla et al (2019).
%
%  Note: This script uses the mcmcstat library, but the main functions,
%  mcmcrun and mcmcpred, have been altered to accomodate user preferences.
%  These functions are now mcmcrun_custom and
%  mcmcpred_custom. The original functions still exist in the
%  mcmcstat folder.
%
% Directions:   This script has user-specified options that must be
%               designated before run. To run click on the "Run" button 
%               in the Editor tab of MATLAB. This script takes ~60mins to
%               parameterize a subset of 10 parameters.
%
%% SET USER PREFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clear model data params options
% add mcmcstat folder to path
addpath('mcmcstat');
addpath('helperFuncs');
addpath('validation&comparison_funcs');

% USER CHOICE: NOD vs. wild type mice, type of disease, wave on/off,
% optimizing ICs or not, subset of params to optimize for

% DATASET
%       0       -> Mathews et al.
%       1       -> Li et al. averaged is DEFAULT
%               -> options exist for individual mouse data sets see
%                  getData.m
%       2       -> Simulated data + parameter set 'base/1' ONLY
        dataset = 1;
% --------------------------------------------------
% PRIOR FUNCTION
%       'mse'   -> no prior function do a mean squared error instead
%                  (requires #data pts > # parameters)
%       'prior' -> predetermined (log)prior
        priorinfo = 'mse';
% --------------------------------------------------
% MODELTYPE
%
%       0       -> NOD
%       1       -> Wild type
        modeltype = 0;
% ---------------------------------------------------
% DISEASE
%       0       -> Progressive mice
%       1       -> Acute mice
        disease = 1;
% ---------------------------------------------------
% WAVE 
%       0       -> Apoptotic wave OFF
%       1       -> Apoptotic wave ON
        wave = 1;
% ---------------------------------------------------
% PARAMETER sets and version for mcmcrun: param_subset and parVer
%    'base'/1   -> Paremeterize sensitive (UKF), eta_vary, and
%                  initial conditions (IC)
%    'base+/1'  -> Same set as 'base', but with mu and sigma taken from
%                  previous chainstats run. Just to test logprior function
%    'gluc'/2   -> Paremeterize glucose, eta_vary, eFAST sensitive, time,
%                  IC parameters
%    'all'/3    -> Parameterize params considered sensitive by eFAST and
%                  UKF analysis, time, IC parameters
%    'all+/3'   -> Same set as 'all', but with mu and sigma taken from
%                  post-UKF distributions, initial condition priors taken
%                  from previous DRAM runs (chainstats)
      param_subset = 'all';
      parVer = 3;

% =========================================================================
% SAVING WORKSPACE POST RUN
% set file name to save workspace (at end of algorithm)

% specify unique name head for saving - ex. date/project number
head = '';

filename = filename(head, modeltype, wave, disease, dataset);

%% LOAD, FORMAT, AND PLOT OBSERVED DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getData function formats data and plots
data = getData(dataset, disease);
%% DEFINE INITIAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All parameters are constrained to be positive.
% {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}

% initial conditions for parameters
[lb, ub, init] = makeparambounds; % 41 params
[lbIC, ubIC, initIC] = makeICbounds(0.1); % Initial conditions
[parmu, parsig2] = makepar_musig; % mean and standard deviation

% SELECT PARAMETER SUBSET
switch param_subset
    case 'base' % sensitive, eta_vary, initial conditions parameters
        params = {
        % sensitive parameters (D. Shenker's table)
        {'e1', init(5), lb(5), ub(5)}
        {'e2', init(6), lb(6), ub(6)}
        {'delta_B', init(10), lb(10), ub(10)}
        {'SI', init(15), lb(15), ub(15)}
        {'GI', init(18), lb(18), ub(18)}
        {'mues_r', init(34), lb(34), ub(34)}
        {'mues_e', init(35), lb(35), ub(35)}

        % eta_vary parameters
        {'alpha_eta', init(40), lb(40), ub(40)}
        {'beta_eta', init(41), lb(41), ub(41)}
        {'eta', 0.0100, 0.0100, 0.0300} % taken from An's eFAST parameters
        };
    case 'base+'
        params = {
        % sensitive parameters (D. Shenker's table)
        {'e1', init(5), lb(5), ub(5), parmu(5), parsig2(5)}
        {'e2', init(6), lb(6), ub(6), parmu(6), parsig2(6)}
        {'delta_B', init(10), lb(10), ub(10), parmu(10), parsig2(10)}
        {'SI', init(15), lb(15), ub(15), parmu(15), parsig2(15)}
        {'GI', init(18), lb(18), ub(18), parmu(18), parsig2(18)}
        {'mues_r', init(34), lb(34), ub(34), parmu(34), parsig2(34)}
        {'mues_e', init(35), lb(35), ub(35), parmu(35), parsig2(34)}

        % eta_vary parameters
        {'alpha_eta', init(40), lb(40), ub(40), parmu(39), parsig2(39)}
        {'beta_eta', init(41), lb(41), ub(41), parmu(40), parsig2(40)}
        {'eta', 0.0100, 0.0100, 0.0300, parmu(41), parsig2(41)} % taken from An's eFAST parameters
        };
    case 'gluc' % glucose, eta_vary, An's eFAST sens., time parameters
        params = {
        % Glucose parameters
        {'R0', init(13), lb(13), ub(13)}
        {'EG0', init(14), lb(14), ub(14)}
        {'SI', init(15), lb(15), ub(15)}

 
        % eta_vary parmeters
        {'alpha_eta', init(40), lb(40), ub(40)}
        {'beta_eta', init(41), lb(41), ub(41)}

        % An's sensitive parameters (eFAST)
        {'D_ss', init(24), lb(24), ub(24)}
        {'alpha_B', init(9), lb(9), ub(9)}
        {'Qpanc', init(19), lb(19), ub(19)}
        {'mues_r', init(34), lb(34), ub(34)}
        {'mues_e', init(35), lb(35), ub(35)}

        % time parameters
        {'onset_mu', mu, distCI(1,1), distCI(2,1)}
        {'onset_sigma', sigma, distCI(1,2), distCI(2,2)}
        };
    case 'all' % all 41, time parameters
        params = {
        {'J', init(1), lb(1), ub(1)}
        {'k', init(2), lb(2), ub(2)}
        {'b', init(3), lb(3), ub(3)}
        {'c', init(4), lb(4), ub(4)}
        {'e1', init(5), lb(5), ub(5)}
        {'e2', init(6), lb(6), ub(6)}
        {'fD', init(7), lb(7), ub(7)}
        {'ftD', init(8), lb(8), ub(8)}
        {'alpha_B', init(9), lb(9), ub(9)}
        {'delta_B', init(10), lb(10), ub(10)}
        {'B_conv', init(11), lb(11), ub(11)}
        {'Ghb', init(12), lb(12), ub(12)}
        {'R0', init(13), lb(13), ub(13)}
        {'EG0', init(14), lb(14), ub(14)}
        {'SI', init(15), lb(15), ub(15)}
        {'sigmaI', init(16), lb(16), ub(16)}
        {'deltaI', init(17), lb(17), ub(17)}
        {'GI', init(18), lb(18), ub(18)}
        {'Qpanc', init(19), lb(19), ub(19)}
        {'Qblood', init(20), lb(20), ub(20)}
        {'Qspleen', init(21), lb(21), ub(21)}
        {'mu_PB', init(22), lb(22), ub(22)}
        {'mu_BP', init(23), lb(23), ub(23)}
        {'D_ss', init(24), lb(24), ub(24)}
        {'bDEday', init(25), lb(25), ub(25)}
        {'bIRday', init(26), lb(26), ub(26)}
        {'aEaday', init(27), lb(27), ub(27)}
        {'T_naive', init(28), lb(28), ub(28)}
        {'bpday', init(29), lb(29), ub(29)}
        {'ramday', init(30), lb(30), ub(30)}
        {'baEday', init(31), lb(31), ub(31)}
        {'baRday', init(32), lb(32), ub(32)}
        {'aEmday', init(33), lb(33), ub(33)}
        {'mues_r', init(34), lb(34), ub(34)}
        {'mues_e', init(35), lb(35), ub(35)}
        {'thetaD', init(36), lb(36), ub(36)}
        {'d', init(37), lb(37), ub(37)}
        {'sE', init(38), lb(38), ub(38)}
        {'sR', init(39), lb(39), ub(39)}
        {'alpha_eta', init(40), lb(40), ub(40)}
        {'beta_eta', init(41), lb(41), ub(41)}

        % Vary non-zero initial conditions (IC)
        % IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0];
        {'resting_macrophages', initIC(1), lbIC(1), ubIC(1)}
        {'beta_cells', initIC(2), lbIC(2), ubIC(2)}
        {'glucose', initIC(3), lbIC(3), ubIC(3)}
        {'insulin', initIC(4), lbIC(4), ubIC(4)}
        
%         % time parameters
%         {'onset_mu', mu, distCI(1,1), distCI(2,1)}
%         {'onset_sigma', sigma, distCI(1,2), distCI(2,2)}
        };
    case 'all+'
        params = {
        {'J', init(1), lb(1), ub(1), parmu(1), parsig2(1)}
        {'k', init(2), lb(2), ub(2), parmu(2), parsig2(2)}
        {'b', init(3), lb(3), ub(3), parmu(3), parsig2(3)}
        {'c', init(4), lb(4), ub(4), parmu(4), parsig2(4)}
        {'e1', init(5), lb(5), ub(5), parmu(5), parsig2(5)}
        {'e2', init(6), lb(6), ub(6), parmu(6), parsig2(6)}
        {'fD', init(7), lb(7), ub(7), parmu(7), parsig2(7)}
        {'ftD', init(8), lb(8), ub(8), parmu(8), parsig2(8)}
        {'alpha_B', init(9), lb(9), ub(9), parmu(9), parsig2(9)}
        {'delta_B', init(10), lb(10), ub(10), parmu(10), parsig2(10)}
        {'B_conv', init(11), lb(11), ub(11), parmu(11), parsig2(11)}
        {'Ghb', init(12), lb(12), ub(12), parmu(12), parsig2(12)}
        {'R0', init(13), lb(13), ub(13), parmu(13), parsig2(13)}
        {'EG0', init(14), lb(14), ub(14), parmu(14), parsig2(14)}
        {'SI', init(15), lb(15), ub(15), parmu(15), parsig2(15)}
        {'sigmaI', init(16), lb(16), ub(16), parmu(16), parsig2(16)}
        {'deltaI', init(17), lb(17), ub(17), parmu(17), parsig2(17)}
        {'GI', init(18), lb(18), ub(18), parmu(18), parsig2(18)}
        {'Qpanc', init(19), lb(19), ub(19), parmu(19), parsig2(19)}
        {'Qblood', init(20), lb(20), ub(20), parmu(20), parsig2(20)}
        {'Qspleen', init(21), lb(21), ub(21), parmu(21), parsig2(21)}
        {'mu_PB', init(22), lb(22), ub(22), parmu(22), parsig2(22)}
        {'mu_BP', init(23), lb(23), ub(23), parmu(23), parsig2(23)}
        {'D_ss', init(24), lb(24), ub(24), parmu(24), parsig2(24)}
        {'bDEday', init(25), lb(25), ub(25), parmu(25), parsig2(25)}
        {'bIRday', init(26), lb(26), ub(26), parmu(26), parsig2(26)}
        {'aEaday', init(27), lb(27), ub(27), parmu(27), parsig2(27)}
        {'T_naive', init(28), lb(28), ub(28), parmu(28), parsig2(28)}
        {'bpday', init(29), lb(29), ub(29), parmu(29), parsig2(29)}
        {'ramday', init(30), lb(30), ub(30), parmu(30), parsig2(30)}
        {'baEday', init(31), lb(31), ub(31), parmu(31), parsig2(31)}
        {'baRday', init(32), lb(32), ub(32), parmu(32), parsig2(32)}
        {'aEmday', init(33), lb(33), ub(33), parmu(33), parsig2(33)}
        {'mues_r', init(34), lb(34), ub(34), parmu(34), parsig2(34)}
        {'mues_e', init(35), lb(35), ub(35), parmu(35), parsig2(35)}
        {'thetaD', init(36), lb(36), ub(36), parmu(36), parsig2(36)}
        {'d', init(37), lb(37), ub(37), parmu(37), parsig2(37)}
        {'sE', init(38), lb(38), ub(38), parmu(38), parsig2(38)}
        {'sR', init(39), lb(39), ub(39), parmu(39), parsig2(39)}
        {'alpha_eta', init(40), lb(40), ub(40), parmu(40), parsig2(40)}
        {'beta_eta', init(41), lb(41), ub(41), parmu(41), parsig2(41)}
        {'eta', 0.0100, 0.0100, 0.0300, 0.020334, 0.0065641^2}  
end

%% SUM OF SQAURES FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.ssfun = @(theta, data, modeltype, parVer, wave) T1Dss_r(theta, data,...
    modeltype, parVer, wave);

%% DEFINE INITIAL VARIANCE FOR PRIOR FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We assume having at least some prior information on the
% repeatability of the observation and assign rather non informational
% prior for the residual variances of the observed states. The default
% prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (see for example Gelman et al.).

% default to uniform prior if #data pts < #parameters
% for some reason these are needed in order to save the s2chain (below)
model.S20  = [1]; 
model.N0   = [1]; 

switch priorinfo
    case 'mse'
        if length(params) <= length(data)
        % T1D_fitInitialParams runs a quick fminsearch to get initial parameter and
        % s2 guesses for DRAM. Will only run if #data pts > #parameters
        [mse, ~] = T1D_fitInitialParams(params, data, modeltype, parVer, wave);
        model.sigma2 = mse;
        end
    case 'prior'
    % Establish unique prior function
    % To use model.priorfun - prior function, need mu and sigma for each of the
    % parameters
    % hard-coded as NORMAL pdf -> change accordingly in logprior.m
    model.priorfun = @(params, paramsmu, paramssigma) logprior(params,...
        paramsmu, paramssigma); % produces logvalue
   
end
%% BURN-IN ITERATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First generate an initial chain.

options.nsimu = 1000;
options.method = 'dram';
res = [];
[results, chain, s2chain]= mcmcrun_custom(model,data,params,options,...
                           modeltype,parVer,wave,res);

%% RUN & SAMPLE MCMC SAMPLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then re-run starting from the results of the previous run

options.nsimu = 1e4;
options.method = 'dram';
[results, chain, s2chain] = mcmcrun_custom(model,data,params,options,...
                            modeltype,parVer,wave,results);

%% PLOT CHAIN RESULTS AND DENSITIES - CONVERGENCE CHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chain plots should reveal that the chain has converged and we can
% use the results for estimation and predictive inference.
%cannot plot more than 10 pairs - this has 53
if length(params) <= 10
figure(2); clf
mcmcplot(chain,[],results,'pairs');
end
% plot densities
figure(3); clf
mcmcplot(chain,[],results,'denspanel',2);


% plot chain
figure(4); clf
for i = 1:length(params)
    hold on
    subplot(4,3,i)
    plot(chain(:,i))
    set(gca, 'FontSize', 15)
    %xlabel('iteration of chain')
end
han=axes(figure(4),'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('iteration of chain')
ylabel('parameter value')
set(gca, 'FontSize', 20)
%% COMPUTE STATISTICS FOR RESULTANT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function |chainstats| calculates mean ans std from the chain and
% estimates the Monte Carlo error of the estimates. Number |tau| is
% the integrated autocorrelation time and |geweke| is a simple test
% for a null hypothesis that the chain has converged.

% This saves the chainstats data but it is only written for the parameter
% set "base/1' and "base+/1"
chainstats(chain,results)
cs = chainstats(chain,results);
pName = {'e1'; 'e2'; 'delta_B'; 'SI'; 'GI'; 'mu_R'; 'mu_E'; 'alpha_eta';...
         'beta_eta'; 'eta'};
mean = cs(:,1);
std = cs(:,2);
MC_err = cs(:,3);
tau = cs(:,4);
geweke = cs(:,5);
T = table(pName, mean, std, MC_err, tau, geweke);
writetable(T, strcat(filename, 'chainstats.txt'))

%% PLOT PREDICTION VS. OBSERVED DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to use the |mcmcpred| function we need
% function |modelfun| with input arguments given as
% |modelfun(xdata,theta)|. We construct this as an anonymous function.

% Evaluates ODE from time 0 to 350 days
modelpred = @(modeltype,parVer,waveOn,theta,data)...
    T1Dfun_meanpred(modeltype,parVer,waveOn,theta,data);

% We sample 500 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.
nsample = 500;

% Get 350 data points
meanParams = results.theta;
[t, y] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);
dataPred = [t,y];

out = mcmcpred_custom(results,chain,s2chain,data,modelpred,nsample,...
               modeltype,parVer,wave);


%% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES

% Plot model predictions
figure(5); clf
mcmcpredplot(out); 
hold on
scatter(data(:,1), data(:,2), 25, 'k', 'filled', 'o'); % our data is only from Glucose measurements so we only plot data over 1 of the densities

set(gca, 'FontSize', 15)
xlabel('day from birth');
ylabel('glucose (mg/dl)');

if disease == 0
    title('Glucose model prediction: progressive', 'Fontsize', 15);

else
    title('Glucose model prediction: acute','Fontsize', 15);
end

% Plot prediction of glucose with mean DRAM parameters
meanParams = results.theta;
[t, y] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);
% save predicted data
writematrix([t,y(:,6)], strcat(filename, '_predmodData.csv'));

% plot
figure(6); clf
plot(t,y(:,6), 'LineWidth',2);
hold on
scatter(data(:,1), data(:,2), 25, 'k', 'filled','o');
set(gca, 'FontSize', 15)
xlabel('day from birth');
ylabel('glucose (mg/dl)');

if disease == 0
    legend('Model prediction', 'Li data, progressive', 'Location', 'southeast');
    %title('Glucose model prediction w/mean parameter values: progressive','Fontsize', 15);

else
    legend('Model prediction', 'Li data, acute', 'Location', 'southeast');
    %title('Glucose model prediction w/mean parameter values: acute','Fontsize', 15);
    
end

% Plot prediction of all 12 variables with mean DRAM parameters

% Titles and plotting data
plotTitles{1}='Resting Macrophages (M)';
plotTitles{2}='Activated Macrophages (Ma)';
plotTitles{3}='Apoptotic beta cells (Ba)';
plotTitles{4}='Necrotic beta cells (Bn)';
plotTitles{5}='Healthy beta cells (B)';
plotTitles{6}='Glucose (G)';
plotTitles{7}='Insulin (I)';
plotTitles{8}='Immunogenic DCs (D)';
plotTitles{9}='Tolerogenic DCs (tD)';
plotTitles{10}='Effector T cells (E)';
plotTitles{11}='Regulatory T cells (R)';
plotTitles{12}='Memory T cells (Em)';


figure(7); clf
for i = 1:6
    hold on
    subplot(3,2,i)
    plot(t, y(:,i), 'LineWidth',2)
    set(gca, 'FontSize', 15)
    xlabel('day');
    title(plotTitles{i});
end

figure(8); clf
for i = 7:12
    hold on
    subplot(3,2,i-6)
    plot(t, y(:,i), 'LineWidth',2)
    set(gca, 'FontSize', 15)
    xlabel('day');
    title(plotTitles{i});
end

%% VALIDATION OF MODEL
% VISUAL VALIDATION USING LI ET AL DATA
save(filename)
workspace = strcat(filename,'.mat');
vis_validation_lietal(workspace, modeltype, disease, parVer, wave);

% Assign model a score based on standard error of the estimate
[tmod, ymod] = T1Dfun_r(modeltype,parVer,wave,meanParams,data);
score = error_mse(data(:,2), ymod);
writematrix([tmod,ymod], strcat(filename,'_errorestData.csv'));
%% SAVE VARIABLES
% save workspace
save(filename)


