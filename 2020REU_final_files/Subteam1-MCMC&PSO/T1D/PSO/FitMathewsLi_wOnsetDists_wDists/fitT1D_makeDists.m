%% Fitting T1D Model with PSO: Make distributions
%  Author:       Christina Catlett, using code by B. Shtylla as reference
%  Date:         June 10, 2020 (reference code: April 24, 2020)
%  Desc:         Main function to fit parameters to T1D
%                model. Using mean data for Mathews + Li, fits many times, makes parameter historgrams.
%                Requires Global Optimization Toolbox.

% Reset at beginning of run
clear all;
rng default;

% Load dataset
load('mathewsli_prog.mat');

% Select group of interest
y = mathewsliprog.Glucose;
datalen = length(y);

% Number of times to run PSO
nIter = 50;

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
ICRange = .03;
paramRange = .03; % Superceded by vals in eFAST/variance table

% Declare final matrix of vars (run once to get dimensions)
paramsOut = runPSO_makeDists(y, datalen, type, disease, waveon, wICs, subset, bounds, ICRange, paramRange);

for i = 2:nIter
newParams = runPSO_makeDists(y, datalen, type, disease, waveon, wICs, subset, bounds, ICRange, paramRange);
paramsOut = [paramsOut; newParams];
end

writematrix(paramsOut, 'paramsDist2.csv');
for i = 1:60
    histogram(paramsOut(:,i))
    title(int2str(i))
    pause
end