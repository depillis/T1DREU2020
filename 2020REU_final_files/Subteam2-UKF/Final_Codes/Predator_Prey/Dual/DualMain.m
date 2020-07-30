%MAIN FILE FOR THE DUAL UKF ON THE PREDATOR-PREY SYSTEM
%THE UKF FUNCTIONS ARE ADAPTED FROM CHOW-FERRER AND THEN HAVE BEEN FULLY
%IMPLEMENTED BY DANIEL SHENKER AND RACHEL WANDER




clear all;
clc;
echo off;
%PARAMETERS
%ALPHA - PREY GROWTH RATE
%BETA - PREDATOR GROWTH RATE
%GAMMA - PREDATION RATE
%DELTA - REPRODUCTION RATE OF PREDATORS PER PREY EATEN

T = 91;  %Number of data points
Ny = 2;  %Number of observables
Nx = 2;  %Number of states
Nstates=2; %number states
Npar=4;   %number parameters
Nsubj=1;
truepar = [0.62526 0.6607 0.1896 0.0468];
truepar_generating = [0.62526 0.6607 0.1896 0.0468 0.0001 0.0001 20 20];

%SET UKF PARAMETERS

alpha = 10e-3; % UKF : point scaling parameter
beta = 2; %Prior knowledge of distribution (for gaussian use 2)
kappa = 0;    %second scaling parameter, set to 0 or 3 - L 


rand('state',sum(100*clock));   
randn('state',sum(100*clock));   
HLData = load('HaresLynxData.mat');  %Load dataset
rawData = HLData.Lotka_Volterra_Data; %Draw out the data
  
x(1:2,:) = rawData(:, 2:3)'; %Get just the predator and prey populations
y = zeros(Ny, T); %will hold observables

%Generate measurements
for t = 1:T
    y(:,t) = feval('MeasurementFcn', x(:,t), [], [], [], []); %Call measurement function
end


fprintf('Now estimating the model parameters...')
fprintf('\n')



x0_real = [y(1,1); y(2,1)];
params0 = [0.62526; 0.6607; 0.1896; 0.0468];
P0_x_real = [20 0;0 15]; %initial covariance for states
P0_params_real = [0.5 .001 .001 .001; .001 0.7 .001 .001; .001 .001 0.02 .001; .001 .001 .001 0.02]; %initial covariance for parameters


Q_x_real = []; %initialize structs for process noise
Q_param_real = [];
R_x_real = []; %Initialize structs for measurement noise
R_param_real = [];


%SET UP STATE DATA STRUCTURE AND INFO FOR STATES
InfDS_X_real = [];
InfDS_X_real.spkfParams = [alpha beta kappa];
InfDS_X_real.hfun = @MeasurementFcn;
InfDS_X_real.par = params0;
InfDS_X_real.obsdim = Ny;
InfDS_X_real.Nsubj = Nsubj;
InfDS_X_real.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS_X_real.NxNoPar = 2;	          %Number of state variables
InfDS_X_real.partQflag = 0; % A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS_X_real.Xdim = 2;
InfDS_X_real.ffun = @Lotka_Volterra_Model;
InfDS_X_real.Xbool = 1;

R_x_real.cov = diag([30 15]); %set measurement noise noise
Q_x_real.cov = 180 * [1.3444 0.1594;0.1594 0.4845]; %value adapted from output of Process_Noise_Estimation.m


%SET UP PARAMETER DATA STRUCTURE AND INFO FOR PARAMETERS
InfDS_Param_real = [];
InfDS_Param_real.spkfParams = [alpha beta kappa];
InfDS_Param_real.hfun = @MeasurementFcn;
InfDS_Param_real.par = params0;
InfDS_Param_real.obsdim = Ny;
InfDS_Param_real.Nsubj = Nsubj;
InfDS_Param_real.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS_Param_real.NxNoPar = 2;	          %Number of state variables
InfDS_Param_real.partQflag = 0; % A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS_Param_real.Xdim = 4;
InfDS_Param_real.ffun = @ParamTransitionFunction;
InfDS_Param_real.Xbool = 0;

R_param_real.cov = [30 0;0 15]; %set measurement noise (should be same as R_x_real.cov)
Q_param_real.cov = .01 * P0_params_real; %set using formula from ____

xhat_real = zeros(Nx, T); %will hold state estimates
paramshat_real = zeros(Npar, T); %will hold parameters estimates
xhat_real(:,1) = x0_real; %set first time entry
paramshat_real(:,1) = params0; %set first time entry

%LOOP TO PERFORM THE DUAL UKF
for t = 2:T
   [param_next, P_param_next] = UKF_Params(params0, P0_params_real, Q_param_real, R_param_real, y(:,t), [], [], InfDS_Param_real, x0_real, Nx); %estimate parameters
   %set new values based on output
   params0 = param_next;
   P0_params_real = P_param_next;
   InfDS_X_real.par = params0;
   InfDS_Param_real.par = params0;
   [x_next, Px_next] = UKF(x0_real,P0_x_real,Q_x_real, R_x_real, y(:,t),[],[],InfDS_X_real); %estimate states
   %set new values based on output
   x0_real = x_next;
   P0_x_real = Px_next;
   xhat_real(:, t) = x_next;
   paramshat_real(:,t) = param_next;
end

param_final = paramshat_real(:, end);
prey_error = xhat_real(1,:) - x(1,:);
prey_error = prey_error.^2;
prey_error = sum(prey_error(:, 45:91));
MSE_prey = prey_error / 47;

predator_error = xhat_real(2,:) - x(2,:);
predator_error = predator_error.^2;
predator_error = sum(predator_error(:,45:91));
MSE_predator = predator_error /47;

param_final
MSE_prey
MSE_predator


%CALL FUNCTION TO CREATE REAL DATA FIGURES
Dual_UKF_RealData_Figures(xhat_real,paramshat_real, T, x, truepar, y, param_final);





