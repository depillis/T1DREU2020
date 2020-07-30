%MAIN FILE FOR THE JOINT UKF ON THE PREDATOR-PREY SYSTEM
%THE UKF FUNCTIONS ARE ADAPTED FROM CHOW-FERRER AND THEN HAVE BEEN FULLY
%IMPLEMENTED BY DANIEL SHENKER AND RACHEL WANDER


% Portions of this script was adapted from previous scripts written by:
%       Nando de Freitas      (jfgf@cs.berkeley.edu)
%       Rudolph van der Merwe (rvdmerwe@ece.ogi.edu)
%       Kevin Murphy 
% So please acknowledge their previous work
% The appropriate references and descriptions of the algorithm used herein can be found in
% Chow, Ferrer & Nesselroade (2005). An unscented Kalman filter approach to the estimation of nonlinear dynamical systems models.
% DATE     : 18 July 2005

clear all;
clc;
echo off;

% INITIALISATION AND PARAMETERS:
% ==============================

no_of_runs = 1;           % number of simulation runs
T = 91;  %Number of data points
Ny = 2;  %Number of observables
Nx = 6;  %Number of values in x vector (# states + # parameters)
Nstates = 2; %number states
Npar = 4;   %number parameters
Nsubj = 1; %Number of subjects to run the UKF on, in our case unused
truepar = [0.62526 0.6607 0.1896 0.0468];
truepar_generating = [0.62526 0.6607 0.1896 0.0468 0.0001 0.0001 20 20];


Q = []; %initialize struct


param_final = zeros(Npar, no_of_runs); %matrice for the final parameter values
MSE_prey = zeros(1, no_of_runs); %matrice for the prey mean-squared error
MSE_predator = zeros(1, no_of_runs); %matrice for the predator mean-squared error



%SET UKF PARAMETERS

alpha = 10e-3; % UKF : point scaling parameter
beta = 2;  % UKF : scaling parameter for higher order terms of Taylor series expansion 
kappa = 0;    % UKF : sigma point selection scaling parameter (best to leave this = 0)

%**************************************************************************************

% MAIN LOOP

for j=1:no_of_runs,

  rand('state',sum(100*clock));   % Shuffle the pack!
  randn('state',sum(100*clock));   % Shuffle the pack!
  HLData = load('HaresLynxData.mat');  %Load dataset
  rawData = HLData.Lotka_Volterra_Data;
  
  x(1:2,:) = rawData(:, 2:3)'; %Get just the predator and prey populations
  y = zeros(Ny, T); %set up empty matrix for measurements
  
  %Generate measurements
  for t = 1:T
      y(:,t) = feval('MeasurementFcn', x(:,t), [], [], [], []); %Call measurement function
  end


fprintf('Now estimating the model parameters...')
fprintf('\n')

x0 = [y(1,1); y(2,1); 0.62526; 0.6607; 0.1896; 0.0468];
P0 = [20 .01 .01 .01 .01 .01;
      .01 15 .01 .01 .01 .01;
      .01 .01 0.5 .01 .01 .01;
      .01 .01 .01 0.3 .01 .01;
      .01 .01 .01 .01 0.02 .01;
      .01 .01 .01 .01 .01 0.02];

par0 = truepar; %Set parameters to actual values
%noise to help speed convergence   
Q.cov = diag(diag([repmat(10e-5,Nx)])); % Add a little bit of process
R = []; %Initialize covariance of measurement noise as struct

%Set up struct of information needed to run the UKF
InfDS.spkfParams = [alpha beta kappa]; %UKf scaling parameters
InfDS.hfun = @MeasurementFcn; %The measurement function H
InfDS.par = par0; %True values for model parameters
InfDS.obsdim = Ny; %Number of observables
InfDS.Nsubj = Nsubj;
InfDS.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS.NxNoPar = 2;	          %Number of state variables
InfDS.partQflag = 0;	% A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS.Xdim = 6; %x dimension 3 for 3 states

R.cov = diag([30, 15]); %Set R covariance matrix using known parameter values for noise (determined ad hoc)
%covariance of measurement noise
Q.cov=diag(diag([repmat(10e-5,Nx)])); %Initial covariance of process noise
starttime = cputime; %Initialize starttime (to calculate how long UKF takes to run)

[xhat, Px, xh_, Px_, yh_, inov, Py, KG, Pxall, Pyall, Px_all]= UKF(x0,P0,Q, R, y,[],[],InfDS); %Use UKF to do state estimates
time_ukf = cputime-starttime %Time that UKF took
xhat %print out all state estimations

%FILL IN RUN INFORMATION
param_final(:, j) = xhat(3:6, end); %set the final parameter guesses to be the final parameter estimation

%Calsulate mean squared error for prey
prey_error = xhat(1,:) - x(1,:);
prey_error = prey_error.^2;
prey_error = sum(prey_error(:,45:91));
MSE_prey(:,j) = prey_error / 47;

%Calculate mean squared error for predator
predator_error = xhat(2,:) - x(2,:);
predator_error = predator_error.^2;
predator_error = sum(predator_error(:, 45:91));
MSE_predator(:,j) = predator_error / 47;
end

param_final %Print out final parameters
MSE_prey %Print out prey mean-squared error
MSE_predator %Print out predator mean-squared error


%Create plots for real data
Joint_UKF_Real_Data_Figures(xhat, x, T, truepar, y, param_final)


