% PURPOSE : Fitting the Lorenz system model to noisy time series data using
% particle filter with UKF proposal

% AUTHOR: Sy-Miin Chow (schow@nd.edu)
% Portions of this script was adapted from previous scripts written by:
%       Nando de Freitas      (jfgf@cs.berkeley.edu)
%       Rudolph van der Merwe (rvdmerwe@ece.ogi.edu)
%       Kevin Murphy 
% So please acknowledge their previous work
% The appropriate references and descriptions of the algorithm used herein can be found in
% Chow, Ferrer & Nesselroade (2005). An unscented Kalman filter approach to the estimation of nonlinear dynamical systems models.
% DATE     : 18 July 2005



%DSS - DC's in the pancreas
%Mues_r - effector and regulatory cell interaction
%Mues_e - effector and regulatory cell interaction
%Q_panc - pancreas volume
%SI - rate glucose taken up in proportion to insulin
%GI - insulin saturation point
%alpha_B - beta cell production rate


clear all;
clc;
echo off;

% INITIALISATION AND PARAMETERS:
% ==============================

no_of_runs = 1;           % number of simulation runs
T = 40;  %Number of data points
Ny = 1;  %Number of observables
Nx = 53;  %Number of values in x vector (# states + # parameters)
Nstates = 12; %number states
Npar = 41;   %number parameters
Nsubj = 1; %number of subjects you have data for

%params for ODE
T1D_ODE_Parameters;

Nruns = 1; %number of times to run the filter - just to make sure it runs
file_num = 6; %Mouse you want to use


num_str = int2str(file_num);
file = strcat('dat',num_str,'.csv'); %create file name
X = readtable(file);
X = X{:,:};
X = flip(X);
X = X'; %set up table as wide
numReadings = size(X,2)-1;



truepar = [10^5; 0.0334; sqrt(20000); .72; .194; aEmday/5000; aEmday/5000]; %True parameter values - taken form original parameters file
InitialStateGuess = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
InitialParamGuess = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8); 
                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21];




%turn wave on or off
wave = wave_basal;

Q = []; %initialize struct

%SET UKF PARAMETERS

alpha = 10e-4; % UKF : point scaling parameter
beta = 2;  % UKF : scaling parameter for higher order terms of Taylor series expansion 
kappa = 0;    % UKF : sigma point selection scaling parameter (best to leave this = 0)

%**************************************************************************************

% MAIN LOOP

 


fprintf('Now estimating the model parameters...')
fprintf('\n')



P0_x = diag([repmat(1,Nstates,1)]);
    
P1 = zeros(41,12);
P1 = P1 + .0001;
P2 = zeros(12,41);
P2 = P2 +.0001;
P0_params = diag([10^1 0.0001 0.00001 0.0001 10e-11 10e-10 10e-12 10e-12 ...
                   10e-9 10e-6 1 .1 .1 .0001 .1 .1 ...
                   .2 50 10^-6 10^-6 10^-6 10^-6 10^-6 100 10e-11 ...
                   10e-11 10e-6 .1 .01 10e-6 10^-8 10^-8 10^-9 ...
                   10^-5 10^-5 50 0.0001 0.0001 .05 10 8]);
P0 = [P0_x P2;  P1 P0_params];
par0 = truepar;  
R = []; %Initialize struct


%SET UP DATA STRUCTURE OF MODEL INFO
InfDS.spkfParams = [alpha beta kappa];
InfDS.hfun = @MeasurementFcn;
InfDS.par = InitialParamGuess;
InfDS.obsdim = Ny;
InfDS.Nsubj = Nsubj;
InfDS.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS.NxNoPar = 12;	          %Number of state variables
InfDS.partQflag = 0; % A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS.Xdim = 53;


Q_x = diag([10^3 10^3 10^3 10^3-300 100 50 2 10^2 10^3 10^4 100 200]);
Q_param = .01 * P0_params;
Q.cov=[Q_x zeros(12,41); zeros(41,12) Q_param];
all_param_estimates = zeros(Npar,Nruns);
params_alltime = zeros(Npar, Nruns * (numReadings));



file = strcat('dat',num_str,'.csv'); %create file name
X = readtable(file);
X = X{:,:};
X = flip(X);
X = X'; %set up table as wide
T = X(1,1);
InitialStateGuess = Get_Initial_Conditions(T);


for i = 1:Nruns
    x0 = InitialStateGuess;
    x0(13:53) = InitialParamGuess;
    [g Time] = size(X);
    xhat_fake = zeros(Nx,Time); %allocate memory for output
    xhat_fake(:,1) = x0;
    y_fake = X(2,:);
    y_fake = y_fake(2:Time);
    y_fake = y_fake';
    tspans = zeros(Time-1,1);
    
    if i == 1
        %----ESTIMATE MEASUREMENT NOISE-----%
        tspan_noise = [X(1,1) * 7 X(1,end) * 7]; %solve over span where data is available
        options=odeset('RelTol',1e-12,'AbsTol',1e-12);
        sol = ode15s(@(t, y) ODE_allparams(t, y, f1n, f2n, wave, InitialParamGuess), tspan_noise, x0,options); %Use ODE solver
        x_predict_noise = zeros(Nstates, Time);
        for t=1:numReadings %Generate data for each time point
            x_predict = deval(sol, 7*X(1,t));
            x_predict_noise(:,t) = x_predict(1:12); %Get ODE solution for time 7*t
        end
        x_predict_noise(6,:);
        meas_noise = y_fake - x_predict_noise(6,:);
        noise_var = var(meas_noise);
        R.cov = [noise_var(1)];
    end
    
    for t = 1:numReadings
        tspans(t) = X(1,t+1)-X(1,t);
    end
    tspans= 7.*tspans;
    
    starttime = cputime;
    [xhat_fake, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall, Pyall, Px_all]= UKF_allparams(x0,P0,Q, R, y_fake,[],[],InfDS, tspans); %do the estimation
    time_ukf = cputime-starttime
    all_param_estimates(:,i) = xhat_fake(13:53, end);
    params_alltime(:,numReadings * (i-1) + 1:numReadings * i) = xhat_fake(13:53,:);
    InitialParamGuess = xhat_fake(13:53, end);
   
end
