%DOES PARAMETER ESTIMATION FOR 11 MICE FROM LI ET AL DATASET. WE SIMULATE
%USING THE GIVEN ODE'S UP UNTIL THE FIRST AVAILABLE TIME POINT. THEN THE
%DUAL UKF IS RUN ON EACH MOUSE TO GET PARAMETER ESTIMATES. THE DUAL UKF
%ALGORITHM IS RUN 5 TIMES ON EACH MOUSE AND THEN THE BEST OPTION, BASED ON
%LEAST SQUARED ERROR, IS CHOSEN AS THE FINAL SET OF PARAMETERS FOR THAT
%MOUSE

%IN THIS FILE, WE KEEP MU_R AND MU_E AS THE SAME VALUES


clear all;
clc;
echo off;

%Set general parameters

alpha = 10e-4; % UKF : point scaling parameter
beta = 2;
kappa = 0;

Ny = 1;  %Number of observables
Nx = 12;  %Number of states
Nstates=12; %number states
Npar = 41;
figure_params = zeros(Npar + 1, 11); %this will hold the best fit parameter values for each mouse
fig_2_3_4_5_Parameters; %set up default parameters
wave = wave_basal; %set wave value

Nruns = 5; %number of times to run the filter

%LOOPING OVER ALL 11 MICE DATASETS
for w = 1:11
file_num = w; %Mouse you want to use


%---LOAD THE DATASET----%
num_str = int2str(file_num);
file = strcat('dat',num_str,'.csv'); %create file name
X = readtable(file);
X = X{:,:};
X = flip(X);
X = X'; %set up table as wide
numReadings = size(X,2);
%---------------------%


%Set initial parameters
InitialParamGuess = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8);
                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21; 0.018];
P0_param_initial = diag([200 0.05 0.005 0.05 10e-10 10e-10 10e-12 10e-12 10e-4 10e-3 10 2 2 .02 5 5 10 175 10^-3 ...
                      10^-3 10^-3 10^-3 10^-3 175 10e-11 10e-11 10e-4 1 4 10e-3 10^-5 10^-5 10^-5 10^-3 ...
                      200 1 1 1 5 15 .0001]);
    
            
all_param_estimates = zeros(Npar + 1, Nruns); %this will hold all the final parameter estimates for each run
params_alltime = zeros(Npar, Nruns * numReadings); %this will hold the parameter estimates after each reading
%--------RUN LOOP FOR HOW MANY RUNS YOU WANT ON EACH MOUSE ---------------------------------------%
for i = 1:Nruns
    %-----LOAD DATASET------%
    num_str = int2str(file_num);
    file = strcat('dat',num_str,'.csv'); %create file name
    X = readtable(file);
    X = X{:,:};
    X = flip(X);
    X = X'; %set up table as wide
    %-------------------------%
    
    %get first and last times available
    first_time = X(1,1);
    last_time = X(1,end);
    
    InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0]; %Set initial state values
    
    truepar = InitialParamGuess;
    T = first_time;
    tspan = [0 T]; %span for simulating up to first reading available
    options=odeset('RelTol',1e-12,'AbsTol',1e-12); %set tolerance levels (helps ODE solver)
    sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, truepar), tspan, InitialState,options); %Simulate ODE's up until first time point
    initial_state = deval(sol,first_time); %update initial states after simulation
    
   
    
    numReadings = size(X,2); %number of available readings
    X_final = X; %set name for notation
   
    
    %-------DO PARAMETER ESTIMATES------------%
  
    
    x0 = initial_state; %initial state estimate
    T = numReadings;  %Number of data points
    params0 = InitialParamGuess; %initial parameter estimate
   
    %set up vector of observables
    y = zeros(Ny, T);
    for l = 1:T
        y(:,l) = feval('MeasurementFcn2', X_final(2,l), [], [], [], []); %Call measurement function
    end
    

    
    %----ESTIMATE MEASUREMENT NOISE-----%
    tspan_noise = [first_time * 7 last_time * 7]; %solve over span where data is available
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, truepar), tspan_noise, x0,options); %Use ODE solver
    x_predict_noise = zeros(Nstates, T);
    for t=1:numReadings %Generate data for each time point
        x_predict_noise(:,t) = deval(sol, 7*X(1,t)); %Get ODE solution for time 7*t
    end
    meas_noise = y - x_predict_noise(6,:);
    noise_var = var(meas_noise); %variance of the measurement noise
    
    %------------------------------------%
    
    %Set initial covariances
    P0_x = diag([repmat(1,Nstates,1)]); %state covariance
    P0_params = P0_param_initial;
    
                  
                  
    %--------SET UP STRUCTS WITH NOISE AND GENERAL INFO BOTH FOR PARAMETER AND STATE UKF FUNCTIONS--------%             
    Q_x = [];
    Q_param = [];
    R_x = []; %Initialize structs for R 
    R_param = [];
    %SET UP STATE DATA STRUCTURE AND INFO FOR STATES
    InfDS_X_real = [];
    InfDS_X_real.spkfParams = [alpha beta kappa];
    InfDS_X_real.hfun = @MeasurementFcn;
    InfDS_X_real.par = params0;
    InfDS_X_real.obsdim = Ny;
    
    InfDS_X_real.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
    InfDS_X_real.NxNoPar = 2;	          %Number of state variables
    InfDS_X_real.partQflag = 0; % A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
    InfDS_X_real.Xdim = 12;
    InfDS_X_real.Xbool = 1;
    R_x.cov = diag(noise_var); %set noise
    Q_x.cov = diag([10^3 10^3 10^3 10^3-300 100 50 2 10^2 10^3 10^4 100 200]);
    
    InfDS_Param = [];
    InfDS_Param.spkfParams = [alpha beta kappa];
    InfDS_Param.hfun = @MeasurementFcn;
    InfDS_Param.par = params0;
    InfDS_Param.obsdim = Ny;
    
    InfDS_Param.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
    InfDS_Param.NxNoPar = 2;	          %Number of state variables
    InfDS_Param.partQflag = 0; % A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
    InfDS_Param.Xdim = 41;
    InfDS_Param.Xbool = 0;
    
    
    
    %parameter noise covariances
    R_param.cov = diag(noise_var);
    Q_param.cov = .01 * P0_params;
    %-------------------------------------------------------------------%
    
    
    
    xhat = zeros(Nx, T); %to hold state estimates
    paramshat = zeros(Npar, T); %to hold parameter estimates
    %initialize with first guesses
    xhat(:,1) = x0;
    paramshat(:,1) = params0;
    for t = 1:T
        if (t == 1)
            tspan = 0; %if first reading then no time span
        else
            tspan = X_final(1,t) - X_final(1,t-1); %calculate time between readings
        end
        [param_next, P_param_next] = UKF_Params(params0, P0_params, Q_param, R_param, y(:,t), [], [], InfDS_Param, x0, Nx, tspan); %run parameter UKF
        %update values with output from parameter UKF
        params0 = param_next;
        P0_params = P_param_next;
        InfDS_X_real.par = params0;
        InfDS_Param.par = params0;
        [x_next, Px_next] = UKF(x0,P0_x,Q_x, R_x, y(:,t),[],[],InfDS_X_real, tspan); %run state UKF
        %update values with output from state UKF
        xhat(:, t) = x_next;
        paramshat(:,t) = param_next;
        x0 = x_next;
        P0_x = Px_next;
    end
    
    P0_param_initial = P0_params; %set for the next iteration
    all_param_estimates(2:42,i) = paramshat(:,end); %put parameters in rows 2 to 42
    final_params = paramshat(:, end); %final parameters for that run
    all_param_estimates(1,i) = determineError(X,final_params, wave); %get LS error for this run and put in first row
    InitialParamGuess = paramshat(:,end); %Update parameter guess for the next run
    params_alltime(:,numReadings * (i-1) + 1:numReadings * i) = paramshat; %add entire set of parameter estimates to set of all

 
end
errors = all_param_estimates(1,:); %get the 5 errors for this mouse
[M, I] = min(errors); %get index of min error
column_min = I;
figure_params(:, file_num) = all_param_estimates(:, column_min); %set final parameters for that mouse as best fit
end
figure_params %print out the final parameters


