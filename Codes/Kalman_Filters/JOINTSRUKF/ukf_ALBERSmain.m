%%%%%% LdeP Comments: %%%%%%%
% Square Root Unscented Kalman Filter - Main calling code.
%	ukf_ALBERSmain.m (this file)
% Other files needed to run this script: 
%	AlbersMeasFcn.m
%	AlbersODEwParams.m
%	AlbersParamVals.m
%	AlbersStateFcn.m
%	PlotMyUKFResults.m
%	srukf.m
%
% LdeP April 2020 - Updated to run the Albers T2D system. One main change to the call of srukf.m is that we need to pass through a time-step argument. This allows us to use built-in MATLAB ODE solvers to find a solution x[k+1] based on x[k].
%
% The Square-Root Unscented Kalman Filter for State and Parameter-Estimation 
% History: 
%	First ukf.m implments UKF written 2008 (YC)
%	Updated to ukf.m to use square-root UKF written 2017 (ZH)
%	Updated to squkf.m to perform joint state-parameter estimation, and to work with ODE system 2020 (LdeP)
%
%%%%%% LdeP Comments: %%%%%%%
%
% First written by Yi Cao, Cranfield University, 2008. Used Julier & Uhlman 2004 paper
% Later updated to square root version By Zhe Hu at City University of Hong Kong, 05/01/2017
% Used Reference: R. van der Merwe and E. Wan. 
% THE SQUARE-ROOT UNSCENTED KALMAN FILTER FOR STATE AND PARAMETER-ESTIMATION, Conference Paper in Acoustics, Speech, and Signal Processing, 1988. ICASSP-88., 1988 International Conference on · February 2001
%
% LdeP April 7, 2020 - Albers implemented with SKUKF working;
%			- Eliminated the use of "svec". It seemed to be introducing errors.
% LdeP April 7, 2020 - Continuing on Albers implemented with SKUKF.
% LdeP April 6, 2020 - Using VdP code to try to build Albers code.
% LdeP April 1, 2020 - Update to try VanDerPol equation as a test application.
%
%%%%%%
%
% Updated by Lisette de Pillis to perform joint estimation (states and parameters)
% Using information from the 2001 Wan and van der Merwe paper, explanation of the problem follows;
% Notation may be changed a bit.
%
% Helpful MATLAB examples of problems using Kalman filters can be found https://www.mathworks.com/help/control/examples.html?category=state-estimation and for an unscented Kalman filter, look here https://www.mathworks.com/help/control/ug/nonlinear-state-estimation-using-unscented-kalman-filter.html
%
%  State of a discrete-time nonlinear dynamic system:
%	x[k+1] = F(x[k], u[k]) + q[k]		% State with process noise q
%	y[k]	= H(x[k]) + r[k]		% Measure with measurement noise r
%	
%	where
%	x[k] represents the UNOBSERVED (TRUE) state of the system at step k
%	u[k] represents exogenous input (usually some kind of nonlinear forcing term). Input u[k] is optional and represents additional inputs to F, for instance system inputs or parameters. Input u[k] can be specified as zero or just represent more function arguments. 
%	q[k] is the process noise (assumed additive in this case).
%	y[k] is the observed measurement signal.
%	r[k] is the measurement noise. 
%	F(): State transition function F returns x[k+1] given x[k], etc. This can be the solution to an ODE system. The user will have to input the time step represented by going from step k to step k+1.
%	H(): Measurement function H returns what part of the state variable vector we can actually measure. In the case of the Albers ODE, we probably want H to return actual values of glucose and the three parameters we'll track. We can add some measurement noise r[k] and put that into our measure vector y[k].
%
% 	We often assume that process and measurement noises are standard Gaussian noises, and their covariance matrices are all known.
%
%  Parameter estimation (AKA "System Identification") assumes that our system is described by
%	y[k] = G(x[k], w)
%  where x[k] is the input, y[k] is the output, and G() is a nonlinear map parameterized by parameter vector w. This usually requires that we provide some "training pairs" of known inputs with desired outputs, {x[k], d[k]}.  The "error of the machine" is then 
%	e[k] = d[k] - G(x[k],w). This can be rewritten to look a lot like our state-space problem, and becomes
%	w[k+1] 	= w[k] + qp[k]
%	d[k]	= G(x[k], w[k]) + e[k]
%	
%	Since w[] is the vector of parameters, qp[k] (q[k] above) is the process noise at step k. 
%	The choice of d[k] corresponds to a nonlinear observation on w[k]. 
%%%
%	The SR-UKF process itself has these elements:
%	P: State covariance matrix.
%	S: Square root of the covariance matrix -> P = S S'. Can be found with a Cholesky factorization.
%	In UKF, the full state-covariance matrix P is propagated.
%	In SR-UKF, the square root matrix S is propagated.
%
%%%%%%%%%%%%%%%%%%%%%%%

n_numStates	=9;      %dimension of state, 6 states, 3 parameters
m_numMeas	=4;      %dimension of measurement - LdeP We'll start with measuring only glucose state G, and 3 parameters E, Vi, Ti.


% LdeP STATE FUNCTION
% Input current state x. Iterate forward from time step k to time step k+1. Return updated x.

% LdeP MEASUREMENT FUNCTION
% Input current state xk, return yk=xk(1) for Glucose only, or yk = [xk(1) xk(7) xk(8) xk(9)] for Glucose and three parameters.

% LdeP Initial state guess - Update for Albers
% Initial conditions taken from Albers code
% I_p = 200
% I_i = 200
% G = 12000 <-- LdeP This is total glucose (mg)
TotalGlucose = 12000; % LdeP Introducing this so we can modify this and GlucosePerDL will follow
% h_1 = 0.1
% h_2 = 0.2
% h_3 = 0.1
% LdeP Append paramter ICs
E = 0.2;  % (= 0.2) l/min : insulin exchange rate between remote and plasma compartments
Vi = 11;  % (= 11) l : interstitial volume
ti = 100; %  (= 100) min : time constant for remote insulin degradation
initialStateGuessODE = [200; 200; TotalGlucose; 0.1; 0.2; 0.1; E; Vi; ti]; %LdeP This is the same set of initial values used in the Albers code. The parameter value are set to their actual value used for simulations.

% System has 10L of glucose = 100dL
% We compute total mg of glucose in system
% but plot mg dL^(-1)
% So we have to divide glucose by 100 before plotting
GlucosePerDL = TotalGlucose/100; %LdeP 

%NumberOfStates = numel(initialStateGuessODE); %Ldep keep track of the number of states you have (9 states means we are also tracking parameters)

% SPECIFY INITIAL and FINAL TIME IN MINUTES
timeInitial = 0; %LdeP ODE is time invariant, so we can start at 0
% timeFinal = 12960; % t = 9 days worth of minutes = 12960 minutes
% timeFinal=500; %LdeP Try a different length of time. Albers went to over 1200.
timeFinal=1440; %LdeP 24 hours in minutes

% SOLVE FOR TRUE SOLUTIONS before setting up the UKF (we can test UKF against these "true" solutions)
% LdeP Get an actual good true solution xTrue through ODE45 function call
[timeSteps,xTrue]=ode45(@AlbersODEwParams,[timeInitial:1:timeFinal],initialStateGuessODE); %LdeP Force timestep 1
xTrueGlucosePerDL=xTrue(:,3)/100; % LdeP Convert Glucose element to be mg per dL

% SET UP UKF INITIAL VALUES and NOISE STANDARD DEVIATION VALUES
initialStateGuessUKF = [200; 200; GlucosePerDL; 0.1; 0.2; 0.1; E*0.25; Vi*1.5; ti*0.5]; %LdeP State vector with glucose scaled to mg/dL; Parameter value initial state guesses are deliberately scaled away from their "actual" values, to simulate a reasonable guess and see how well we do with getting parameter values to converge to their true value.
%%
% LdeP Update to look more like Albers example: 
% Recall: sqrt(Variance) = StDev
% R = 0.2;  % Variance of the measurement noise v[k],  original
R = 2.25;  % Variance of the measurement noise 
sqrtR = sqrt(R); % Standard deviation of measurement noise

% LdeP We are currently using only one common process noise value.
% LdeP Once you get this working, see if you can update this to match the
% collection of process noise values in the most recent AlbersUKF_ParState.m code, 
% which is: diag([0.02 0.1 0.04 0.2 0.5 0.01 0.0001 0.0001 .0001]); %LdeP These are some random small process noise values
q = 0.1;     % Standard deviation of process noise (using baseline code value)
% q = 0; %LdeP seeing if we can get a more accurate "true" measure
Rs=sqrtR*eye(m_numMeas);	% stdev matrix of measurement noise
Qs=q*eye(n_numStates);		% stdev matrix of process noise

% LdeP Original from baseline code
% r=0.1;    %std of measurement
% q=0.1;    %std of process 
% Rs=r*eye(m);  % stdev of measurement  
% Qs=q*eye(n); 	% stdev matrix of process
%%

%%%%% Recall our system setup:
%  State of a discrete-time nonlinear dynamic system:
%       x[k+1] = F(x[k], u[k]) + q[k]           % State with process noise q
%       y[k]    = H(x[k]) + r[k]                % Measure with measurement noise r
%
% LdeP x becomes our state with some process noise (standard deviation q)
rng('default'); % LdeP Matlab said use this to fix the random number generator for reproducible results
x = initialStateGuessUKF + q*randn(n_numStates,1); % Initialize state vector to just be initial state with process noise; We represent x=F(x[k],u[k])+q[k]. It's xhat[k|k-1], where k=1;

S = eye(n_numStates);	% initial square root of state covariance

%%%%% Time stepping
% Recall that above we already specified timeInitial and timeFinal (in minutes)
% LdeP April 2020 - note that if dt is our chosen sample time-step size, and number of computational steps Nsteps=20, then timeFinal = 20*dt is where our simulation will end.
%
% LdeP Time step dt - a future iteration may adapt to have irregular sample time steps
%
dt = 1;   % LdeP Pass this through to AlbersStateFcn.m
Nsteps  = timeFinal/dt; % Total number of computational steps we'll take.


%% Allocation memory by creating zero vectors
xV = zeros(n_numStates,Nsteps); % UKF estimate of state vector, LHS of our discrete time equation
zV = zeros(m_numMeas,Nsteps);	% Measured values of system  (LdeP like yMeas in AlbersUKF code)

%LdeP yTrue taken from previous AlbersUKF code 
yTrue = [xTrueGlucosePerDL'; xTrue(:,7)'; xTrue(:,8)'; xTrue(:,9)']; % Glucose G is the only state we currently measure, and we also assume we measure E, Vi, ti 
zV = yTrue(:,1:Nsteps) + (sqrtR*randn(m_numMeas,Nsteps));

for k=1:Nsteps
    % z = h(s) + r*randn(m_numMeas,1);                     % measurments
    % Moving the computation of z to outside the loop - it should be a
    % close approximation to true values.
  	e(:,k) = zV(:,k)-AlbersMeasFcn(x,dt); % Calculate error at step k between measure z and measure based on last computed state
  	% LdeP putting the computation of zV outside of this loop
    	% zV(:,k)  = z;	% save measurment; LdeP: z should look like the true solution plus a little noise with stdev r
  	% [x, S] = ukf(f,x,S,h,z,Qs,Rs);           % ukf 
  	% LdeP: AlbersStateFcn returns state x advanced by time step dt
  	% LdeP: input x is our UKF computed state at current time, step k;
  	% LdeP: output x is our UKF computed state at next step k+1, time advanced by time step dt
  	% LdeP: input S is the "a priori" estimated the square root of state covariance
  	% LdeP: output S is the "a posteriori" square root of the state covariance
  	% LdeP: input z is our measurement (current) at step k
  	% LdeP: input Qs is process noise standard deviation
  	% LdeP: input Rs is measurement noise standard deviation
    	% LdeP: input dt - we will use ODE45, so a time argument must be sent in.
  	[x, S] = srukf(@AlbersStateFcn,x,S,@AlbersMeasFcn,zV(:,k),Qs,Rs,dt);  % LdeP: ukf call to Albers measurement function, returns elements 3,7,8,9 of state variable x (glucose and 3 parameters)

  	xV(:,k) = x;% save estimate at iteration step k; x is the UKF computed
		    % guess of our state function given the measured data, 
		    % and the noise matrices
    k        
    
end


%%%%%%%%
% Plots - the script below assumes the parameters have the same names as above. In future, we can change the plot script to a function and pass values.
%%%%%%%%
PlotMySRUKFResults;

