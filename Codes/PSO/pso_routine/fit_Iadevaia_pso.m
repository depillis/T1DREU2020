%Author: B. Shtylla
%Last Updated: 04/24/2020

%This code implements the matlab built-in particle swarm optimization
%routine in order to fit an ODE model. 
%To use this code you need access to the matlab global optimization
%package found here: https://www.mathworks.com/products/global-optimization.html


%The structure of the code is as follows: 
%-Input the parameter vector, P as a global parameter that is read in the other
%helper functions. 
%-Read the data points from csv files (requires both data values and time
%points)
%-Declare an objective function, using the squared_error function
%   -Error function requires parameter vector: P_test; time: data time points;
% points: data values. All of these are used to compute a squared
% difference.
%   -Objective function minimizes the squared error_function by varying the
% parameter vector. The objective came from this Github repo:
% https://github.com/ExpectozJJ/Particle-Swarm-Optimization-with-Passive-Congregation-


%This code parametrizes the GLP-1 QSP (reduced) model published by Iadevaia et al
%(2010), paper is found in the references folder. The code of the methods is 
%found in the supplementary files. For ease we use their fitted sets in
%sfile_5.xls as initial guess for our fitting. We fix the initial
%conditions- meaning the last few entries in P are not changed. 



clear all;
close all;

%declare parameter vector as a global quantity to be updated in the 
%appropriate location
global P;

%load the data points
data_points = csvread('data/synthetic_data_points_Nsample15.csv');
%load the data times
data_time = csvread('data/synthetic_data_time_Nsample15.csv');


%Read the inital guess for the parameter vector
param=xlsread('data/sfile_5.xls');
P_test=param(1:97,1);%Pick the first set of parameters



%Set the number of variables to be equal to the unknowns-not sure about
%this (seems to need to be a multiple of 10)
nvars = length(P_test);

%Set the lower and upper bound for the parameter intervals by expanding
%around the initial values
lb = P_test-0.3*P_test(P_test>0); %Only subtract if initial values are greater than zero
ub = P_test+0.5*P_test;

%Set up the solver options; swarm size has the number of particles
%Display progress and parallelize the search
options = optimoptions('particleswarm','SwarmSize',10,'MaxIterations',400,...
    'UseParallel',true,'Display','iter','DisplayInterval',20);



%Declare the objective function
objfun = @(k) squared_error(k,data_time,data_points); 

%apply the optimization routine output solutions in pso_param
[pso_param,fval,exitflag,pso_output] = particleswarm(objfun,nvars,lb,ub,options);

%save the output
csvwrite('pso_param_set2_Nsample15.csv',pso_param);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%helper functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function error = squared_error(P_test,time, points)
%   Calculates the squared error between a set of data points and a
% series of points representing a function output from ODE solver.
%
% Inputs:
%   time   - A 1D array containing the time points of the data
%   points - A 1D array of the data values
%   P      - Parameter and (possibly) initial conditions being varied by the pso
%   algorithm
% Output:
%   error = the sum of the squared residuals of the data vs. the function

      %Assign the parameter values using P
      global P;
      P=P_test;
      pso_parameters;
      
      %Declare initial conditions for ODE solver
      pso_init_vec;

      % Solve the ODE
      [fn_time, fn_points]=ode15s(@IGFRODE,time,init_vec);
    
    
    
      resid = (fn_points-points).*(fn_points-points);
 
    % Square the error values and sum them up
      error = sum(sum(resid));
    
end