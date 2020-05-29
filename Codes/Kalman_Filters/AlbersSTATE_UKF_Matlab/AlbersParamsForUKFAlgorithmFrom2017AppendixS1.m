% Author: LdeP
% Date: March 30, 2020
% Script to enter the Albers 2017 Dual UKF parameters (for the algorithm)
% See Albers2017S1Appendix.pdf

nx = 6; % Number of model states for the ultradian model that we are using
nw = 3; % Number of estimated model parameters (they fit E, V_i, t_i)
alphax = 0.4; % State sigma point constant
alphaw = 0.3; % Parameter sigma point constant
betax = 2; % State sigma point constant (beta = 2 if often used)
betaw = 2; % Parameter sigma point constant
kx = 3-nx; % State sigma point constant (problem: negative kx can lead to numerical issues, maybe don't use this)
kw = 0; % Parameter sigma point constant
lambdatilde = 0.9975; % Forget factor for parameter covariance
qv = 0.4; % State sigma point constant for process noise
qn = .01; % State sigma point constant for measurement noise

% Below, the vall and vobs values are those reported in S1
vall = [82.24 ;191.57 ;11235.00; 82.02; 82.27; 81.86]; % Ultradian model average state values given in Albers2017S1Appendix.pdf
vobs = vall(3);% Ultradian average observed state value (glucose total)

% Note: If we run the sample code for the glucose ultradian model (no
% filtering), and initial conditions given by
% vall_init = [200 200 12000 .1 .2 .1], then the average values of the
% solutions in X, given by taking mean(X(2:end,:)) give
% vall_average = [86.84   205.51   12247.88   86.78    86.72   86.66]
% so it's not clear where the vall in Albers Appendix S1 came from.

%LdeP Since the UKF algorithm is fed glucose data that has been scaled to
%per dL, we must divide GLUCOSE by 100:
vallscaled=vall; vallscaled(3)=vall(3)/100;
vobsscaled = vobs/100;
Rv = diag((qv*vallscaled).^2); %Assumed process noise covariance
Rn = diag((qn*vobsscaled).^2); %Assumed measurement noise covariance
