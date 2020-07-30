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

%MODIFIED BY DANIEL SHENKER AND RACHEL WANDER MAY 2020 TO USE FOR SOLELY
%STATE ESTIMATION

clear all;
clc;
echo off;

% INITIALISATION AND PARAMETERS:
% ==============================

no_of_runs = 1;           % number of simulation runs
T = 2000; %end time value
Ny = 3; %num observables
Nx = 3; %one x variable for each state
Nstates=3; %same as Nx 
Npar=6;  %number of parameters being estimated (sigma, rho, beta)
Nsubj=1; %number of subjects (i.e. number of sets of states)
truepar = [10 28 8/3 26 34 32]; %true values of parameters which lines up with table 2 in Chow paper [sigma rho beta Theta]
Q = []; %initializing empty struct
alpha = 10e-4; % UKF : point scaling parameter
beta = 2;  % UKF : scaling parameter for higher order terms of Taylor series expansion 
kappa = 0;    % UKF : sigma point selection scaling parameter (best to leave this = 0)

%**************************************************************************************

% MAIN LOOP

for j=1:no_of_runs,

  rand('state',sum(100*clock));   % Randomizing for use in the noise
  randn('state',sum(100*clock));   % Randomizing for use in the noise
  

% GENERATE THE DATA:
% ==================

x = zeros(Nx,T); %One row for each state, one column for each time
y = zeros(Ny,T); %One row for each observable, one column for each time
measureNoise = zeros(Ny,T); %Matrix for measurement noise
measureNoise(:,1) = [sqrt(truepar(4))*randn(1,1); sqrt(truepar(5))*randn(1,1); sqrt(truepar(6))*randn(1,1)]; %Initial measurement noise  
x(:,1) = [-5; 20; -5];  % Initial state guesses (3 of them)
y(:,1) = feval('LormeasO',x(:,1),[],[],[],[]) + measureNoise(:,1); %Initial observable values

for t=2:T, %Generate data for each time point
  measureNoise(:,t) = [sqrt(truepar(4))*randn(1,1); sqrt(truepar(5))*randn(1,1); sqrt(truepar(6))*randn(1,1)];  %Generate measurement noise using random numbers  
  x(:,t) = feval('LordynO',x(:,t-1),[],[],t,truepar(1:3)'); %Solve the system using Runge-Kutta
  y(:,t) = feval('LormeasO',x(:,t),[],[],[],[]) + measureNoise(:,t); %Solve for observables  
end;  

% PLOT THE GENERATED DATA:
% ========================
figure(1) %Plots states over time
p1=plot(1:T,x(1,:),'k','Linewidth',2);hold on;
p2=plot(1:T,x(2,:),'r-.','Linewidth',2);
p3=plot(1:T,x(3,:),'b+');hold off;
ylabel('Data','fontsize',15);
xlabel('Time','fontsize',15);
legend('x_{1t}','x_{2t}','x_{3t}');

figure(2) %Plots relationship between x1 and x2
subplot(121)
plot(x(1,:),x(2,:),'k','Linewidth',2)
ylabel('x_{2t}','fontsize',24);
xlabel('x_{1t}','fontsize',24);

figure(2) %Plots relationship between x1 and x3
subplot(122)
plot(x(1,:),x(3,:),'k','Linewidth',2)
ylabel('x_{3t}','fontsize',24);
xlabel('x_{1t}','fontsize',24);

figure(3) %Plots relationship between y1 and y2
subplot(121)
plot(y(1,:),y(2,:),'k','Linewidth',2)
ylabel('y_{2t}','fontsize',24);
xlabel('y_{1t}','fontsize',24);

figure(3) %Plots relationship between y1 and y3
subplot(122)
plot(y(1,:),y(3,:),'k','Linewidth',2)
ylabel('y_{3t}','fontsize',24);
xlabel('y_{1t}','fontsize',24);


fprintf('Now estimating the states...')
fprintf('\n')

%Initial state vector and covariance matrix
x0 = [y(1,1); y(2,1); y(3,1)]; %Initial guesses based on observable values

P0 = diag([repmat(1,3,1)]); %Initial 3x3 covariance matrix
par0 = truepar; %Set parameters to actual values
Q.cov=diag(diag([repmat(10e-4,Nx)]));  % Add a little bit of process noise to help speed convergence  
R = []; %Initialize covariance as struct
xhat = zeros(Nx,T); %Initialize matrix of state estimates with same dimensions as x
xhat(:,1) = x0; %First state estimate is first guess

%Set up struct of information needed to run the UKF
InfDS.spkfParams = [alpha beta kappa]; %UKf scaling parameters
InfDS.ffun = @LordynO; %The transition function F
InfDS.hfun = @LormeasO; %The measurement function H
InfDS.par = par0; %True values for model parameters
InfDS.obsdim = Ny; %Number of observables
InfDS.Nsubj = Nsubj;
InfDS.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS.NxNoPar = 3;	          %Number of state variables
InfDS.partQflag = 0;	% A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS.Xdim = 3; %x dimension 3 for 3 states

R.cov = diag(truepar(4:6)); %Set R covariance matrix using known parameter values
Q.cov=diag(diag([repmat(10e-4,Nx)]));
starttime = cputime;
[xhat, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall, Pyall, Px_all]= ukfaddmultiP(x0,P0,Q, R, y,[],[],InfDS); %Use UKF to do state estimates
time_ukf = cputime-starttime %Time that UKF took
xhat; %Print the state estimates

end

%Plot the state estimates
figure(4) %relationship between x1 and x2 state estimates
subplot(121)
plot(xhat(1,:),xhat(2,:),'k','Linewidth',2)
ylabel('xhat_{2t}','fontsize',24);
xlabel('xhat_{1t}','fontsize',24);

figure(4) %relationship between x1 and x3 state estimates
subplot(122)
plot(xhat(1,:),xhat(3,:),'k','Linewidth',2)
ylabel('xhat_{3t}','fontsize',24);
xlabel('xhat_{1t}','fontsize',24);


figure(5) %Plots overlay of predicted and actual state 1
p1= plot(1:T,x(1,:),'b+','Linewidth',2);hold on;
p2 = plot(1:T, xhat(1,:), 'm', 'Linewidth', 2); hold off
ylabel('x1 state','fontsize',15);
xlabel('Time','fontsize',15);
legend('x_{1t}','xhat_{1t}');

figure(6) %Plots overlay of predicted and actual state 2
p1= plot(1:T,x(2,:),'b+','Linewidth',2);hold on;
p2 = plot(1:T, xhat(2,:), 'r', 'Linewidth', 2); hold off
ylabel('x2 state','fontsize',15);
xlabel('Time','fontsize',15);
legend('x_{2t}','xhat_{2t}');

figure(7) %Plots overlay of predicted and actual state 3
p1= plot(1:T,x(3,:),'r+','Linewidth',2);hold on;
p2 = plot(1:T, xhat(3,:), 'b', 'Linewidth', 2); hold off
ylabel('x3 state','fontsize',15);
xlabel('Time','fontsize',15);
legend('x_{3t}','xhat_{3t}');


figure(8) %Plots overlay of predicted and actual states for all 3 at once
p1=plot(1:T,x(1,:),'ko','Linewidth',2);hold on;
p2=plot(1:T,x(2,:),'r-.','Linewidth',2);
p3=plot(1:T,x(3,:),'b+');
p4 = plot(1:T, xhat(1,:), 'k', 'Linewidth', 2);
p5 = plot(1:T, xhat(2,:), 'r', 'Linewidth', 2);
p6 = plot(1:T, xhat(3,:), 'b', 'Linewidth', 2); hold off
ylabel('Data','fontsize',15);
xlabel('Time','fontsize',15);
legend('x_{1t}','x_{2t}','x_{3t}', 'xhat_{1t}', 'xhat_{2t}', 'xhat_{3t}');


%Calculate error as difference between data and estimates
x_error = x - xhat; %use matrix subtraction to get error
error_norm = vecnorm(x_error); %take columnwise norm

figure(9) %Plot norm of the error over time
p1 = plot(1:T, error_norm(1,:), 'b', 'Linewidth', 2);
ylabel('Error Norm', 'fontsize', 15);
xlabel('Time', 'fontsize', 15);
