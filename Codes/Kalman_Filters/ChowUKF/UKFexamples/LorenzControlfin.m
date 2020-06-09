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

clear all;
clc;
echo off;

% INITIALISATION AND PARAMETERS:
% ==============================

no_of_runs = 1;           % number of simulation runs
T = 2000;
Ny = 3;
Nx = 6;
Nstates=3;
Npar=6;
Nsubj=1;
truepar = [10 28 8/3 26 34 32];
%Q = zeros(Nx,Nx);
Q = [];

alpha = 10e-4; % UKF : point scaling parameter
beta = 2;  % UKF : scaling parameter for higher order terms of Taylor series expansion 
kappa = 0;    % UKF : sigma point selection scaling parameter (best to leave this = 0)

%**************************************************************************************

% MAIN LOOP

for j=1:no_of_runs,

  rand('state',sum(100*clock));   % Shuffle the pack!
  randn('state',sum(100*clock));   % Shuffle the pack!
  

% GENERATE THE DATA:
% ==================

x = zeros(Nx,T);
y = zeros(Ny,T);
measureNoise = zeros(Ny,T);
measureNoise(:,1) = [sqrt(truepar(4))*randn(1,1); sqrt(truepar(5))*randn(1,1); sqrt(truepar(6))*randn(1,1)];   
x(:,1) = [-5; 20; -5;truepar(1:3)'];                         % Initial state.
y(:,1) = feval('LormeasO',x(:,1),[],[],[],[]) + measureNoise(:,1);

for t=2:T,
  measureNoise(:,t) = [sqrt(truepar(4))*randn(1,1); sqrt(truepar(5))*randn(1,1); sqrt(truepar(6))*randn(1,1)];    
  x(:,t) = feval('LordynO',x(:,t-1),[],[],t,[]);     
  y(:,t) = feval('LormeasO',x(:,t),[],[],[],[]) + measureNoise(:,t);      
end;  

% PLOT THE GENERATED DATA:
% ========================
figure(1)
p1=plot(1:T,x(1,:),'k','Linewidth',2);hold on;
p2=plot(1:T,x(2,:),'r-.','Linewidth',2);
p3=plot(1:T,x(3,:),'b+');hold off;
ylabel('Data','fontsize',15);
xlabel('Time','fontsize',15);
legend('x_{1t}','x_{2t}','x_{3t}');

figure(2)
subplot(121)
plot(x(1,:),x(2,:),'k','Linewidth',2)
ylabel('x_{2t}','fontsize',24);
xlabel('x_{1t}','fontsize',24);

figure(2)
subplot(122)
plot(x(1,:),x(3,:),'k','Linewidth',2)
ylabel('x_{3t}','fontsize',24);
xlabel('x_{1t}','fontsize',24);

figure(3)
subplot(121)
plot(y(1,:),y(2,:),'k','Linewidth',2)
ylabel('y_{2t}','fontsize',24);
xlabel('y_{1t}','fontsize',24);

figure(3)
subplot(122)
plot(y(1,:),y(3,:),'k','Linewidth',2)
ylabel('y_{3t}','fontsize',24);
xlabel('y_{1t}','fontsize',24);



%printf('\n')
%printf('Press a key to continue')  
%printf('\n')
%printf('\n')
fprintf('Now estimating the model parameters...')
fprintf('\n')

%Initial state vector and covariance matrix
x0 = [y(1,1); y(2,1); y(3,1);6;20;2];
P0 = diag([repmat(1,6,1)]);
par0=[30 40 35];
Q.cov=diag(diag([repmat(10e-4,Nx)]));  % Add a little bit of process noise to help speed convergence  
R = [];
xhat = zeros(Nx,T);
xhat(:,1) = x0;


InfDS.spkfParams = [alpha beta kappa];
InfDS.ffun = @LordynO;
InfDS.hfun = @LormeasO;
InfDS.par = par0;
InfDS.obsdim = Ny;
InfDS.Nsubj = Nsubj;
InfDS.pNoiseAdaptpar = [];    %Parameters to shape process noise to speed convergence. Set to empty matrix if don't want to shape process noise.
InfDS.NxNoPar = 3;	          %Number of state variables
InfDS.partQflag = 0;% A flag--if only estimating part of the state covariance matrix, set to 1; otherwise set to 0 
InfDS.Xdim = 6;

starttime = cputime;
options=optimset('Display','final','TolX',1e-7,'MaxIter',10000,'MaxFunEvals',10000);
[finalpar,fval,exitflag,output] =fminsearch('Lor_loglik',par0,options,y,x0,P0,InfDS,Q, R);
time_ML = cputime-starttime


R.cov = diag(finalpar);
Q.cov=diag(diag([repmat(10e-4,Nx)]));
starttime = cputime;
[xhat, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall, Pyall, Px_all]= ukfaddmultiP(x0,P0,Q, R, y,[],[],InfDS); 
time_ukf = cputime-starttime


hess = hessian('Lor_loglik',finalpar',y,x0,Px,InfDS, Q, R);


SEpar = sqrt(diag(inv(hess)));
finalpar2 = finalpar';
Pxall = diag(Px(:));


%filename=['C:\D\ParticleFilter\upf_demos\Lorenz\LorT2000'];
%filename2 = [sprintf('%s',filename),sprintf('%s','par.txt')];
%filename3 = [sprintf('%s',filename),sprintf('%s','states.txt')];
%filename4 = [sprintf('%s',filename),sprintf('%s','Px.txt')];
%filename5 = [sprintf('%s',filename),sprintf('%s','SEpar.txt')];
%filename6 = [sprintf('%s',filename),sprintf('%s','Truestates.txt')];
%filename7 = [sprintf('%s',filename),sprintf('%s','Obsy.txt')];
%filename8 = [sprintf('%s',filename),sprintf('%s','timeML.txt')];
%filename9 = [sprintf('%s',filename),sprintf('%s','timeukf.txt')];



%dlmwrite(filename2,finalpar,'\t');
%dlmwrite(filename3,xhat,'\t');
%dlmwrite(filename4,Pxall,'\t');
%dlmwrite(filename5,SEpar,'\t');
%dlmwrite(filename6,x,'\t');
%dlmwrite(filename7,y,'\t');
%dlmwrite(filename8,time_ML,'\t');
%dlmwrite(filename9,time_ukf,'\t');
end

figure(4)
subplot(121)
plot(xhat(1,:),xhat(2,:),'k','Linewidth',2)
ylabel('xhat_{2t}','fontsize',24);
xlabel('xhat_{1t}','fontsize',24);

figure(4)
subplot(122)
plot(xhat(1,:),xhat(3,:),'k','Linewidth',2)
ylabel('xhat_{3t}','fontsize',24);
xlabel('xhat_{1t}','fontsize',24);
