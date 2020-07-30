% PURPOSE : Fitting the Lorenz system model to noisy multiple-subject time series data 

% AUTHOR: Sy-Miin Chow (schow@nd.edu)
% Portions of this script was adapted from previous scripts written by:
%       Nando de Freitas      (jfgf@cs.berkeley.edu)
%       Rudolph van der Merwe (rvdmerwe@ece.ogi.edu)
%       Kevin Murphy (?)
% So please acknowledge their previous work
% DATE     : 18 July 2005

clear all;
clc;
echo off;
%path('../ukf',path);

% INITIALISATION AND PARAMETERS:
% ==============================

no_of_runs =10;           % number of simulation runs
T = 300;
Ny = 3;
Nx = 6;
Nstates=3;
Nsubj=100;
truepar = [10 28 8/3 26 34 32];
ukf_stateEst = zeros(no_of_runs ,Nx+1);
ukf_PxEst = zeros(no_of_runs ,Nx+1);
temp_ukf = zeros(Nx,T);
temp_par = zeros(1,3);
temp_SE = zeros(1,3);
Index_ukf = [0 0 0 0 0 0];
Index_par = [0 0 0];
Index_SE = [0 0 0];
temp_ukfPx = zeros(Nx,T);
time_ukf = zeros(no_of_runs,1);
Xh_ukfPxFin = zeros(Nx,T);
finalpar = zeros(no_of_runs,3);
time_ML = zeros(no_of_runs,1);
exiflag=zeros(no_of_runs,1);
SEpar = zeros(no_of_runs,3);
Pxout = zeros(no_of_runs,Nx);


alpha = 10e-4; % UKF : point scaling parameter
beta = 2;  % UKF : scaling parameter for higher order terms of Taylor series expansion 
kappa = 0;    % UKF : sigma point selection scaling parameter (best to leave this = 0)

%**************************************************************************************

% MAIN LOOP

for k=1:no_of_runs,

  rand('state',sum(100*clock));   % Shuffle the pack!

% GENERATE THE DATA:
% ==================

x = zeros(Nx*Nsubj,T);
y = zeros(Ny*Nsubj,T);
processNoise = zeros(Nx,T);
measureNoise = zeros(Ny,T);

for j=1:Nsubj,
measureNoise(:,1) = [sqrt(truepar(4))*randn(1,1); sqrt(truepar(5))*randn(1,1); sqrt(truepar(6))*randn(1,1)];  
x(1+(j-1)*Nx:(j*Nx),1) = [-5+ sqrt(25)*randn(1,1) ; 20+sqrt(25)*randn(1,1);-5+sqrt(25)*randn(1,1);truepar(1:3)'];  % Give everyone slightly different initial states
%x(1+(j-1)*Nx:(j*Nx),1) = [-5+ sqrt(1)*randn(1,1) ; 20+sqrt(1)*randn(1,1);-5+sqrt(1)*randn(1,1);truepar(1:3)'];  % Give everyone slightly different initial states

y(1+(j-1)*Ny:(j*Ny),1) = feval('LormeasO',x(1+(j-1)*Ny:(j*Ny),1),[],[],[],truepar) + measureNoise(:,1);
end


for j=1:Nsubj,
for t=2:T,
 measureNoise(:,t) = [sqrt(truepar(4))*randn(1,1); sqrt(truepar(5))*randn(1,1); sqrt(truepar(6))*randn(1,1)];  
 x(1+(j-1)*Nx:(j*Nx),t) = feval('LordynO',x(1+(j-1)*Nx:(j*Nx),t-1),[],[],t,truepar);     
  y(1+(j-1)*Ny:(j*Ny),t) = feval('LormeasO',x(1+(j-1)*Nx:(j*Nx),t),[],[],[],truepar)+measureNoise(:,t);      
end;  
end;


xavg = zeros(Nx,T);
yavg = zeros(Ny,T);
for j=1:Nsubj,
measureNoise(:,1) = [sqrt(truepar(4))*randn(1,1); sqrt(truepar(5))*randn(1,1); sqrt(truepar(6))*randn(1,1)];  
xavg(:,1) = [-5 ; 20;-5;truepar(1:3)'];  % Give everyone slightly different initial states
yavg(:,1) = feval('LormeasO',x(:,1),[],[],[],truepar) + measureNoise(:,1);
end

for t=2:T,
measureNoise(:,t) = [sqrt(truepar(4))*randn(1,1); sqrt(truepar(5))*randn(1,1); sqrt(truepar(6))*randn(1,1)];  
  xavg(:,t) = feval('LordynO',x(1+(j-1)*Nx:(j*Nx),t-1),[],[],t,truepar)+processNoise(:,t);     
  yavg(:,t) = feval('LormeasO',x(1+(j-1)*Nx:(j*Nx),t),[],[],[],truepar)+measureNoise(:,t);      
end;

% PLOT THE GENERATED DATA:
% ========================

figure(1)
clf
subplot(121)
for j=1:Nsubj,
x1 = x(1+(j-1)*Nx:j*Nx,:);
p1=plot(1:T,x1(1,:),'k-');hold on;
end
p2=plot(1:T,xavg(1,:),'w-*','LineWidth',3);hold off;

ylabel('x_{1t}','fontsize',24);
xlabel('Time','fontsize',24);

figure(1)
subplot(122)
for j=1:Nsubj,
ytemp1 = y(1+(j-1)*Ny,:);
ytemp2 = y(2+(j-1)*Ny,:);
ytemp3 = y(j*Ny,:);

p1=plot(1:T,ytemp1,'k-');hold on;
end
p2=plot(1:T,xavg(1,:),'w-*','LineWidth',3);hold off;

ylabel('Noisy x_{1t}','fontsize',24);
xlabel('Time','fontsize',24);


fprintf('\n')
fprintf('\n')
fprintf('Now estimating the model parameters...')
fprintf('\n')



figure(2)
subplot(221)
for j=1:Nsubj,
x1 = x(1+(j-1)*Nx:j*Nx,:);
p1=plot(x1(1,:),x1(2,:),'k-');hold on;
end
p2=plot(xavg(1,:),xavg(2,:),'w-o','LineWidth',2);hold off;
ylabel('x_{2t}','fontsize',24);
xlabel('x_{1t}','fontsize',24);


figure(2)
subplot(222)
for j=1:Nsubj,
x1 = x(1+(j-1)*Nx:j*Nx,:);
p1=plot(x1(1,:),x1(3,:),'k-');hold on;
end
p2=plot(xavg(1,:),xavg(3,:),'w-o','LineWidth',2);hold off;
ylabel('x_{3t}','fontsize',24);
xlabel('x_{1t}','fontsize',24);

figure(2)
subplot(223)
for j=1:Nsubj,
y1 = y(1+(j-1)*Ny:j*Ny,:);
p1=plot(y1(1,:),y1(2,:),'k-');hold on;
end
p2=plot(yavg(1,:),yavg(2,:),'w-o','LineWidth',2);hold off;
ylabel('y_{2t}','fontsize',24);
xlabel('y_{1t}','fontsize',24);


figure(2)
subplot(224)
for j=1:Nsubj,
y1 = y(1+(j-1)*Ny:j*Ny,:);
p1=plot(y1(1,:),y1(3,:),'k-');hold on;
end
p2=plot(yavg(1,:),yavg(3,:),'w-o','LineWidth',2);hold off;
ylabel('y_{3t}','fontsize',24);
xlabel('y_{1t}','fontsize',24);


x00 = zeros(Nx*Nsubj,1);
P00 = diag([repmat(1,6,1)]);
par0=[30 40 35];
for j=1:Nsubj,
x00(1+(j-1)*Nx:j*Nx) = [y(1+(j-1)*Ny,1); y(2+(j-1)*Ny,1); y(3+(j-1)*Ny,1);6;20;2];
end


Arg.type = 'gmm';              % GMM noise source
Arg.cov_type = 'full';         % GSPF use square-root covariance matrices
Arg.dim = Nx;                    % process noise dimension
Arg.M = 1;                     % single coponent
Arg.weights = [1];             % component weight
Arg.mu = repmat(0,Nx,1);
Arg.cov = diag(diag([repmat(10e-2,Arg.dim)]));
pNoise = gennoiseds(Arg);      % generate process noise data structure
pNoise.adaptMethod = [];                     % setup PE process noise adaptation method
pNoise.adaptParams = []; 


Arg2.type = 'gmm';              % GMM noise source
Arg2.cov_type = 'full';         % GSPF use square-root covariance matrices
Arg2.dim = Ny;          % process noise dimension
Arg2.M = 1;                     % single coponent
Arg2.weights = [1];             % component weight
Arg2.mu = repmat(0,Ny,1);
Arg2.cov = diag(truepar(4:6));
oNoise = gennoiseds(Arg2);      % generate process noise data structure


InfDS2.statedim=Nx;
InfDS2.obsdim=Ny;
InfDS2.U1dim=0;
InfDS2.U2dim=0;
InfDS2.Xdim = Nx;
InfDS2.Vdim=Nx;
InfDS2.Ndim=Ny;
InfDS2.spkfParams  = [10e-4  2 0];
InfDS2.ffun=@LordynO;
InfDS2.hfun = @LormeasO; 
InfDS2.estimateType='GMMmean';
InfDS2.likelihood=@directLorlik;
InfDS2.par=par0;
InfDS2.Nsubj = Nsubj;
InfDS2.pNoiseAdaptpar=[];
InfDS2.NxNoPar = 3;
InfDS2.Nsubj2 = 1;
InfDS2.partQflag = 0;

%Perform ukf estimation with ML
%=====================================


starttime = cputime;
options=optimset('Display','final','TolX',1e-5,'MaxIter',1000,'MaxFunEvals',10000);
[finalpar(k,:),fval,exitflag(k),output] =fminsearch('Lor_loglikMultiP2',par0,options,y,x00,P00,InfDS2,pNoise,oNoise);
time_ML(k) = cputime-starttime;


%Pxall = zeros(Nx,Nx,T);
%Pxall(:,:,1) = P00;
oNoise.cov=diag(finalpar(k,:));
%Xh_ukf=zeros(Nx,T);
%Xh_ukf(:,1) = xavg(:,1);
starttime = cputime;
[Xh_ukf, Px_ukf, xh_, Px_, yh_, inov, Py, KG, lik,Pxall]= ukfaddmultiPtemp(x00,P00,pNoise, oNoise, y,[],[],InfDS2,0,1,x); 



ukf_stateEst(k,1) = k;
ukf_stateEst(k,2:end) = Xh_ukf(:,end)';
ukf_PxEst(k,1) =k;
ukf_PxEst(k,2:end) = diag(Px_ukf)';


for j=1:Nx,
    if ~isempty(Xh_ukf(j,end))
    temp_ukf(j,:) = temp_ukf(j,:)+Xh_ukf(j,:);
    Xh_ukf(j,end)
    Index_ukf(j) = Index_ukf(j)+1;end
end

hess = hessian('Lor_loglikMultiP2',finalpar(k,:)',y,x00,P00,InfDS2,pNoise,oNoise);

SEpar(k,:) = abs(sqrt(diag(inv(hess))))';

for j=1:3,
    if ~isempty(finalpar(j))
    temp_par(j) = temp_par(j)+finalpar(k,j);
    Index_par(j) = Index_par(j)+1;end
    if ~isempty(SEpar(k,j))
    temp_SE(j) = temp_SE(j) + SEpar(k,j);
    Index_SE(j) = Index_SE(j)+1;end
end

Pxhold_ukf = zeros(Nx,T);
for t=1:T,
    Pxhold_ukf(:,t) = diag(Pxall(:,:,t));
end
Xh_ukfPxFin = Xh_ukfPxFin+Pxhold_ukf;
Pxout(k,:) = diag(Px_ukf)';

  fprintf('ML : t = ',k);k
  fprintf('\n');
end %end of no_of_runs


Index_ukf2=repmat(Index_ukf',1,T);
Index_par=repmat(Index_par',1);
Xh_ukfFin = temp_ukf./Index_ukf2;

Xh_parFin = temp_par./Index_par';
Xh_SEFin = temp_SE./Index_SE;
Xh_ukfPxFin = Xh_ukfPxFin./Index_ukf2;

filename=['C:\D\ParticleFilter\upf_demos\LorenzSum\MultiP\Final50runs\Subs\PartT300N100runs20'];
filename2 = [sprintf('%s',filename),sprintf('%s','ukfstate.txt')];
filename3 = [sprintf('%s',filename),sprintf('%s','MLstate.txt')];
filename4 = [sprintf('%s',filename),sprintf('%s','AllukfPx.txt')];
filename5 = [sprintf('%s',filename),sprintf('%s','AllMLSE.txt')];
filename6 = [sprintf('%s',filename),sprintf('%s','timeML.txt')];
filename7 = [sprintf('%s',filename),sprintf('%s','exitflag.txt')];
filename8 = [sprintf('%s',filename),sprintf('%s','Allukfstate.txt')];
filename9 = [sprintf('%s',filename),sprintf('%s','AllMLstate.txt')];

%dlmwrite(filename2,Xh_ukfFin,'\t');
%dlmwrite(filename3,Xh_parFin,'\t');
%dlmwrite(filename4,Pxout,'\t');
%dlmwrite(filename5,SEpar,'\t');
%dlmwrite(filename6,time_ML,'\t');
%dlmwrite(filename7,exitflag,'\t');
%dlmwrite(filename8,ukf_stateEst,'\t');
%dlmwrite(filename9,finalpar,'\t');



%Xh_ukfFin = importdata('C:\D\ParticleFilter\upf_demos\LorenzSum\MultiP\50runs\T1000N20ukfstate.txt');
%ML_ukfFin = importdata('C:\D\ParticleFilter\upf_demos\LorenzSum\MultiP\50runs\MultiPT200N200MLstate.txt');




figure(4)
subplot(121)
p3=plot(1:T,y((1+(Nsubj-1)*Ny),:),'g*','Linewidth',2);hold on;
p1=plot(1:T,x((1+(Nsubj-1)*Nx),:),'k','Linewidth',7);
p2=plot(1:T,Xh_ukfFin(1,:),'m-.','Linewidth',2);hold off;
ylabel('x_{1t}','fontsize',24);
xlabel('Time','fontsize',24);
legend([p1 p2 p3],'True x1','ukf x1','noisy x1');
title('UKF','fontsize',24);


figure(4)
subplot(122)
p3=plot(1:T,y((3+(Nsubj-1)*Ny),:),'g*','Linewidth',2);hold on;
p1=plot(1:T,x((3+(Nsubj-1)*Nx),:),'k','Linewidth',7);
p2=plot(1:T,Xh_ukfFin(3,:),'m-.','Linewidth',2);hold off;
ylabel('x_{3t}','fontsize',24);
xlabel('Time','fontsize',24);
legend([p1 p2 p3],'True x3','ukf x3','noisy x3');
title('UKF','fontsize',24);


figure(1)
subplot(221)
p1=plot(1:T,x((4+(Nsubj-1)*Nx),:),'k','Linewidth',2);hold on;
p2=plot(1:T,Xh_ukfFin(4,:),'b-.','Linewidth',2);hold off;
ylabel('Sigma','fontsize',24);
xlabel('Time','fontsize',24);
legend([p1 p2 p3],'True sigma','ukf sigma');

figure(1)
subplot(222)
p1=plot(1:T,x((5+(Nsubj-1)*Nx),:),'k','Linewidth',2);hold on;
p2=plot(1:T,Xh_ukfFin(5,:),'b-.','Linewidth',2);hold off;
ylabel('rho','fontsize',24);
xlabel('Time','fontsize',24);
legend([p1 p2],'True rho','ukf rho');

figure(1)
subplot(223)
p1=plot(1:T,x((6+(Nsubj-1)*Nx),:),'k','Linewidth',2);hold on;
p2=plot(1:T,Xh_ukfFin(6,:),'b-.','Linewidth',2);hold off;
ylabel('Gamma','fontsize',24);
xlabel('Time','fontsize',24);
legend([p1 p2],'True gamma','ukf gamma');


figure(2)
subplot(221)
p2=plot(1:T,Xh_ukfPxFin(1,:),'b-.','Linewidth',2);
ylabel('Variance x_{1t}','fontsize',24);
xlabel('Time','fontsize',24);

figure(2)
subplot(222)
p2=plot(1:T,Xh_ukfPxFin(3,:),'b-.','Linewidth',2);
ylabel('Variance x_{3t}','fontsize',24);
xlabel('Time','fontsize',24);

figure(2)
subplot(223)
p2=plot(1:T,Xh_ukfPxFin(5,:),'b-.','Linewidth',2);
ylabel('Variance rho','fontsize',24);
xlabel('Time','fontsize',24);

figure(2)
subplot(224)
p2=plot(1:T,Xh_ukfPxFin(6,:),'b-.','Linewidth',2);
ylabel('Variance beta','fontsize',24);
xlabel('Time','fontsize',24);
