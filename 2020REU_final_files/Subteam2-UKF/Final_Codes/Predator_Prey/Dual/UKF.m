%THE STATES UKF FUNCTION FOR PREDATOR PREY MODEL DUAL UKF


function [xhat,Phat] = UKF(stateAll, Pstate,Q,R, obsAll, U1, U2, InferenceDS)
%extract information from the inputs
P_orig = Pstate;
ffun = InferenceDS.ffun;
hfun = InferenceDS.hfun;
alpha = InferenceDS.spkfParams(1);
beta = InferenceDS.spkfParams(2);
kappa = InferenceDS.spkfParams(3);
par = InferenceDS.par;
Odim = InferenceDS.obsdim;
Xdim = InferenceDS.Xdim;
Nsubj = InferenceDS.Nsubj;
partQflag = InferenceDS.partQflag;
Q = Q.cov;
R = R.cov;

Vdim  = size(Q,2);                                    % extract process noise dimension
Ndim  = size(R,2);                                   % extract observation noise dimension




if (nargin ~= 8) error(' [ ukf ] Number of input arguments does not match.'); end %check number of inputs


L = Xdim;                                                % number of states
nsp = 2*L+1;                                              % number of sigma-points
kappa = alpha^2*(L+kappa)-L;                              % compound scaling parameter

W = [kappa 0.5 0]/(L+kappa);                              % sigma-point weights
                                                          %[Initial sigma point mean weight,
                                                          % other sigma points mean and covariance
                                                          % weight,
                                                          %Initial sigma point covaraince weight]
W(3) = W(1) + (1-alpha^2) + beta;

Sqrt_L_plus_kappa = sqrt(L+kappa);


%Cholesky Factorication of Q - Process noise covariance
Sv = [];
if (~isempty(Q) & partQflag==1) Sv=Q;end
if (~isempty(Q) & partQflag==0),
[Sv h] = chol(Q);
    if h==0 Sv= Sv'; 
    else Sv=zeros(Vdim,Vdim);end;
end

%Cholesky factorization of R - Measurement noise covariance
[Sn h] = chol(R);
if h==0 Sn = Sn';else Sn=zeros(Ndim,Ndim);end
%--------------------------------------- Loop over all input vectors --------------------------------------------
for j=1:Nsubj,
obs  = obsAll((1+(j-1)*Odim):(Odim*j),:);

xh   = zeros(Xdim,1); %will hold final estimates
xh_  = zeros(Xdim,1); %will hold prior estimates
yh_  = zeros(Odim,1); %will hod prior observable estimate
inov = zeros(Odim,1); %will hold innovation
state=stateAll((1+(j-1)*Xdim):(Xdim*j)); %initial state

%Cholesky factorization of covariance matrix
[Sx h] = chol(Pstate);
if h==0 Sx=Sx'; else Sx = diag([repmat(.001,Xdim,1)]);end

%------------------------------------------------------
    % Projection Step       %
Z    = cvecrep([state], nsp);            % copy needed for possible angle components section
Sz  = [Sx];
sSz  = Sqrt_L_plus_kappa * Sz;
sSzM = [sSz -sSz];
Z(:,2:nsp) = Z(:,2:nsp) + sSzM;           % build sigma-point set
    %-- Calculate predicted state mean and covariance
  
tspan = [0 1]; %timespan to project forward with ODEs
X_ = zeros(Xdim, nsp);
options=odeset('RelTol',1e-12,'AbsTol',1e-12); %set tolerances for ODE solver


%loop through all sigma points, pass them through the transformation
%function (by solving thr ODE) and record them in a matrix
for z = 1:nsp
[~, solutions] = ode45(@(t, y) ffun(t, y, par), tspan, Z(:,z), options); %solve ode with ode45 to project forward 1 time point
X_temp = solutions(end,:);        
X_(:,z) = X_temp';
end
X_bps = X_; %set of outputs



xh_(:,1) = W(1)*X_(:,1) + W(2)*sum(X_(:,2:nsp),2); %calculate prior estimate
temp1 = X_ - cvecrep(xh_(:,1),nsp);
Px_ = W(3)*temp1(:,1)*temp1(:,1)' + W(2)*temp1(:,2:nsp)*temp1(:,2:nsp)'; %calculate prior for covariance matrix
 
%Add process noise to prior covariance matrix calculation
if ~isempty(Q)
  Px_ = Px_ + Q;
end;
    
Y_ = feval(hfun, X_bps, [], [],1, par); %propagate through observation model
yh_(:,1) = W(1)*Y_(:,1) + W(2)*sum(Y_(:,2:nsp),2); %prior observable estimate

temp2 = Y_ - cvecrep(yh_(:,1),nsp);
%Add measurement noise to covariance matrix of observations
Py  = W(3)*temp2(:,1)*temp2(:,1)' + W(2)*temp2(:,2:nsp)*temp2(:,2:nsp)'; %variance of observable
  if ~isempty(R)
      Py = Py + R;end
Pxy = W(3)*temp1(:,1)*temp2(:,1)' + W(2)*temp1(:,2:nsp)*temp2(:,2:nsp)'; %covariance between states and observables

    %------------------------------------------------------
    % MEASUREMENT UPDATE


KG = Pxy * inv(Py); %calculate kalman gain
inov(:,1) = obs(:,1) - yh_(:,1); %calculate innovation
xh(:,1) = xh_(:,1) + KG*inov(:,1); %Use the prior prediction state vector, innovation and Kalman Gain 
                                       %to calculate the posterior prediction
                                       %state vector
Px = Px_ - KG*Py*KG';%Use the prior prediction covariance matrix, observations covariance matrix, and Kalman Gain
                          %to calculate the posterior prediction for the
                          %covariance matrix
   
state = xh(:,1);
Pstate = Px;
for l = 1:Xdim
    if state(l,1) < 0 %check that states did not go negative
        state(l,1) = stateAll(l,1);
        Pstate = P_orig;
    end
end
xhat = state; %set output
Phat = Pstate;


end

