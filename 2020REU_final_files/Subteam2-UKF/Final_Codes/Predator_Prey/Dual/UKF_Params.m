%FUNCTION FOR PARAMETER UKF PART OF DUAL UKF


function [param_hat,P_param_hat] = UKF_Params(param0, P_param_0,Q,R, obsAll, U1, U2, InferenceDS, Xprev, numLatent)
%extract information from input
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



if (nargin ~= 10) error(' [ ukf ] Number of input arguments does not match.'); end %check number of input arguments


L = Xdim;                                                %number of parameters
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

obs  = obsAll(1:2,:); %get observation
xh   = zeros(Xdim,1); %will hold posterior estimate
xh_  = zeros(Xdim,1); %will hold prior estimate
yh_  = zeros(Odim,1); %will hold prior estimate of observable
inov = zeros(Odim,1); %will hold innovation
state=param0((1+(j-1)*Xdim):(Xdim*j));

Px_ = P_param_0; %set initial covariance matrix of parameters
sSzV=zeros(Xdim,nsp);
if ~isempty(Q)
  Px_ = Px_ + Q; %update covariance using process noise
 
    %---Calculate predicted observation and covariance
  sSv2 = Sqrt_L_plus_kappa*Sv;
  sSzV2 = [sSv2 -sSv2];
  sSzV(:,2:nsp) = sSzV(:,2:nsp) + sSzV2;
end;
    
%Cholesky factorization of covariance matrix
[Sx h] = chol(Px_);
if h==0 Sx=Sx'; 
else
    Sx = diag([repmat(.01,Xdim,1)]);
end

%------------------------------------------------------
    % Projection Step       %
Z    = cvecrep([state], nsp);            % copy needed for possible angle components section
Sz  = [Sx];
sSz  = Sqrt_L_plus_kappa * Sz;
sSzM = [sSz -sSz];
Z(:,2:nsp) = Z(:,2:nsp) + sSzM;           % build sigma-point set
Z(:,2:nsp) = Z(:,2:nsp) + sSzV(:,2:nsp);
    %-- Calculate predicted state mean and covariance

X_ = param0; %prior estimate is just initial parameters because transition matrix is identity
xh_(:,1) = X_;

temp1 = Z - cvecrep(xh_(:,1),nsp);
tspan = [0 1]; %time span for solving ODE
options=odeset('RelTol',1e-12,'AbsTol',1e-12); %set tolerances for ode sovler
X_newParam = zeros(numLatent, nsp); %will hold output from ode solver

%Loop over all sigma vectors (parameters) and pass each one into the ODE
%solver with the same initial condition
for z = 1:nsp
    [~, solutions] = ode45(@(t,y) Lotka_Volterra_Model (t, y, Z(:,z)), tspan, Xprev, options); %solve ode for each set of parameters in sigma point set
    X_temp = solutions(end,:);
    X_newParam(:,z) = X_temp';
end
X_bps = X_newParam;

Y_ = feval(hfun, X_bps, [], [],1, par); %propagate through observation model
yh_(:,1) = W(1)*Y_(:,1) + W(2)*sum(Y_(:,2:nsp),2); %calculate prior observation estimate

temp2 = Y_ - cvecrep(yh_(:,1),nsp);
Py  = W(3)*temp2(:,1)*temp2(:,1)' + W(2)*temp2(:,2:nsp)*temp2(:,2:nsp)';
  if ~isempty(R)
      Py = Py + R;end
Pxy = W(3)*temp1(:,1)*temp2(:,1)' + W(2)*temp1(:,2:nsp)*temp2(:,2:nsp)'; %covariance between parameters and observable

    %------------------------------------------------------
    % MEASUREMENT UPDATE

KG = Pxy * inv(Py); %calculate kalman gain
inov(:,1) = obs(:,1) - yh_(:,1); %calculate innovation
xh(:,1) = xh_(:,1) + KG*inov(:,1); %calculate posterior estimate of parameters
Px = Px_ - KG*Py*KG'; %posterior estimate of parameter covaraince
   
state = xh(:,1);
Pstate = Px;
for l = 1:Xdim
    if (state(l,1) < 0) || (state(l,1) > 2) %make sure parameters stay in biological range
        state(l,1) = param0(l,1);
    end
end
param_hat = state; %set output
P_param_hat = Pstate;
    %----------------------------------------------------------------------
   
end   



end





