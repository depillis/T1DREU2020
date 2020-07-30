%THIS IS THE UKF FUNCTION FOR THE PARAMETERS IN THE DUAL SETUP


function [param_hat,P_param_hat] = UKF_Params(param0, P_param_0,Q,R, obsAll, U1, U2, InferenceDS, Xprev, numLatent, time_diff)

%------SET UP VALUES TAKEN FROM INFO STRUCT---------%
hfun = InferenceDS.hfun;
alpha = InferenceDS.spkfParams(1);
beta = InferenceDS.spkfParams(2);
kappa = InferenceDS.spkfParams(3);
par = InferenceDS.par;
Odim = InferenceDS.obsdim;
Xdim = InferenceDS.Xdim;

partQflag = InferenceDS.partQflag;
Q = Q.cov;
R = R.cov;

%------------------------------------%

fig_2_3_4_5_Parameters;
wave = wave_basal; %set wave

Vdim  = size(Q,2);                                    % extract process noise dimension
Ndim  = size(R,2);                                   % extract observation noise dimension


if (nargin ~= 11) error(' [ ukf ] Number of input arguments does not match.'); end



L = Xdim;                                                % number of parameters
nsp = 2*L+1;                                              % number of sigma-points
kappa = alpha^2*(L+kappa)-L;                              % compound scaling parameter

W = [kappa 0.5 0]/(L+kappa);                              % sigma-point weights
W(3) = W(1) + (1-alpha^2) + beta;

Sqrt_L_plus_kappa = sqrt(L+kappa);



Sv = [];
if (~isempty(Q) & partQflag==1) Sv=Q;end
if (~isempty(Q) & partQflag==0),
[Sv h] = chol(Q);
    if h==0 Sv= Sv'; 
    else Sv=zeros(Vdim,Vdim);end;
end


[Sn h] = chol(R);
if h==0 Sn = Sn';else Sn=zeros(Ndim,Ndim);end
%--------------------------------------- Loop over all input vectors --------------------------------------------
for j=1:1,

obs  = obsAll((1+(j-1)*Odim):(Odim*j),:);
xh   = zeros(Xdim,1);
xh_  = zeros(Xdim,1);
yh_  = zeros(Odim,1);
inov = zeros(Odim,1);
state=param0((1+(j-1)*Xdim):(Xdim*j));

Px_ = P_param_0; %Initial guess of parameter covariance
sSzV=zeros(Xdim,nsp);
if ~isempty(Q)
  Px_ = Px_ + Q;
  S2 = [Sv];
    %---Calculate predicted observation and covariance
  sSv2 = Sqrt_L_plus_kappa*Sv;
  sSzV2 = [sSv2 -sSv2];
  size(sSzV(:,2:nsp))
  size(sSzV2)
  sSzV(:,2:nsp) = sSzV(:,2:nsp) + sSzV2;
end;
    

    
[Sx h] = chol(Px_);
if h==0 Sx=Sx'; 
else
   Sx = diag([200 0.05 0.005 0.05 10e-10 10e-10 10e-12 10e-12 10e-4 10e-3 10 2 2 .02 5 5 10 175 10^-3 ...
                      10^-3 10^-3 10^-3 10^-3 175 10e-11 10e-11 10e-4 1 4 10e-3 10^-5 10^-5 10^-5 10^-3 ...
                      200 1 1 1 5 15 .001]);
end


%------------------------------------------------------
    % Projection Step      %
Z    = cvecrep([state], nsp);            % copy needed for possible angle components section
Sz  = [Sx];
sSz  = Sqrt_L_plus_kappa * Sz;
sSzM = [sSz -sSz];
Z(:,2:nsp) = Z(:,2:nsp) + sSzM;           % build sigma-point set
Z(:,2:nsp) = Z(:,2:nsp) + sSzV(:,2:nsp);
    %-- Calculate predicted state mean and covariance

X_ = param0;
xh_(:,1) = X_;

temp1 = Z - cvecrep(xh_(:,1),nsp);
tspan = [0 7 * time_diff];
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
X_newParam = zeros(numLatent, nsp);


for z = 1:nsp
    if (time_diff == 0)
        X_temp = Xprev';
    else
        [~, solutions] = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, Z(:,z)), tspan, Xprev, options); %solve ode with ode45    
        X_temp = solutions(end,:);
    end
    X_newParam(:,z) = X_temp';
end
X_bps = X_newParam;


Y_ = feval(hfun, X_bps, [], [],1, par); %propagate through observation model
yh_(:,1) = W(1)*Y_(:,1) + W(2)*sum(Y_(:,2:nsp),2);
temp2 = Y_ - cvecrep(yh_(:,1),nsp);

Py  = W(3)*temp2(:,1)*temp2(:,1)' + W(2)*temp2(:,2:nsp)*temp2(:,2:nsp)';
  if ~isempty(R)
      Py = Py + R;end
Pxy = W(3)*temp1(:,1)*temp2(:,1)' + W(2)*temp1(:,2:nsp)*temp2(:,2:nsp)';

    %------------------------------------------------------
    % Update Step



KG = Pxy * inv(Py); %calculate kalman gain
inov(:,1) = obs(:,1) - yh_(:,1); %calculate innovation (difference between observable and prior estimate)
xh(:,1) = xh_(:,1) + KG*inov(:,1); %calculate posterior estimate
    
Px = Px_ - KG*Py*KG'; %update covariance matrix
   
state = xh(:,1);
Pstate = Px;

for l = 1:Xdim
    if (state(l,1) < 0)  %make sure no parameter goes negative
        state(l,1) = param0(l,1);
    end
end

param_hat = state;
P_param_hat = Pstate;
    %----------------------------------------------------------------------
   
end   



end





