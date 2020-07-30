%THIS IS THE UKF FUNCTION FOR STATES IN THE DUAL SETUP

function [xhat,Phat] = UKF(stateAll, Pstate,Q,R, obsAll, U1, U2, InferenceDS, time_diff)
%EXTRACT INFORMATION FOR STRUCT
P_orig = Pstate;
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

fig_2_3_4_5_Parameters; %load parameters to get wave value
wave = wave_basal; %set wave
Vdim  = size(Q,2);                                    % extract process noise dimension
Ndim  = size(R,2);                                   % extract observation noise dimension

if (nargin ~= 9) error(' [ ukf ] Number of input arguments does not match.'); end



L = Xdim;                                                % number of states
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
state=stateAll((1+(j-1)*Xdim):(Xdim*j));

    

[Sx h] = chol(Pstate);
if h==0 Sx=Sx'; else Sx = diag([repmat(.001,Xdim,1)]);end

%------------------------------------------------------
    % Projection Step      %
Z    = cvecrep([state], nsp);            % copy needed for possible angle components section
Sz  = [Sx];
sSz  = Sqrt_L_plus_kappa * Sz;
sSzM = [sSz -sSz];
Z(:,2:nsp) = Z(:,2:nsp) + sSzM;           % build sigma-point set
    %-- Calculate predicted state mean and covariance
  
tspan = [0 7 * time_diff];
X_ = zeros(Xdim, nsp);
options=odeset('RelTol',1e-12,'AbsTol',1e-12);

for z = 1:nsp
    if (time_diff == 0) %if no time difference don't solve ODE
        X_temp = Z(:,z);
    else
        [~, solutions] = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, par), tspan, Z(:,z), options); %solve ode with ode15s    
        X_temp = solutions(end,:);   
    end
    X_(:,z) = X_temp';
end
X_bps = X_;



xh_(:,1) = W(1)*X_(:,1) + W(2)*sum(X_(:,2:nsp),2);
temp1 = X_ - cvecrep(xh_(:,1),nsp);
Px_ = W(3)*temp1(:,1)*temp1(:,1)' + W(2)*temp1(:,2:nsp)*temp1(:,2:nsp)';
sSzV=zeros(Xdim,nsp);
if ~isempty(Q)
  Px_ = Px_ + Q;
  S2 = [Sv];
    %---Calculate predicted observation and covariance
  sSv2 = Sqrt_L_plus_kappa*Sv;
  sSzV2 = [sSv2 -sSv2];
  sSzV(:,2:nsp) = sSzV(:,2:nsp) + sSzV2;
end;
    

Y_ = feval(hfun, X_bps, [], [],1, par); %propagate through observation model
yh_(:,1) = W(1)*Y_(:,1) + W(2)*sum(Y_(:,2:nsp),2); %guess for observable

temp2 = Y_ - cvecrep(yh_(:,1),nsp);
Py  = W(3)*temp2(:,1)*temp2(:,1)' + W(2)*temp2(:,2:nsp)*temp2(:,2:nsp)';
  if ~isempty(R)
      Py = Py + R;end
Pxy = W(3)*temp1(:,1)*temp2(:,1)' + W(2)*temp1(:,2:nsp)*temp2(:,2:nsp)';

    %------------------------------------------------------
    % Update Step


KG = Pxy * inv(Py); %calculate kalman gain
inov(:,1) = obs(:,1) - yh_(:,1); %calculate innovation (difference between observable and estimate)
xh(:,1) = xh_(:,1) + KG*inov(:,1); %calculate posterior estimate
    
Px = Px_ - KG*Py*KG'; %calculate posterior covariance estimate
   
state = xh(:,1);
Pstate = Px;
for l = 1:Xdim
    if state(l,1) < 0 %check that states don't go negative
        state(l,1) = stateAll(l,1);
        Pstate = P_orig;
    end
end
xhat = state;
Phat = Pstate;
    %----------------------------------------------------------------------
   


end   



end

