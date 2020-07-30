%This is the main function which performs the UKF estimates


function [xh, Px, xh_, Px_, Y_, inov, Py, KG,negloglik,Pxall,Pyall,Px_all] = ukfaddmultiP(stateAll, Pstate,Q,R, obsAll, U1, U2, InferenceDS)

%Use the InferenceDS struct to set variables (meanings described in
%LorenzControlfin.m
ffun = InferenceDS.ffun;
hfun = InferenceDS.hfun;
alpha = InferenceDS.spkfParams(1);
beta = InferenceDS.spkfParams(2);
kappa = InferenceDS.spkfParams(3);
par = InferenceDS.par;
Odim = InferenceDS.obsdim;
Xdim = InferenceDS.Xdim;
Nsubj = InferenceDS.Nsubj;
pnoiseAdaptpar = InferenceDS.pNoiseAdaptpar;
partQflag = InferenceDS.partQflag;
Q = Q.cov;
R = R.cov;
NxNoPar = InferenceDS.NxNoPar;

Vdim  = size(Q,2);                                    % extract process noise dimension
Ndim  = size(R,2);                                   % extract observation noise dimension
U1dim = size(U1,1);                                  % Not used for state estimation
U2dim = size(U2,1);                                  % Not used for state estimation

NT = size(obsAll,2);                                           % number of input vectors
negloglik= 0;                                          %Initialize likelihood (used for parameter estimation)
Pxall = zeros(Xdim,Xdim,NT);                            %Initialize matrix which will hold ALL instances of posterior state covariances                       
Pyall = zeros(Odim,Odim,NT);                            %Initialize matrix which will hold ALL instances of observables covariances
Px_all = zeros(Xdim,Xdim,NT);                           %Initialize matrix which will hold ALL instances of prior state covariances

if (nargin ~= 8) error(' [ ukf ] Number of input arguments does not match.'); end

xhAll   = zeros(Xdim*Nsubj,NT);             %Will hold all xhat values for all subjects                    
inovAll = zeros(Odim*Nsubj,NT);             %Will hold all innovation values for all subjects

L = Xdim;                                                % augmented state dimension
nsp = 2*L+1;                                              % number of sigma-points
kappa = alpha^2*(L+kappa)-L;                              % compound scaling parameter

W = [kappa 0.5 0]/(L+kappa);                              % sigma-point weights
W(3) = W(1) + (1-alpha^2) + beta;

Sqrt_L_plus_kappa = sqrt(L+kappa);

UU1=zeros(0,nsp);
UU2=zeros(0,nsp);

Sv = [];
if (~isempty(Q) & partQflag==1) Sv=Q;end
if (~isempty(Q) & partQflag==0),
[Sv h] = chol(Q);                       %Perform choletsky factorization to take care of divergence issues
    if h==0 Sv= Sv'; 
    else Sv=zeros(Vdim,Vdim);end;
end

personlik=zeros(Nsubj,1);
[Sn h] = chol(R);                      %Perform choletsky factorization to take care of divergence issues
if h==0 Sn = Sn';else Sn=zeros(Ndim,Ndim);end
%--------------------------------------- Loop over all input vectors --------------------------------------------
for j=1:Nsubj,
obs  = obsAll((1+(j-1)*Odim):(Odim*j),:);
xh   = zeros(Xdim,NT);              %Initialize matrix of posterior estimates
xh_  = zeros(Xdim,NT);              %Initialize matrix of prior estimates
yh_  = zeros(Odim,NT);              %Initialize matrix of prior estimates of observables
inov = zeros(Odim,NT);
state=stateAll((1+(j-1)*Xdim):(Xdim*j));

    
for i=1:NT,                 %Loop over all time values
    
    if U1dim UU1 = cvecrep(U1(:,i),nsp); end
    if U2dim UU2 = cvecrep(U2(:,i),nsp); end

    [Sx h] = chol(Pstate); %Perform choletsky factorization
    if h==0 Sx=Sx'; else Sx = diag([repmat(.001,Xdim,1)]);end

    %------------------------------------------------------
    % TIME UPDATE       %
    Z    = cvecrep([state], nsp);            % copy needed for possible angle components section
    Sz  = [Sx];
    sSz  = Sqrt_L_plus_kappa * Sz;
    sSzM = [sSz -sSz];
    Z(:,2:nsp) = Z(:,2:nsp) + sSzM;           % build sigma-point set
    %-- Calculate predicted state mean and covariance
  
    X_ = feval(ffun, Z, UU1,[],i,par);  % propagate state sigma-points through process model
    X_bps = X_;

    xh_(:,i) = W(1)*X_(:,1) + W(2)*sum(X_(:,2:nsp),2);   %Calculate prior rediction using weighted mean
    temp1 = X_ - cvecrep(xh_(:,i),nsp);
    Px_ = W(3)*temp1(:,1)*temp1(:,1)' + W(2)*temp1(:,2:nsp)*temp1(:,2:nsp)';  %Calculate prior covariance using weights
    sSzV=zeros(Xdim,nsp);
    if ~isempty(Q)
    Px_ = Px_ + Q;
    S2 = [Sv];
    %---Calculate predicted observation and covariance
    sSv2 = Sqrt_L_plus_kappa*Sv;
    sSzV2 = [sSv2 -sSv2];
    sSzV(:,2:nsp) = sSzV(:,2:nsp) + sSzV2;
    end;
    
    X_bps(:,2:nsp) = X_bps(:,2:nsp) + sSzV(:,2:nsp);
    Y_ = feval(hfun, X_bps, UU2, [],i, par);    % propagate through observation model
    yh_(:,i) = W(1)*Y_(:,1) + W(2)*sum(Y_(:,2:nsp),2);  %calculate prior estimate of states using weights
    temp2 = Y_ - cvecrep(yh_(:,i),nsp);
    Py  = W(3)*temp2(:,1)*temp2(:,1)' + W(2)*temp2(:,2:nsp)*temp2(:,2:nsp)'; %calculate observables covariance using weights
    if ~isempty(R)
        Py = Py + R;end
    Pxy = W(3)*temp1(:,1)*temp2(:,1)' + W(2)*temp1(:,2:nsp)*temp2(:,2:nsp)';

    %------------------------------------------------------
    % MEASUREMENT UPDATE

    KG = Pxy / Py; %Calculate Kalman Gain (Chow equation 37)
    inov(:,i) = obs(:,i) - yh_(:,i); %Calculate innovation (difference between observable and prediction
    xh(:,i) = xh_(:,i) + KG*inov(:,i); %Calculate posterior estimate for states (Chow equation 38)
    Px = Px_ - KG*Py*KG'; %Calculate posterior covariance (Chow equation 39)
   
    state = xh(:,i);
    Pstate = Px;
    %----------------------------------------------------------------------
    

lik = .5*Odim*log(2*pi)+.5*log(det(Py))+.5*(inov(:,i)'*inv(Py)*inov(:,i));
personlik(j)=personlik(j)+lik;
negloglik=negloglik+lik; %Add term to likelihood for this time (used for parameter estimation)

if ~isempty(pnoiseAdaptpar) 
    Q = diag(max(pnoiseAdaptpar(1)*diag(Q),pnoiseAdaptpar(2))); 
    Sv = chol(Q);end %Use choletsky factorization

Pxall(:,:,i) = Px;   %Add current prior covariance matrix to set of all
Pyall(:,:,i) = Py;   %Add current observables covariance matrix to set of all
Px_all(:,:,i) = Px_; %Add current posterior covariance matrix to set of all
end   %--- for loop for t=1:T
xhAll  (((1+(j-1)*Xdim):(j*Xdim)),:)   = xh(:,1:NT); %Add subjects state estimate values to set of all
inovAll(((1+(j-1)*Odim):(j*Odim)),:) = inov; %Add subjects innovations to set of all


end %for loop for subject 1:Nsubj

