function [xh, Px, xh_, Px_, Y_, inov, Py, KG,negloglik,Pxall,Pyall,Px_all] = UKF(stateAll, Pstate,Q,R, obsAll, U1, U2, InferenceDS)


%Sets up parameters from InferenceDs Object (assigned in main file)
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
Vdim  = size(Q,2);                                    
Ndim  = size(R,2);                                   
U1dim = size(U1,1);
U2dim = size(U2,1);

NT = size(obsAll,2);                                           % number of input vectors
negloglik= 0;
Pxall = zeros(Xdim,Xdim,NT);
Pyall = zeros(Odim,Odim,NT);
Px_all = zeros(Xdim,Xdim,NT);

if (nargin ~= 8) error(' [ ukf ] Number of input arguments does not match.'); end %checking the number of arguments put in 

xhAll   = zeros(Xdim*Nsubj,NT);                          %setting up vector for all estimates
inovAll = zeros(Odim*Nsubj,NT);                          %setting up vector for all innovations (error)

L = Xdim;                                                % augmented state dimension
nsp = 2*L+1;                                              % number of sigma-points
kappa = alpha^2*(L+kappa)-L;                              % compound scaling parameter

W = [kappa 0.5 0]/(L+kappa);     % sigma-point weights - 
                                 %[Initial sigma point mean weight,
                                 % other sigma points mean and covariance
                                 % weight,
                                 %Initial sigma point covaraince weight]
                                                         
W(3) = W(1) + (1-alpha^2) + beta;

Sqrt_L_plus_kappa = sqrt(L+kappa);

UU1=zeros(0,nsp);
UU2=zeros(0,nsp);

%Cholesky Factorication of Q - Process noise covariance
Sv = [];
if (~isempty(Q) & partQflag==1) Sv=Q;end
if (~isempty(Q) & partQflag==0),
[Sv h] = chol(Q);
    if h==0 Sv= Sv'; 
    else Sv=zeros(Vdim,Vdim);end;
end

personlik=zeros(Nsubj,1);

%Cholesky factorization of R - Measurement noise covariance
[Sn h] = chol(R);
if h==0 Sn = Sn';else Sn=zeros(Ndim,Ndim);end
%--------------------------------------- Loop over all input vectors --------------------------------------------
for j=1:Nsubj,
obs  = obsAll((1+(j-1)*Odim):(Odim*j),:); %Set up vector for observations
xh   = zeros(Xdim,NT); %Set up vector for posterior estimations
xh_  = zeros(Xdim,NT); %Set up vector for prior estimations
yh_  = zeros(Odim,NT); %Set up vector for prior obsevation predictions
inov = zeros(Odim,NT); %Set up vector for innovations (error)
state=stateAll((1+(j-1)*Xdim):(Xdim*j),:); %Set up  initial state guess

%Loop over all time points    
for i=1:NT,
    
    
    if U1dim UU1 = cvecrep(U1(:,i),nsp); end
    if U2dim UU2 = cvecrep(U2(:,i),nsp); end
    
    %Cholesky factorization of covariance matrix
    [Sx h] = chol(Pstate);
    if h==0 Sx=Sx'; 
    else
        Sx = diag([repmat(.001,Xdim,1)]);
    end

    %------------------------------------------------------
    % Projection Step       %
    Z    = cvecrep([state], nsp);  %create vector of sigma points
    Sz  = [Sx];      %Set equal to covariance matrix
    sSz  = Sqrt_L_plus_kappa * Sz; %multiply by scaling compontent to create matrix for creating sigma points
    sSzM = [sSz -sSz]; %Create matrix of postitive and negative for sigma points abovea nd below
    Z(:,2:nsp) = Z(:,2:nsp) + sSzM;           % build sigma-point set by adding and subtracting the transformed matrix
    %-- Calculate predicted state mean and covariance
  
    tspan = [0:1];
    X_ = zeros(Xdim, nsp); %create matrix for transformed sigma points
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    %loop through all sigma points, pass them through the transformation
    %function (by solving thr ODE) and record them in a matrix
    for z = 1:nsp
        [~, solutions] = ode45(@(t, y) Lotka_Volterra_Model(t, y, Z(3:6,z)), tspan, Z(:,z), options); %solve ode with ode45
        X_temp = solutions(end,:);
        X_(:,z) = X_temp';
    end
    X_bps = X_;

    
    xh_(:,i) = W(1)*X_(:,1) + W(2)*sum(X_(:,2:nsp),2); %Set prior state estimate to be weighted sum of sigma points
    temp1 = X_ - cvecrep(xh_(:,i),nsp); %Difference between weighted mean and sigma points
    Px_ = W(3)*temp1(:,1)*temp1(:,1)' + W(2)*temp1(:,2:nsp)*temp1(:,2:nsp)'; %Set prior covariance estimate to weighted sum of differences 
                                                                             %between sigma points and weighted mean
    sSzV=zeros(Xdim,nsp); %Create
    %Add process noise to prior covariance matrix calculation
    if ~isempty(Q)
    Px_ = Px_ + Q;
    S2 = [Sv];
    %---Calculate predicted observation and covariance
    sSv2 = Sqrt_L_plus_kappa*Sv;
    sSzV2 = [sSv2 -sSv2];
    sSzV(:,2:nsp) = sSzV(:,2:nsp) + sSzV2; %Calculate sigma point set for predicted covariance
    end;
    
    X_bps(:,2:nsp) = X_bps(:,2:nsp) + sSzV(:,2:nsp);
    Y_ = feval(hfun, X_bps, UU2, [],i, par);    % propagate through observation model to get estimate of observations
    yh_(:,i) = W(1)*Y_(:,1) + W(2)*sum(Y_(:,2:nsp),2); %create weighted mean of observations from sigma points
    temp2 = Y_ - cvecrep(yh_(:,i),nsp); %Calculate differences between weighted mean and sigma points
    Py  = W(3)*temp2(:,1)*temp2(:,1)' + W(2)*temp2(:,2:nsp)*temp2(:,2:nsp)'; %Calculate covariance matrix of observations
    %Add measurement noise to covariance matrix of observations
    if ~isempty(R)
        Py = Py + R;end
    Pxy = W(3)*temp1(:,1)*temp2(:,1)' + W(2)*temp1(:,2:nsp)*temp2(:,2:nsp)'; %Calculate covariance matrix between states and observations

    %------------------------------------------------------
    % Update Step

    
    KG = Pxy * inv(Py); %Calculate Kalman Gain
    inov(:,i) = obs(:,i) - yh_(:,i); %Calculate innovation (error)
    xh(:,i) = xh_(:,i) + KG*inov(:,i); %Use the prior prediction state vector, innovation and Kalman Gain 
                                       %to calculate the posterior prediction
                                       %state vector
    
   %Make sure parameters are biologically consistent (greater than 0)
    if (xh(5,i) < 0) && i > 0
        xh(5,i) = xh(5,i-1);
    end
    if (xh(6,i) < 0) && i > 0
        xh(6,i) = xh(6,i-1);
    end
    
    
    
    
    Px = Px_ - KG*Py*KG'; %Use the prior prediction covariance matrix, observations covariance matrix, and Kalman Gain
                          %to calculate the posterior prediction for the
                          %covariance matrix
   
    %Update state and Pstate                      
    state = xh(:,i);
    Pstate = Px;
    %----------------------------------------------------------------------
    



if ~isempty(pnoiseAdaptpar) 
    Q = diag(max(pnoiseAdaptpar(1)*diag(Q),pnoiseAdaptpar(2))); 
    Sv = chol(Q);end

%Update matrices of all results with results from this time step
Pxall(:,:,i) = Px; 
Pyall(:,:,i) = Py;
Px_all(:,:,i) = Px_;
end   %--- for loop for t=1:T
xhAll  (((1+(j-1)*Xdim):(j*Xdim)),:)   = xh(:,1:NT);
inovAll(((1+(j-1)*Odim):(j*Odim)),:) = inov;


end %for loop for subject 1:Nsubj



