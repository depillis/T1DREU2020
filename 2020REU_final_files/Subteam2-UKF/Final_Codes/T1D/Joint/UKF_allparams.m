function [xh, Px, xh_, Px_, Y_, inov, Py, KG,negloglik,Pxall,Pyall,Px_all] = ukfaddmultiP(stateAll, Pstate,Q,R, obsAll, U1, U2, InferenceDS, tspans)


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
NxNoPar = InferenceDS.NxNoPar;


%params for ODE
T1D_ODE_Parameters;

%Set wave 
wave = wave_basal;

Vdim  = size(Q,2);                                    % extract process noise dimension
Ndim  = size(R,2);                                   % extract observation noise dimension
U1dim = size(U1,1);
U2dim = size(U2,1);

[NT y] = size(obsAll);% number of input vectors
NT
negloglik= 0;
Pxall = zeros(Xdim,Xdim,NT);
Pyall = zeros(Odim,Odim,NT);
Px_all = zeros(Xdim,Xdim,NT);

if (nargin ~= 9) error(' [ ukf ] Number of input arguments does not match.'); end %checking the number of arguments put in 

xhAll   = zeros(Xdim*Nsubj,NT); %setting up vector for all estimates
inovAll = zeros(Odim*Nsubj,NT);  %setting up vector for all innovations (error)

L = Xdim;   % augmented state dimension                                             % augmented state dimension
nsp = 2*L+1;  % number of sigma-points                                            % number of sigma-points
kappa = alpha^2*(L+kappa)-L; % compound scaling parameter                             % compound scaling parameter

W = [kappa 0.5 0]/(L+kappa);    % sigma-point weights - 
                                %[Initial sigma point mean weights other sigma points mean and covariance weight, Initial sigma point covaraince weight]               
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
obs  = obsAll'; %Set up vector for observations
xh   = zeros(Xdim,NT); %Set up vector for posterior estimations
xh_  = zeros(Xdim,NT); %Set up vector for prior estimations
yh_  = zeros(Odim,NT);  %Set up vector for prior obsevation predictions
inov = zeros(Odim,NT); %Set up vector for innovations (error)
%(1+(j-1)*Xdim)
%(Xdim*j)
%size(stateAll)
%state=stateAll((1+(j-1)*Xdim):(Xdim*j),:);
state = stateAll; %Set up  initial state guess
state = state';
prevstate = state;

%Loop over all time points     
for i=1:NT,
    i
    
    if U1dim UU1 = cvecrep(U1(:,i),nsp); end
    if U2dim UU2 = cvecrep(U2(:,i),nsp); end

    %Cholesky factorization of covariance matrix
    [Sx h] = chol(Pstate);
    if h==0 Sx=Sx'; 
    %else Sx = diag([repmat(.001,Xdim,1)]); 
    else
    fprintf('not PD');
    %Sx = diag([repmat(.01,Xdim,1)]);
    %Sx = [ diag([repmat(1,12,1)]) zeros(12,7); zeros(7, 12) diag([10^5 10^-2 5 10^-1 10^-1 10^-(5.5) 10^-(5.5)])];
    Sx = [diag([repmat(.001,12,1)]) zeros(12,41); zeros(41, 12) diag([10^1 0.0001 0.00001 0.0001 10e-11 10e-10 10e-12 10e-12 10e-9 10e-6 1 .1 .1 .0001 .1 .1 .2 50 10^-6 ...
                      10^-6 10^-6 10^-6 10^-6 100 10e-11 10e-11 10e-6 .1 .01 10e-6 10^-8 10^-8 10^-9 10^-6 10^-5 ...
                      50 0.0001 0.0001 .05 10 8])];
    end

    %------------------------------------------------------
     % Projection Step       %
    Z    = cvecrep([state], nsp); %create vector of sigma points
    Sz  = [Sx];  %Set equal to covariance matrix
    sSz  = Sqrt_L_plus_kappa * Sz; %multiply by scaling compontent to create matrix for creating sigma points
    sSzM = [sSz -sSz]; % build sigma-point set by adding and subtracting the transformed matrix
    %size(Z(:,2:nsp))
    %size(Z(:,2:nsp))
    %size(sSzM)
    Z(:,2:nsp) = Z(:,2:nsp) + sSzM;           % build sigma-point set
    %-- Calculate predicted state mean and covariance
  
    tspan = [0:tspans(i)]; %create tspan for solving system
    X_ = zeros(Xdim, nsp); %create matrix for transformed sigma points
    options=odeset('RelTol',1e-5,'AbsTol',1e-5);
    %loop through all sigma points, pass them through the transformation
    %function (by solving thr ODE) and record them in a matrix
    for z = 1:nsp
        [~, solutions] = ode15s(@(t, y) ODE_allparams2(t, y, f1n, f2n, wave, Z(13:53,z)), tspan, Z(:,z), options); %solve ode with ode15s
        X_temp = solutions(end,:);
        X_(:,z) = X_temp';
    end
    X_bps = X_;
    
    
    xh_(:,i) = W(1)*X_(:,1) + W(2)*sum(X_(:,2:nsp),2); %Set prior state estimate to be weighted sum of sigma points
    temp1 = X_ - cvecrep(xh_(:,i),nsp); %Difference between weighted mean and sigma points
    Px_ = W(3)*temp1(:,1)*temp1(:,1)' + W(2)*temp1(:,2:nsp)*temp1(:,2:nsp)'; %Set prior covariance estimate to weighted sum of differences 
                                                                             %between sigma points and weighted mean
    sSzV=zeros(Xdim,nsp);
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
    yh_(:,i) = W(1)*Y_(:,1) + W(2)*sum(Y_(:,2:nsp),2);  %create weighted mean of observations from sigma points
    temp2 = Y_ - cvecrep(yh_(:,i),nsp); %Calculate differences between weighted mean and sigma points
    Py  = W(3)*temp2(:,1)*temp2(:,1)' + W(2)*temp2(:,2:nsp)*temp2(:,2:nsp)'; %Calculate covariance matrix of observations
     %Add measurement noise to covariance matrix of observations
    if ~isempty(R)
        Py = Py + R;
    end
    Pxy = W(3)*temp1(:,1)*temp2(:,1)' + W(2)*temp1(:,2:nsp)*temp2(:,2:nsp)'; %Calculate covariance matrix between states and observations

    %------------------------------------------------------
    % Update Step

    KG = Pxy * inv(Py); %Calculate Kalman Gain
    inov(:,i) = obs(:,i) - yh_(:,i); %Calculate innovation (error)
    xh(:,i) = xh_(:,i) + KG*inov(:,i); %Use the prior prediction state vector, innovation and Kalman Gain 
                                       %to calculate the posterior prediction
                                       %state vector
    %Make sure atates and parameters are biologically consistent (greater than 0)
    if xh(1,i) < 0 
        xh(1,i) = prevstate(1);
    end
    
    if (xh(2,i) < 0 || xh(2,i) > 6e6) 
        xh(2,i) = prevstate(2);
    end
    
    if (xh(3,i) < 0 || xh(3,i) > 1e8) 
        xh(3,i) = prevstate(3);
    end
    
    if (xh(4,i) < 0 || xh(4,i) > 1e7) 
        xh(4,i) = prevstate(4);
    end
    
     if xh(6,i) < 0 || xh(6,i) > 1000
        xh(6,i) = prevstate(6);
    end
    
    if (xh(8,i) < 0 || xh(8,i) > 1e6) 
        xh(8,i) = prevstate(8);
    end
    
    if (xh(9,i) < 0 || xh(9,i) > 1e6) 
        xh(9,i) = prevstate(9);
    end
    
    if (xh(10,i) > 1e8 || xh(10,i) < 0) 
        xh(10,i) = prevstate(10);
    end
    
    if (xh(11,i) > 1e5 || xh(11,i) < 0)
        xh(11,i) = prevstate(11);
    end
    
    if (xh(12,i) < 0 || xh(12,i) > 1e5) 
        xh(12,i) = prevstate(12);
    end
    
     if xh(13,i) < 0
        xh(13,i) = prevstate(13);
     end
    
    if xh(14,i) < 0
        xh(14,i) = prevstate(14);
    end
    
     if xh(15,i) < 0
        xh(15,i) = prevstate(15);
     end
    
     if xh(16,i) < 0
        xh(16,i) = prevstate(16);
     end
    
     if xh(17,i) < 0
        xh(17,i) = prevstate(17);
     end
    
    if xh(18,i) < 0
        xh(18,i) = prevstate(18);
    end
    
    if xh(19,i) < 0
        xh(19,i) = prevstate(19);
    end
    
    if xh(20,i) < 0
        xh(20,i) = prevstate(20);
    end
    
    if xh(21,i) < 0
        xh(21,i) = prevstate(21);
    end
    
    if xh(22,i) < 0
        xh(22,i) = prevstate(22);
    end
    
    if xh(23,i) < 0
        xh(23,i) = prevstate(23);
    end
    
    if xh(24,i) < 0
        xh(24,i) = prevstate(24);
    end
    
    if xh(25,i) < 0
        xh(20,i) = prevstate(20);
    end
    
    if xh(25,i) < 0
        xh(25,i) = prevstate(25);
    end
    
    if xh(26,i) < 0
        xh(26,i) = prevstate(26);
    end
    
    if xh(27,i) < 0
        xh(27,i) = prevstate(27);
    end
    
    if xh(28,i) < 0
        xh(28,i) = prevstate(28);
    end
    
    if xh(29,i) < 0
        xh(29,i) = prevstate(29);
    end
    
    if xh(30,i) < 0
        xh(30,i) = prevstate(30);
    end
    
    if xh(31,i) < 0
        xh(31,i) = prevstate(31);
    end
    
    if xh(32,i) < 0
        xh(32,i) = prevstate(32);
    end
    
    if xh(33,i) < 0
        xh(33,i) = prevstate(33);
    end
    
    if xh(34,i) < 0
        xh(34,i) = prevstate(34);
    end
    
    if xh(35,i) < 0
        xh(35,i) = prevstate(35);
    end
    
    if xh(36,i) < 0
        xh(36,i) = prevstate(36);
    end
    
    if xh(37,i) < 0
        xh(37,i) = prevstate(37);
    end
    
    if xh(38,i) < 0
        xh(38,i) = prevstate(38);
    end
    
    if xh(39,i) < 0
        xh(39,i) = prevstate(39);
    end
    
    if xh(40,i) < 0
        xh(40,i) = prevstate(40);
    end
    
    if xh(41,i) < 0
        xh(41,i) = prevstate(41);
    end
    
    if xh(42,i) < 0
        xh(42,i) = prevstate(42);
    end
    
    
    if xh(43,i) < 0
        xh(43,i) = prevstate(43);
    end
    
    if xh(44,i) < 0
        xh(44,i) = prevstate(44);
    end
    
    if xh(45,i) < 0
        xh(45,i) = prevstate(45);
    end
    
    if xh(46,i) < 0
        xh(46,i) = prevstate(46);
    end
    
    if xh(47,i) < 0
        xh(47,i) = prevstate(47);
    end
    
    if xh(48,i) < 0
        xh(48,i) = prevstate(48);
    end
    
    if xh(49,i) < 0
        xh(49,i) = prevstate(49);
    end
    
    if xh(50,i) < 0
        xh(50,i) = prevstate(50);
    end
    
    if xh(51,i) < 0
        xh(51,i) = prevstate(51);
    end
    
    if xh(52,i) < 0
        xh(52,i) = prevstate(52);
    end
    
    if xh(53,i) < 0
        xh(53,i) = prevstate(53);
    end
    
    
    
    
    
    Px = Px_ - KG*Py*KG'; %Use the prior prediction covariance matrix, observations covariance matrix, and Kalman Gain
                          %to calculate the posterior prediction for the
                          %covariance matrix
   
    %Update state and Pstate   
    state = xh(:,i);
    Pstate = Px;
    prevstate = state;
    %----------------------------------------------------------------------
    

lik = .5*Odim*log(2*pi)+.5*log(det(Py))+.5*(inov(:,i)'*inv(Py)*inov(:,i));
personlik(j)=personlik(j)+lik;
negloglik=negloglik+lik;

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

