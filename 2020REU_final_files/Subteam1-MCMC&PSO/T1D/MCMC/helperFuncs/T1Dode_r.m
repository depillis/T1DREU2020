%% T1D ODE Model 
%  Author:       M. Watanabe edited from Shtylla et al. (2019)
%  Date:         June 2020
%  Desc:         The ODE model that describes the single compartment type 1
%                diabetes model.

function dy = T1Dode_r(t,y,f1t,f2t,wave,params,version)
dy = zeros(12,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
initParams;
eta=0.018; %The value of eta in eta_vary nontrivially affects NOD outcomes with wave

if version == 1 % D. Shenker sensitive, eta_vary (& IC) parameters
    e1          = params(1);
    e2          = params(2);
    delta_B     = params(3);
    SI          = params(4);
    GI          = params(5);
    mues_r      = params(6);
    mues_e      = params(7);

    % eta_vary parameters
    alpha_eta = params(8);
    beta_eta = params(9);
    eta = params(10);

elseif version == 2 % glucose, eta_vary, eFAST sensitive parameters
    % Gluose parameters
    R0          = params(1);
    EG0         = params(2);
    SI          = params(3);

    % eta_vary parameters (impacts apoptic wave)
    alpha_eta   = params(4);
    beta_eta    = params(5);

    % eFAST Sensitive parameters
    D_ss        = params(6);
    alpha_B     = params(7);
    Qpanc       = params(8);
    mues_r      = params(9); 
    mues_e      = params(10);

elseif version == 3 % all 41 parameters
    J = params(1);
    k = params(2);
    b = params(3);
    c = params(4);
    e1 = params(5);
    e2 = params(6);
    fD = params(7);
    ftD = params(8);
    alpha_B = params(9);
    delta_B = params(10);
    B_conv = params(11);
    Ghb = params(12);
    R0 = params(13);
    EG0 = params(14);
    SI = params(15);
    sigmaI = params(16);
    deltaI = params(17);
    GI = params(18);
    Qpanc = params(19);
    Qblood = params(20);
    Qspleen = params(21);
    mu_PB = params(22);
    mu_BP = params(23);
    D_ss = params(24);
    bDEday = params(25);
    bIRday = params(26);
    aEaday = params(27);
    T_naive = params(28);
    bpday = params(29);
    ramday = params(30);
    baEday = params(31);
    baRday = params(32);
    aEmday = params(33);
    mues_r = params(34);
    mues_e = params(35);
    thetaD = params(36);
    d = params(37);
    sE = params(38);
    sR = params(39);
    alpha_eta = params(40);
    beta_eta = params(41);
    eta = params(42);

else
   msg='Invalid parameter set: enter 1, 2 or 3';
     error(msg) 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define T1D model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define state variables
M = y(1); % cells/ml
Ma = y(2); % cells/ml
Ba = y(3); % cells/ml
Bn = y(4); % cells/ml

% Added state variables
B = y(5); % mg
G = y(6); % mg/dl
I = y(7); % mu U

D = y(8); % cells/ml
tD = y(9); % cells/ml

E = y(10); % cells/ml
R = y(11); % cells/ml
Em = y(12); % cells/ml

% Compute the apoptotic wave
W = wave*.1*B*exp(-((t-9)/9)^2); % Modified from Maree et al 2006
%W = wave*.1*exp(-((t-9)/9)^2); % Modified from Maree et al 2006
%alpha_eta=0.1; %This value seems to control the time to sick
%beta_eta=22; % This also might affect the time to sick
%eta=0.018; %The value of eta in eta_vary nontrivially affects NOD outcomes with wave

%eta=0.009+0.004775*randn;
eta_vary=eta+2*eta*(1 + tanh(alpha_eta*(t - beta_eta*7)));
%eta_vary=eta;
% W=0; % For testing with no wave


%M: Resting macrophages
dy(1) = J + (k+b)*Ma -c*M - f1t*M.*Ba - f1t*M.*Bn - e1*M.*(M+Ma);
%Ma: Activated macrophages
dy(2) = f1t*M.*Ba + f1t*M.*Bn - k*Ma - e2*Ma.*(M+Ma);
%Ba: Apoptotic beta cells
dy(3) = W*B_conv/Qpanc + eta_vary*((sE*E)^2./(1+(sE*E)^2 + (sR*R)^2)).*B*B_conv/Qpanc...
    - f1t*M.*Ba-f2t*Ma.*Ba - d*Ba ...
    + delta_B*B*B_conv/Qpanc - ftD*(D_ss-D).*Ba - fD*D.*Ba;
%Bn: Necrotic beta cells
dy(4) = d*Ba - f1t*M.*Bn - f2t*Ma.*Bn - ftD*(D_ss-D).*Bn - fD*D.*Bn;


%B: Healthy beta cells
%dy(5) = (alpha_B*G.^2./(G.^2+Ghb^2)).*B - delta_B*B - W ...
 %        - eta*((sE*E).^2./(1+(sE*E).^2+ (sR*R).^2)).*B; % modified saturation
     
dy(5) = (alpha_B*G.^2./(G.^2+Ghb^2)).*B - delta_B*B - W ...
         - eta_vary*((sE*E).^2./(1+(sE*E).^2+ (sR*R).^2)).*B; % modified saturation
     
     
%G: Glucose
dy(6) = R0 - (EG0 + (SI*I)).*G;%G
%I: Insulin
dy(7) = sigmaI*(G.^2/(G.^2 + GI.^2))*B - deltaI*I;


%D: Immunogenic DCs
dy(8) = ftD*Bn.*(D_ss-D) - bDEday*E.*D - mu_PB*D;%
%tD: Tolerogenic DCs
dy(9) = ftD*(D_ss-tD-D)*Ba - mu_PB*tD - ftD*Bn*tD - bIRday*R.*tD;


%E: Effector T cells
dy(10) = aEaday*(T_naive/Qspleen - E)+ bpday*(D.*E)./(thetaD + D)...
            - ramday*E + baEday*D.*Em - mues_e.*E*R;
%R: Regulatory T cells
dy(11) = aEaday*(T_naive/Qspleen - R) + bpday*(tD.*R)./(thetaD + tD)...
            - ramday*R + baRday*tD.*Em - mues_r.*E.*R;
%Em: Memory T cells
dy(12) = ramday*(E+R) - (aEmday + baEday*D + baRday*tD).*Em;


end