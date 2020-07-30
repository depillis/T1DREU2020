%% T1D ODE System Definition (Mathews & Li Data)
%  Desc:         Single compartment ODE model for T1D

function dy = T1Dsys_noDist(t, y, f1t, f2t, wave, params)
dy = zeros(12, 1);
%Get parameter values
D_ss         = params(1);
EG0          = params(2);
GI           = params(3);
Ghb          = params(4); 
J            = params(5); 
Qpanc        = params(7); 
Qspleen      = params(8); 
R0           = params(9); 
SI           = params(10); 
T_naive      = params(11); 
aEaday       = params(12); 
aEmday       = params(13); 
alpha_B      = params(14); 
alpha_eta    = params(15); 
b            = params(16); 
bDEday       = params(17); 
bIRday       = params(18); 
baEday       = params(19); 
baRday       = params(20); 
beta_eta     = params(21); 
bpday        = params(22); 
c            = params(23); 
d            = params(24); 
deltaI       = params(25); 
delta_B      = params(26); 
e1           = params(27); 
e2           = params(28); 
fD           = params(33); 
ftD          = params(34); 
k            = params(35); 
mu_PB        = params(37); 
mues_e       = params(38); 
mues_r       = params(39); 
ramday       = params(40); 
sE           = params(41); 
sR           = params(42); 
sigmaI       = params(44); 
thetaD       = params(45); 
eta          = params(46);

% Conversion factor
B_conv = 2.59*10^5;

%Define state variables
M  = y(1); % cells/ml
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

eta_vary = eta+2*eta*(1 + tanh(alpha_eta*(t - beta_eta*7)));

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
dy(5) = (alpha_B*G.^2./(G.^2+Ghb^2)).*B - delta_B*B - W ...
         - eta_vary*((sE*E).^2./(1+(sE*E).^2+ (sR*R).^2)).*B; % modified saturation
          
%G: Glucose
dy(6) = R0 - (EG0 + (SI*I)).*G; 
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