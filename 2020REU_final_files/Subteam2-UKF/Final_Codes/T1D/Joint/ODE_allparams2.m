% The ODE model that describes the single compartment model

function dy = ODE_allparams(t,y,f1t,f2t,wave, par)
dy = zeros(53,1);

%Get parameter values
T1D_ODE_Parameters;
J = par(1);
k = par(2);
b = par(3);
c = par(4);
e1 = par(5);
e2 = par(6);
fD = par(7);
ftD = par(8);
alpha_B = par(9);
delta_B = par(10);
B_conv = par(11);
Ghb = par(12);
R0 = par(13);
EG0 = par(14);
SI = par(15);
sigmaI = par(16);
deltaI = par(17);
GI = par(18);
Qpanc = par(19);
Qblood = par(20);
Qspleen = par(21);
mu_PB = par(22);
mu_BP = par(23);
D_ss = par(24);
bDEday = par(25);
bIRday = par(26);
aEaday = par(27);
T_naive = par(28);
bpday = par(29);
ramday = par(30);
baEday = par(31);
baRday = par(32);
aEmday = par(33);
mu = par(34);
%eta = par(35);
mues_r = par(35); %unused placeholder
thetaD = par(36);
d = par(37);
sE = par(38);
sR = par(39);
alpha_eta = par(40);
beta_eta = par(41);

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
eta=0.018; %The value of eta in eta_vary nontrivially affects NOD outcomes with wave
%eta = 0.01;

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
            - ramday*E + baEday*D.*Em - mu.*E*R;
%R: Regulatory T cells
dy(11) = aEaday*(T_naive/Qspleen - R) + bpday*(tD.*R)./(thetaD + tD)...
            - ramday*R + baRday*tD.*Em - mu.*E.*R;
%Em: Memory T cells
dy(12) = ramday*(E+R) - (aEmday + baEday*D + baRday*tD).*Em;

dy(13:53) = 0;


end