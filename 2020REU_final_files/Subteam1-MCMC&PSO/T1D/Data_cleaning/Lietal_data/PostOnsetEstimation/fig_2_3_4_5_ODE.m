% The ODE model that describes the single compartment model

function dy = fig_2_3_4_5_ODE(t,y,f1t,f2t,wave, par)
dy = zeros(12,1);

%{
D_ss = par(1);
alpha_B = par(2);
GI = par(3);
SI = par(4);
Qpanc = par(5);
mues_e = par(6);
mues_r = par(7);
%}


%Get parameter values
fig_2_3_4_5_Parameters;
 %{
R0 = par(1);
EG0 = par(2);
SI = par(3);
%}
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

%dy(13:19) = 0;


end