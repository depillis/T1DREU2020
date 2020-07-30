%% Create upper/lower bounds for PSO
%  Author:       C. Catlett
%  Date:         June 18, 2020
%  Desc:         Use eFAST analysis by A. Do to assign upper/lower bounds
%                to params; if no info on bounds available, use a default
%                of +/- a percentage (user input)

function [lb, ub, orig] = makebounds_eFAST(percent)
load('eFAST.mat');
table = eFAST;
D_ss         = table.D_ss;
EG0          = table.EG0;
GI           = table.GI;
Ghb          = table.Ghb;
J            = table.J; 
Qblood       = [3; 3*(1-percent); 3*(1+percent)];
Qpanc        = table.Qpanc; 
Qspleen      = [.1; .1*(1-percent); .1*(1+percent)];
R0           = table.R0; 
SI           = table.SI; 
T_naive      = table.T_naive; 
aEaday       = table.aEaday; 
aEmday       = table.aEmday; 
alpha_B      = table.alpha_B; 
alpha_eta    = table.alpha_eta; 
b            = table.b; 
bDEday       = table.b_DE; 
bIRday       = table.b_IR; 
baEday       = table.baEday; 
baRday       = table.baRday; 
beta_eta     = table.beta_eta; 
bpday        = table.bpday; 
c            = table.c; 
d            = table.d; 
deltaI       = table.deltaI; 
delta_B      = table.delta_B; 
e1           = table.e1; 
e2           = table.e2; 
f1ns         = [0.8948; 0.8948*(1-percent); 0.8948*(1+percent)];
f1s          = [2; 2*(1-percent); 2*(1+percent)];
f2ns         = [.98; .98*(1-percent); .98*(1+percent)];
f2s          = [5; 5*(1-percent); 5*(1+percent)];
fD           = table.fD; 
ftD          = table.ftD; 
k            = table.k; 
mu_BP        = table.mu_BP; 
mu_PB        = table.mu_PB; 
mues_e       = table.mu_E; 
mues_r       = table.mu_R; 
ramday       = table.ramday; 
sE           = table.sE; 
sR           = table.sR;
scale_factor = [.0623; .0604; .0635];
sigmaI       = table.sigmaI; 
thetaD       = table.thetaD; 
eta          = table.eta_basal;

lb = [D_ss(2), EG0(2), GI(2), Ghb(2), J(2), Qblood(2), Qpanc(2), Qspleen(2), ...
    R0(2), SI(2), T_naive(2), aEaday(2), aEmday(2), alpha_B(2), alpha_eta(2), ...
    b(2), bDEday(2), bIRday(2), baEday(2), baRday(2), beta_eta(2), bpday(2), ...
    c(2), d(2), deltaI(2), delta_B(2), e1(2), e2(2), f1ns(2), f1s(2), f2ns(2), ...
    f2s(2), fD(2), ftD(2), k(2), mu_BP(2), mu_PB(2), mues_e(2), mues_r(2), ...
    ramday(2), sE(2), sR(2), scale_factor(2), sigmaI(2), thetaD(2), eta(2)];

ub = [D_ss(3), EG0(3), GI(3), Ghb(3), J(3), Qblood(3), Qpanc(3), Qspleen(3), ...
    R0(3), SI(3), T_naive(3), aEaday(3), aEmday(3), alpha_B(3), alpha_eta(3), ...
    b(3), bDEday(3), bIRday(3), baEday(3), baRday(3), beta_eta(3), bpday(3), ...
    c(3), d(3), deltaI(3), delta_B(3), e1(3), e2(3), f1ns(3), f1s(3), f2ns(3), ...
    f2s(3), fD(3), ftD(3), k(3), mu_BP(3), mu_PB(3), mues_e(3), mues_r(3), ...
    ramday(3), sE(3), sR(3), scale_factor(3), sigmaI(3), thetaD(3), eta(3)];

orig = [D_ss(1), EG0(1), GI(1), Ghb(1), J(1), Qblood(1), Qpanc(1), Qspleen(1), ...
    R0(1), SI(1), T_naive(1), aEaday(1), aEmday(1), alpha_B(1), alpha_eta(1), ...
    b(1), bDEday(1), bIRday(1), baEday(1), baRday(1), beta_eta(1), bpday(1), ...
    c(1), d(1), deltaI(1), delta_B(1), e1(1), e2(1), f1ns(1), f1s(1), f2ns(1), ...
    f2s(1), fD(1), ftD(1), k(1), mu_BP(1), mu_PB(1), mues_e(1), mues_r(1), ...
    ramday(1), sE(1), sR(1), scale_factor(1), sigmaI(1), thetaD(1), eta(1)];
end
