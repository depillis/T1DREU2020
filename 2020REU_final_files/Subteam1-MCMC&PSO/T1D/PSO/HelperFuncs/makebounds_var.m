%% Create upper/lower bounds for PSO
%  Author:       C. Catlett
%  Date:         June 24, 2020
%  Desc:         Use variance values from D. Shenker, R. Wander to assign upper/lower bounds
%                to params; if no info on bounds available, use a default
%                of +/- a percentage (user input). Assumes range is 95%
%                confidence interval if param was distibuted normally.

function [lb, ub, orig] = makebounds_var(percent)
load('paramVariance.mat');
table = paramVariance;
D_ss         = makeInterval(1, table.D_ss, percent);
EG0          = makeInterval(2, table.EG0, percent);
GI           = makeInterval(3, table.GI, percent);
Ghb          = makeInterval(4, table.Ghb, percent);
J            = makeInterval(5, table.J, percent); 
Qblood       = makeInterval(6, table.Qblood, percent); 
Qpanc        = makeInterval(7, table.Qpanc, percent); 
Qspleen      = makeInterval(8, table.Qspleen, percent);
R0           = makeInterval(9, table.R0, percent); 
SI           = makeInterval(10, table.SI, percent); 
T_naive      = makeInterval(11, table.T_naive, percent); 
aEaday       = makeInterval(12, table.aEaday, percent); 
aEmday       = makeInterval(13, table.aEmday, percent); 
alpha_B      = makeInterval(14, table.alpha_B, percent); 
alpha_eta    = makeInterval(15, table.alpha_eta, percent); 
b            = makeInterval(16, table.b, percent); 
bDEday       = makeInterval(17, table.bDEday, percent); 
bIRday       = makeInterval(18, table.bIRday, percent); 
baEday       = makeInterval(19, table.baEday, percent); 
baRday       = makeInterval(20, table.baRday, percent); 
beta_eta     = makeInterval(21, table.beta_eta, percent); 
bpday        = makeInterval(22, table.bpday, percent); 
c            = makeInterval(23, table.c, percent); 
d            = makeInterval(24, table.d, percent); 
deltaI       = makeInterval(25, table.deltaI, percent); 
delta_B      = makeInterval(26, table.delta_B, percent); 
e1           = makeInterval(27, table.e1, percent); 
e2           = makeInterval(28, table.e2, percent); 
f1ns         = [0.8948; 0.8948*(1-percent); 0.8948*(1+percent)];
f1s          = [2; 2*(1-percent); 2*(1+percent)];
f2ns         = [.98; .98*(1-percent); .98*(1+percent)];
f2s          = [5; 5*(1-percent); 5*(1+percent)];
fD           = makeInterval(33, table.fD, percent); 
ftD          = makeInterval(34, table.ftD, percent); 
k            = makeInterval(35, table.k, percent); 
mu_BP        = makeInterval(36, table.mu_BP, percent); 
mu_PB        = makeInterval(37, table.mu_PB, percent); 
mues_e       = makeInterval(38, table.mues_e, percent); 
mues_r       = makeInterval(39, table.mues_r, percent); 
ramday       = makeInterval(40, table.ramday, percent); 
sE           = makeInterval(41, table.sE, percent); 
sR           = makeInterval(42, table.sR, percent);
scale_factor = [.0623; .0604; .0635];
sigmaI       = makeInterval(44, table.sigmaI, percent); 
thetaD       = makeInterval(45, table.thetaD, percent); 
eta          = [0.02; .01; .03];

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

% Returns vector of original value, lb, ub based on 95% confidence interval
% If negative, default to eFAST minimum
function bounds = makeInterval(ind, var, percent)
    [lb_temp,~,orig_temp] = makebounds_eFAST(percent);
    lb = orig_temp(ind)-2*sqrt(var);
    % Do not allow negative values
    if lb < 0
        lb = lb_temp(ind);
    end
    bounds = [orig_temp(ind); lb; orig_temp(ind)+2*sqrt(var)];
end