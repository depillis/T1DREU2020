%From Table 1
J = 5*10^4; % cells ml^-1d^-1 : normal resting macrophage influx
k = 0.4; %d^-1: Ma deactivation rate
b = 0.09 ; %d^-1: recruitment rate of M by Ma
c = 0.1; %d^-1: macrophage egress rate
e1= 1*10^(-8); %cell^-1d^-1: anti-crowding terms for macrophages
e2= 1*10^(-8);
scale_factor = .0623; % works for scale factor between .0604 and .0635 (to that many sig figs)

%Scale factors for clearance rates
f1ns = 0.8948;%0.9155;
f1s = 2;
f2ns = 0.98;%1.0527;%Making this less than 1 moves away from the 250 boundary and smoothes the kink
f2s = 5;


%Clearance rates
%Macrophages
f1  = scale_factor*f1s*10^(-5); %ml cell^-1d^-1 %wt basal phagocytosis rate per M
f1n = scale_factor*f1ns*10^(-5); %ml cell^-1d^-1 %nod basal phagocytosis rate per M 
f2  = scale_factor*f2s*10^(-5);%ml cell^-1d^-1 %wt activated phagocytosis rate per Ma
f2n = scale_factor*f2ns*10^(-5); %ml cell^-1d^-1 % nod activated phagocytosis rate per Ma 
%DC conversion factors
DCtoM = 5.49 * 10^(-2);
tDCtoM = 3.82 * 10^(-1);
%DC clearance rates
fD = f2s*DCtoM*scale_factor*10^(-5);
ftD = f2s*tDCtoM*scale_factor*10^(-5);


% Added parameters
%alpha_B = 0.0334;
delta_B = 1/60;
% B_tot = 7.76*10^7; % Measured in cells, computed from Maree and Paredes
B_conv  = 2.59*10^5; % Units: cell/mg, Computed using the Maree value for density of cells in healthy pancreas and Topp value for mg beta cells in healthy pancreas
Ghb     = 90;
R0      = 864; %mg per dl
EG0     = 1.44; %per day
%SI      = .72; % ml per muU per day
sigmaI  = 43.2; %muU per ml per day per mg
deltaI  = 432; %per day
%GI      = sqrt(20000); %mg per dl
%Qpanc   = .194; % ml
Qblood  = 3.00; % ml
Qspleen = 0.10; % ml

mu_PB = 0.51;         % per day, original value
mu_BP = 0.1;          % per day
%D_ss = 10^5;


bDEday = 0.487e-5;  % ml/cell/day, per capita elimination rate of DC by CTL
bIRday = 0.487e-5;    %ml/cell/day, per capita elimination rate of iDC by Tregs
aEaday = .1199;     % per day; BUT I DON'T THINK THIS CAN BE RIGHT!  [ESTIMATED, =ln(2)/(5.78 days)]
T_naive = 370;   % cells, number of naive cells contributing to primary clonal expansion
bpday = 12;         % per day, maximal expansion factor of activated CTL
ramday = 0.01;      % per day, reversion rate of active CTL to memory CTL
baEday =  10^(-3); %BSH distinguish activation rate for effector memory cells by DC
baRday =  10^(-3); %BSH distinguish activation rate for Treg memore cell by iDC
aEmday = 0.0100;    % per day; death rate of memory CTL cells [ESTIMATED, =ln(2)/(69 days)]
%mues_r = aEmday/5000 ;     %per day;  rate of CTL-Treg removal equal to death rate of memory CTL
%mues_e = aEmday/5000;
thetaD = 2.12e5; % cells/ml;  Threshold in DC density in the spleen for half maximal proliferation rate of CTL [7.5e2,1.2e4]
d    = 0.5 ; % d^(-1)  rate of apoptotic cell decay into a necrotic cell

sE     = 1; %saturation terms for effector killing
sR     = 36; %control of trefs over effector killing
%eta    =eta_test;%0.01;%0.013;% effectiveness of effectors
%Parameters for eta_vary
eta = 0.01;
alpha_eta=0.11; %This value seems to control the jump slope
beta_eta=21;%22 % This afects the time to sick
wave_basal=0.75;