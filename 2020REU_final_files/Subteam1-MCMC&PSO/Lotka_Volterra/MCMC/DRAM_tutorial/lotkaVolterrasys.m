%% Lotka-Volterra Equations
% Christina Catlett (June 4, 2020)
% 
% Set up of Lotka-Volterra system; used anonymously in ODE solver
% (lotkaVolterra.fun)

function ydot = lotkaVolterrasys(t,y,params)
y1  = y(1);     % Hares
y2  = y(2);     % Lynx

% Load parameters
alpha = params(1);
beta  = params(2);
gamma = params(3);
delta = params(4);

% Differential equations
dy1 = alpha * y1 - beta * y1 * y2;
dy2 = -gamma * y2 + delta * y2 * y1;

ydot = [dy1; dy2];
end