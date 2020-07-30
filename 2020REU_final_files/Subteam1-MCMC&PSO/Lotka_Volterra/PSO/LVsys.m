%% Lotka Volterra ODE System
%  Author:      Christina Catlett
%  Date:        June 4, 2020
%  Desc:        Set up of Lotka-Volterra system; used anonymously in ODE solver (LVfun.m) 

function ydot = LVsys(t,y,params)
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