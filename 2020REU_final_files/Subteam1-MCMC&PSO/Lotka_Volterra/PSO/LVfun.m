%% ODE Solver for Lotka-Volterra
%  Author:      N. Tania, edited by Christina Catlett
%  Date:        May 30th, 2018; edited June 3, 2020
%  Desc:        Solve the LV ODEs using the given parameters (params) and initial conditions
%               (from ydata) for series of times (time); Used in fitLV.m

function y = LVfun(time, theta, ydata)
% Extract initial conditions
% Note: ydata does not include time column
y0 = ydata(1,:);
[~,y] = ode45(@(t,y) LVsys(t, y, theta), time, y0);

