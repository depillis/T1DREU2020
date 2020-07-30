%% ODE Solver
% N. Tania (May 30, 2018), edited by Christina Catlett (June 3, 2020)
%
% Solve the ODE using the given parameters (theta) and initial conditions
% (from ydata) for times (time)

function y = lotkaVolterrafun(time, theta, ydata)
% extract initial conditions
% ydata does not include time column
y0 = ydata(1,:);
[~,y] = ode45(@(t,y) lotkaVolterrasys(t, y, theta), time, y0);
