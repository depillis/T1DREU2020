function yk = AlbersMeasFcn(xk,dt)
% Author: CLe
% Date: November 21, 2019
% Summary:  Example measurement function for discrete 
%           time nonlinear state estimators
% Inputs:   xk = states at time k, x[k]
%           dt = time step to take (LdeP - necessary for call to ODE45)
% Outputs:  yk = y[k], measurements at time k
%
% The measurement is the first state with multiplicative noise

% yk = [xk(1); xk(2); xk(3); xk(4); xk(5); xk(6)];
% Only send out glucose portion
% yk = xk(3);
% LdeP Only send out glucose portion and parameters
 yk = [xk(3); xk(7); xk(8); xk(9)];

end
