function yk = AlbersNoiseFcn(xk,vk)
% Author: CLe
% LdeP Converting back to actual "non-additive" measurement noise form. 
% Date: 
% Summary:  Example measurement function for discrete time nonlinear
%           state estimators with non-additive measurement noise.
% Inputs:   xk = states at time k, x[k]
%           vk = measurement noise at time, k v[k]
% Outputs:  yk = y[k], measurements at time k
%
% The measurement is the third state (glucose) with multiplicative noise
yk = xk(3)*(1+vk);
%LdePTest
%yk = xk(3)-vk;
end
