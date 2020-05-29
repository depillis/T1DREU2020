function yk = AlbersNoiseFcn(xk)
% Author: CLe
% LdeP - trying out making this "additive" noise ... so we don't specify vk.In
% this case, I believe it is assumed that vk is just some Gaussian (with
% mean 0), and that just gets added.
% LdeP Update - the noise here is actually additive in this case, but the
% function calls are treating it as non-additive.
% Date: 
% Summary:  Example measurement function for discrete time nonlinear
%           state estimators with non-additive measurement noise.
% Inputs:   xk = states at time k, x[k]
%           vk = measurement noise at time, k v[k]
% Outputs:  yk = y[k], measurements at time k
%
% The measurement is the third state (glucose) with "additive" noise - here
% we don't add anything special, the vk is automatically generated as a Gaussian (I think) by the
% UKF algorithm.
 yk = xk(3);

%LdePTest
%yk = xk(3)-vk;
%yk = xk(3)*(1+vk);
end
