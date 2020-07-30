function p = log_multimodal_target(x)

% log_multimodal_target.m
% 
% Copyright (c) Blazej Miasojedow, Eric Moulines and Matti Vihola 2012
% 
% This file is part of the implementation of the Adaptive Parallel 
% tempering algorithm (hereafter "APT"); see http://arxiv.org/abs/1205.1076.
% 
% APT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% APT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with APT.  If not, see <http://www.gnu.org/licenses/>.

% The means of the bivariate 20-mode target
m=[2.18, 5.76 ; 3.25 ,3.47; 5.41, 2.65; 4.93, 1.50; ...
 8.67, 9.59 ; 1.70, 0.50;  2.70, 7.88;  1.83, 0.09; ...
 4.24, 8.48;  4.59, 5.60;  4.98, 3.70;  2.26, 0.31; ...
8.41, 1.68;  6.91, 5.81;  1.14, 2.39;   5.54, 6.86; ...
 3.93, 8.82;  6.87, 5.40;  8.33, 9.50;  1.69, 8.11]';

% The variance of the mixture components
sigma2 = 0.1^2;
% The variance of the dims d>=3
sigma_xtra_dims2 = sigma2;

% Evaluate the log-density of each of the Gaussians
l_dens = -0.5/sigma2*sum((m-repmat(x(1:2),1,20)).^2, 1) ...
         -0.5/sigma_xtra_dims2*sum(x(3:end).^2);

% Prevent underflow by normalised log-sum
l_max = max(l_dens);
p = l_max + log(sum(exp(l_dens-l_max)));
