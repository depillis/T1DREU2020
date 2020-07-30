function param = pt_default_param(d, L)
% PT_DEFAULT_PARAM  Default parameters for adaptive parallel tempering
%
% Usage: param = pt_default_param(d, L)
%
% In:
%   d -- Dimension (default=10).
%   L -- Number of temperature levels (default=12).

% pt_default_param.m
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


if nargin<1
  d = 10;
end
if nargin<2
  L = 12;
end

param = struct(...
'L', L, ...              % Number of temperature levels
'd', d, ...              % Dimension
'alpha_sw', 0.234, ...   % Desired mean acceptance probability of a swap
'alpha_rwm', 0.234, ...  % Desired mean acceptance probability of RWM step
'rwm_exp', 0.6, ...      % Exponent of random-walk adaptation step size
'rwm_fixed_p', 0, ...    % The probability of drawing from the fixed (initial) proposal
'sw_exp_rm', 0.6, ...    % Exp. of temperature adaptation step size
'theta', zeros(1, L), ... % Initial log-scalings of RWM
'R', eye(d), ...         % Initial Cholesky factor
'x', zeros(d,1), ...     % Starting point
'd_log_T', ones(1,L-1), ... % Initial log diff. of temperatures
'verbose', true, ...     % Whether to show statistics during the run
'rwm_proposal_gen', 'randn', ...   % Unscaled proposal variate generation in RWM
'swap_strategy', 'one_random_swap', ... % The swap strategy function
'loc_shape_adapt', 'mean_cov_adapt', ... % The mean-covariance update function
'separate_shape_adaptation', true,  ... % Whether to use a separate covariance 
...                      % for each temperature
'ram_adapt', false, ...  % Use the robust AM adaptation
'save_memory', false, ... % Save memory by keeping only the coolest chain
'save_temp_adapt', false, ... % Instead of final estimates, record 
...                      % the trajectories of temperature adaptation.
'nthin', 1 ...           % Thinning: Save only every N:th item.
);
