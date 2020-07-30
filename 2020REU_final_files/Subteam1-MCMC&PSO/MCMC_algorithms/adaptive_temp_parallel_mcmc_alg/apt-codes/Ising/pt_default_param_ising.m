function param = pt_default_param_ising(L)
% PT_DEFAULT_PARAM_ISING  Default parameters for adaptive parallel tempering
% for spatial imaging example
% Usage: param = pt_default_param_ising(L)
%
% In:
%  
%   L -- Number of temperature levels (default=12).

% pt_default_param_ising.m
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
  L = 12;
end
global IceFloe;
param = struct(...
'L', L, ...              % Number of levels
'alpha_sw', 0.234, ...    % Desired mean acceptance probability of a swap
'sw_exp_rm', 0.6, ...   % Exp. of temperature adaptation step size
'x', IceFloe, ...     % Starting point
'd_log_T', ones(1,L-1), ... % Initial log diff. of temperatures
'verbose', true, ...     % Whether to show statistics during the run
'swap_strategy', 'one_random_swap_ising', ... % The swap strategy: 
...                      % 'one_random_swap', 'random_order_sweep', or 'sweep_swap'
'save_memory', true, ... % Save memory by keeping only mean posterior probabilities
'save_temp_adapt', false, ... % Instead of final estimates, record 
...                      % the trajectories of temperature adaptation.
'nthin', 1, ...           % Thinning: Save only every N:th item.
'burn_in' , 10^4 ...        % burn in time , for save memory mode only
);
