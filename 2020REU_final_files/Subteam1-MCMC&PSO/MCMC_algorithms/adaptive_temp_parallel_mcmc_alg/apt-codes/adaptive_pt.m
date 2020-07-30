function [X, m, R, theta, beta, stats] = adaptive_pt(target, param, N)
% ADAPTIVE_PT  Adaptive parallel tempering algorithm
%  
% Usage: [X, m, R, theta, beta, stats] = adaptive_pt(target, param, N)
%
% In:
%   target -- Function returning log-density value at a given point.
%             (function handle or string)
%   param  -- Parameters of the algorithm; see PT_DEFAULT_PARAM
%   N      -- Number of iterations
%
% Out:
%   X      -- Simulated variables in 3D array (dim*levels*N)
%   m      -- The final proposal means (vectors of length dim)
%   R      -- The final proposal covariances (dim*dim matrices)
%   theta  -- The final scaling factors (vector of length L)
%   beta   -- The final inverse temperatures (vector of length L)
%   stats  -- Acceptance rate statistics etc.

% adaptive_pt.m
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

% Shorthand notations
L = param.L;
d = param.d;
nthin = param.nthin;
alpha_sw = param.alpha_sw;
alpha_rwm = param.alpha_rwm;
rwm_exp = param.rwm_exp;
sw_exp_rm = param.sw_exp_rm;
verbose = param.verbose;
separate_shape_adaptation = param.separate_shape_adaptation | param.ram_adapt;
ram_adapt = param.ram_adapt;
save_memory = param.save_memory;
save_temp_adapt = param.save_temp_adapt;

% The initial points for all levels
if size(param.x, 2) == 1
  x = repmat(param.x, 1,  L);
else
  x = param.x;
end

% Initial values for the mean, covariance and scalings
if separate_shape_adaptation
  m = x;
  if size(param.R,3) == 1
    R0 = repmat(param.R, [1, 1, L]);
  else
    R0 = param.R;
  end
else
  m = mean(x,2);
  R0 = param.R;
end
R = R0;

theta0 = param.theta;
theta = theta0;

% Initial temperature differences
d_log_T = param.d_log_T;
beta = inv_temperatures(d_log_T, L);
if save_temp_adapt
  Beta = zeros(L, N+1);
  Beta(:,1) = beta';
end

% Initialise some auxiliary variables
C_mean = zeros(d);
m_mean = zeros(d,1);
alpha = zeros(1,L);

% Evaluate the initial values of the log-densities
e = zeros(1,L);
for ell=1:L
  e(ell) = feval(target, x(:,ell));
  if ~isfinite(e(ell))
    warning('Initial value of the log-density infinite')
  end
end

% Array containing the (thinned) simulated samples
if save_memory
  X = zeros(d,1,floor(N/nthin));
else
  X = zeros(d,L,floor(N/nthin));
end

% Acceptance statistics for RWM and switch:
acc = zeros(1,L);
acc_sw = zeros(1,L-1);
N_sw_stat = zeros(1,L-1);

% Auxiliary variables for the mean acceptance probability 
% of switch between levels
alpha_sw_est = zeros(1,L-1);
N_sw = zeros(1,L-1);

tic
if verbose, t=toc; end
for k = 1:N
  gamma = (k+1)^(-rwm_exp);
  
  % Random-walk Metropolis moves:
  u = feval(param.rwm_proposal_gen, d, L);
  fixed_component = (rand(1,L) <= param.rwm_fixed_p);
  ell_R = 1;
  for ell=1:L
    if separate_shape_adaptation, ell_R = ell; end
    if fixed_component(ell)
      S = exp(theta0(ell))*R0(:,:,ell_R)';
    else
      S = exp(theta(ell))*R(:,:,ell_R)';
    end
    x_(:,ell) = x(:,ell) + S*u(:,ell);
  end
  
  ru = rand(L,1);
  for ell = 1:L
    e_ = feval(target, x_(:,ell));
    alpha(ell) = min(1,exp(beta(ell)*(e_-e(ell))));
    if ru(ell) <= alpha(ell)
      % accept
      e(ell) = e_;
      x(:,ell) = x_(:,ell);
      acc(ell) = acc(ell) + 1;
    end
  end
  
  alpha_diff = (~fixed_component).*(alpha-alpha_rwm);
  if ram_adapt
    R = ram_adapt_shape(R, gamma, d, L, u, alpha_diff);
  else  
    % Adapt the scales of the level
    theta = theta + gamma*alpha_diff;
    % Do the location-shape adaptation:
    [m, R] = feval(param.loc_shape_adapt, m, R, x, gamma, d, L, ...
    separate_shape_adaptation);
  end
  
  % Swap states between temperatures feval(param.swap_temperature
  [x, e, acc_sw_, alpha_sw_, ell_sw] = feval(param.swap_strategy, ...
  x, e, beta, L);
  alpha_sw_est = alpha_sw_est + alpha_sw_;
  N_sw = N_sw + ell_sw;
  acc_sw = acc_sw + acc_sw_;
  
  % Adaptation of the temperature schedule
  gamma = (k+1)^(-sw_exp_rm);    
  d_log_T = d_log_T + ...
            ((L-1)/sum(ell_sw))*gamma*ell_sw.*(alpha_sw_ - alpha_sw);
  beta = inv_temperatures(d_log_T, L);
  
  if save_temp_adapt
    Beta(:, k+1) = beta';
  end
  
  % Record the statistics
  N_sw_stat = N_sw_stat + N_sw;
  N_sw(:) = 0;
  alpha_sw_est(:) = 0;
  
  if rem(k, param.nthin) == 0
    if save_memory
      X(:,1,k/param.nthin) = x(:,1);
    else
      X(:,:,k/param.nthin) = x;
    end
  end

  if verbose & toc-t>1 
    t = toc;
    fprintf('\rElapsed %.0fs, remaining %.0fs; accepted %.2f%% (%.2f%%)  ', ...
            t, t/k*(N-k), mean(acc./k)*100, mean(acc_sw./N_sw_stat)*100 );
  end

end
N_sw_stat = N_sw_stat + N_sw;

stats = struct('acc_rwm', acc/N, 'acc_swaps', acc_sw./N_sw_stat, ...
'time_elapsed', toc);
if verbose
  fprintf('\n');
  fprintf('Accepted: '); fprintf('%3.2f%% ',stats.acc_rwm*100); fprintf('\n');
  fprintf('Swaps:    '); fprintf('%3.2f%% ',stats.acc_swaps*100); fprintf('\n');
end

if save_temp_adapt
  beta = Beta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the exponents of the distributions (inverse temperatures)
function beta = inv_temperatures(d_log_T, L)
T = ones(1,L);
for k = 1:L-1
  T(k+1) = T(k) + exp(d_log_T(k));
end
beta = 1./T;

