function [X, beta, stats] = adaptive_pt_ising(param, N)
% ADAPTIVE_PT  Adaptive parallel tempering algorithm for spatail imaging
% example
%  
% Usage: [X, beta, stats] = adaptive_pt_ising( param, N)
%
% In:
%   param  -- Parameters of the algorithm; see PT_DEFAULT_PARAM_ISING
%   N      -- Number of iterations
%
% Out:
%   X      -- Simulated variables in 4D array (40*40*levels*N)
%   beta   -- The final inverse temperatures (vector of length L)
%   stats  -- Acceptance rate statistics etc.

% adaptive_pt_ising.m
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
nthin = param.nthin;
alpha_sw = param.alpha_sw;
sw_exp_rm = param.sw_exp_rm;
verbose = param.verbose;
NN=param.burn_in;
save_memory = param.save_memory;
save_temp_adapt = param.save_temp_adapt;

% The initial points for all levels
x=zeros(40,40,L);
for ell= 1:L
    x(:,:,ell)=param.x;
end



% Initial temperature differences
d_log_T = param.d_log_T;
beta = inv_temperatures(d_log_T, L);
if save_temp_adapt
  Beta = zeros(L, N+1);
  Beta(:,1) = beta';
end

% Initialise some auxiliary variables
alpha = zeros(1,L);

% Evaluate the initial values of the log-densities
e = zeros(1,L);
for ell=1:L
  e(ell) = log_ising_target(x(:,:,ell));
  if ~isfinite(e(ell))
    warning('Initial value of the log-density infinite')
  end
end

% Array containing the (thinned) simulated samples
if save_memory
  X = zeros(40,40,L);
else
  X = zeros(40,40,L,floor(N/nthin));
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
%  gamma = (k+1)^(-rwm_exp);
  
  % Metropolis whithin Gibbs moves:
  u = randi(40,2,L);
  for ell=1:L
    x_(:,:,ell) = x(:,:,ell);
    if x(u(1,ell),u(2,ell),ell)==1
    x_(u(1,ell),u(2,ell),ell)=0;
    else
    x_(u(1,ell),u(2,ell),ell)=1;    
  
    end
  end
  
  ru = rand(L,1);
  for ell = 1:L
    e_ = e(ell)+update_ising_target(x(:,:,ell),u(1,ell),u(2,ell));
    alpha(ell) = min(1,exp(beta(ell)*(e_-e(ell))));
    if ru(ell) <= alpha(ell)
      % accept
      e(ell) = e_;
      x(:,:,ell) = x_(:,:,ell);
      acc(ell) = acc(ell) + 1;
    end
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
     if k>NN
      X = X+1/(k-NN)*(x-X);
     end
    end
    else
      X(:,:,:,k/param.nthin) = x;
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

