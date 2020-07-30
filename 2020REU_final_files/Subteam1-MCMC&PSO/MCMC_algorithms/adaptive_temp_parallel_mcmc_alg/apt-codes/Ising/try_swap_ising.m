function [x, e, acc_sw, alpha] = try_swap_ising( ...
                         x, e, beta, ell);
% Try to swap the levels ell <-> ell+1

% try_swap_ising.m
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

% The acceptance probability:
alpha = min(1, exp((beta(ell)-beta(ell+1))*(e(ell+1)-e(ell))));

% An accept-reject for the swap
if rand <= alpha
  % Swap states
  x_ = x(:,:,ell);
  x(:,:,ell) = x(:,:,ell+1);
  x(:,:,ell+1) = x_;
  
  % Swap also the energies accordingly
  e_ = e(ell);
  e(ell) = e(ell+1);
  e(ell+1) = e_;
  
  % Update statistics
  acc_sw = 1;
else
  acc_sw = 0;
end
