function [x, e, acc_sw, alpha_sw, ell_sw] = one_random_swap( ...
                         x, e, beta, L);

% one_random_swap.m
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

ell_sw = zeros(1,L-1);
alpha_sw = ell_sw;
acc_sw = alpha_sw;

% Pick random integer in [1,L-1]
ell = min(floor(rand*(L-1)+1),L-1);

[x, e, acc_sw(ell), alpha_sw(ell)] = try_swap(x, e, beta, ell);
ell_sw(ell) = 1;
