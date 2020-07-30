% test script to try the adaptive parallel tempering

% test.m
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

% The name of the target distribution (file without .m)
target = 'log_multimodal_target';

% Parameters in a struct
param = pt_default_param(2, 5);

% Number of iterations
N = 10e3;

% Run the PT
[X, m, R, theta, beta, stats] = adaptive_pt(target, param, 10e3);

% Produce a figure similar to Figure 1 in the paper
visualise_bimodal_results(X, m, R, theta, 0, [3 2], [-5 15 -5 15]);

