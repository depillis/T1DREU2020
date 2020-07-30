function R = ram_adapt_shape(R, gamma, d, L, u, alpha_diff)
% Robust AM adaptation for each of the levels

% ram_adapt_shape.m
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

u_norm_sq = sum(u.^2,1);

for ell = 1:L
  R_ = R(:,:,ell);
  const = sqrt(min(0.9,d*gamma)*abs(alpha_diff(ell))/u_norm_sq(ell));
  if alpha_diff(ell) < 0
    R(:,:,ell) = cholupdate(R_, const*R_'*u(:,ell), '-');
  else
    R(:,:,ell) = cholupdate(R_, const*R_'*u(:,ell), '+');
  end
end
