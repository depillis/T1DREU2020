function [m, R] = mean_cov_adapt(m, R, x, gamma, d, L, separate)
% The mean & covariance adaptation

% mean_cov_adapt.m
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

C_mean = zeros(d,d);
m_mean = zeros(d,1);
  
for ell = 1:L    
  if separate
    dx = x(:,ell) - m(:,ell);
    % Mean update
    m(:,ell) = m(:,ell) + gamma*dx;
    % Covariance
    R(:,:,ell) = cholupdate(sqrt(1-gamma)*R(:,:,ell), sqrt(gamma)*dx);
  else   
    dx = x(:,ell)-m;
    m_mean = m_mean + dx/L;
    C_mean = C_mean + dx*dx'/L;
  end
end

if ~separate
  % Mean
  m = (1-gamma)*m + gamma*m_mean;
  % Covariance
  R = chol((1-gamma)*R'*R + gamma*C_mean);
end

