function visualise_bimodal_results(X, m, R, theta, no_common_scaling, layout, limits)
% Usage: visualise_bimodal_results(X, m, R, theta, [no_common_scaling,layout, limits])
%
% Visualise the outcome of adaptive_pt (first two coordinates). 
% Unless the fifth argument is given, there are common coordinate
% axis in all the figures.

% visualise_bimodal_results.m
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

if nargin<5
  no_common_scaling = false;
else
%  no_common_scaling = true;
end

if nargin<2
  vis_adapt = false;
  no_common_scaling = true;
else
  vis_adapt = true;
end

[d,L_x,N] = size(X);
if vis_adapt
  L = max(L_x, size(m,2));
else
  L = L_x;
end

x_all = permute(X(1,:,:), [2 3 1]);
y_all = permute(X(2,:,:), [2 3 1]);

if nargin<7
  limits = [min(x_all(:)), max(x_all(:)), min(y_all(:)), max(y_all(:))];
end

% Determine the layout for the figures
if nargin<6
  m_fig = ceil(sqrt(L_x));
  n_fig = ceil(L_x/m_fig);
else
  m_fig = layout(1);
  n_fig = layout(2);
end

% Ad-hoc extra scaling to make the shapes visible
if vis_adapt
max_ell_sz = [0;0];
m_theta = max(exp(theta));
for ell = 1:L
  R_ = R(1:2, 1:2, min(ell,size(R,3)));
  [Q, D] = eig(R_'*R_);
  % Get the projections to the coordinate axis, and take the dominating 
  % element
  sz_ell = sqrt(max(abs(Q*D), [], 2));
  max_ell_sz = max([max_ell_sz exp(theta(ell))*sz_ell], [], 2);
end
vis_scale = 0.5*min([diff(limits(1:2)); diff(limits(3:4))]./max_ell_sz);

% If the means are equal, use the global center point instead
if size(m,2)==1 | all(all(repmat(m(:,1),1,L)==m))
  m = repmat([mean(limits(1:2)); mean(limits(3:4))], 1, L);
end
end

for ell = 1:L_x
  subplot(n_fig, m_fig, ell);
  plot_single_level(x_all, y_all, limits, ell);
  if vis_adapt
  ind_R = min(ell,size(R,3));
  covariance_ellipsoid(m(1:2,ell), R(1:2,1:2,ind_R), vis_scale*exp(theta(ell)));
  end
  title(['level ' num2str(ell)]);
  if ~no_common_scaling
    axis(limits)
  end
end

% Plot the points of single level
function plot_single_level(x_all, y_all, limits, ell)
plot(x_all(ell,:),y_all(ell,:),'k.', 'markersize', 2);

% Draw covariance ellipsoid with given center, shape and scaling
function covariance_ellipsoid(center, R, theta)
n_ell = 100;
u = linspace(0,2*pi,n_ell);
X = repmat(center,1,n_ell) + theta*R'*[cos(u); sin(u)];
line(X(1,:),X(2,:),'linestyle','-','color','red');
