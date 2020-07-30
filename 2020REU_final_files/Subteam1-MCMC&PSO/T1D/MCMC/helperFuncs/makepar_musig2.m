%% Extract post-UKF parameter means/vars
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         Helper function to extract post-UKF means and variances
function [mu, sigma2] = makepar_musig2
raw = readmatrix('Parameter_Distributions_updated.csv');
mu = raw(:,1);
sigma2 = raw(:,2).^2;
end