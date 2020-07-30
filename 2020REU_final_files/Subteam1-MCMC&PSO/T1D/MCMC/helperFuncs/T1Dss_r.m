%% Sum of Squares - Likelihood
%  Author:       Edited by M. Watanabe, from N. Tania (May 30, 2018)
%  Date:         June 2020
%  Desc:         Sum of squares of between model predictions and observed
%                data. Used in loglikelihood function.

function ss = T1Dss_r(params, data, modeltype, parVer, wave)
% EXTRACT GLUCOSE DATA
ydata=data(:,2);

% SOLVE ODE MODEL
[~, ymodel]= T1Dfun_r(modeltype,parVer,wave,params,data);

% SUM OF SQUARES
ss = sum((ymodel-ydata).^2);
end
