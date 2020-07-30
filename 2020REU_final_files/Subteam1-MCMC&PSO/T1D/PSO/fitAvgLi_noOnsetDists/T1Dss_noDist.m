%% Sum of squares for T1D (Mathews & Li Data)
%  Author:       C. Catlett
%  Date:         June 11, 2020
%  Desc:         Computes sum of squares for T1D system given options
%                (waveon, type, wICs); used as objective function

function ss = T1Dss_noDist(params, obsPts, waveon, type, wICs, time)

% Solve the ODE, find residual for glucose
[~, Y] = T1Dfun_noDist(params, waveon, type, wICs, time); 
res = (Y(:,6) - obsPts).*(Y(:,6) - obsPts);

% Square the error values and sum
ss = sum(res);
end
