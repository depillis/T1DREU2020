%% Sum of squares for T1D
%  Author:       C. Catlett
%  Date:         June 11, 2020
%  Desc:         Computes sum of squares for T1D system given options
%                (waveon, type, wICs) and parameters (params)

function ss = T1Dss_makeDists(params, obsPts, waveon, type, wICs, datalen) 
% Separate pre/post onset data
obsPts_pre = obsPts(1:25);
obsPts_post = obsPts(26:end);

% Solve the ODE, find residual for glucose
[~, Ypre, ~, Ypost] = T1Dfun_makeDists(params, waveon, type, wICs, datalen); 
res_pre = (Ypre(:,6) - obsPts_pre).*(Ypre(:,6) - obsPts_pre);
res_post = (Ypost(2:end,6) - obsPts_post).*(Ypost(2:end,6) - obsPts_post);
 
% Square the error values and sum
ss = sum(res_pre) + sum(res_post);    
end