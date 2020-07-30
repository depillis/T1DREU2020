function Tmse = PSOmse(params)
waveon = 1;
type = 'NOD';
wICs = 1;
% Load dataset
load('liMouse6.mat');

% Select time series, data of interest
t = dat6.VarName1.*7;
y = dat6.VarName2;

[~, Y] = T1Dfun_noDist(params, waveon, type, wICs, t);
modPred = Y(:, 6);
res = (y-modPred).*(y-modPred);
ss = sum(res);
mse = ss./length(t)
end
