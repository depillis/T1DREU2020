%% Create upper and lower bounds for initial conditions (T1D)
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         Create ranges for T1D model initial conditions: resting
%                macrophage, beta cells, glucose, insulin

function [lb, ub, base] = makeICbounds(percent)
% creates ranges for T1D initial conditions
    base = [4.77*10^5, 300, 100, 10];
    
    lb = [base(1)*(1-percent), 300*(1-percent), 100*(1-percent), 10*(1-percent)];
    
    ub = [base(1)*(1+percent), base(2)*(1+percent), base(3)*(1+percent), base(4)*(1+percent)];

end