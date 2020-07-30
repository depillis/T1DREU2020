%% Create upper and lower bounds for T1D parameters (41) 
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         95% confidence bounds for parameters from table 
%                (D. Shenker & R. Wander, June 2020). Returns a vector 
%                containing baseline, upper, and lower bounds of each
%                parameter.
function [lb, ub, base] = makeparambounds
% Read in Daniel/Rachel's table
T = readtable('Parameter_info.xlsx');

    % extract parameter variance and ranges
    % [variance, baseline, lower bnd, upper bnd]
    temp = table(T.Variance, T.efast_base, T.efast_lower, T.efast_upper);
    info = table2array(temp(:,2:end));
    
    % set up to collect ranges and baseline values
    base = zeros(1,1);
    lb = zeros(1,1);
    ub = zeros(1,1);
    
    % find 95% confidence interval for parameter variability (normal
    % assumption)
    % 95% CI = +/- 2 standard deviations
    for i = 1:length(info)
        sd = sqrt(info(i,1));
        base = vertcat(base, info(i,2));
        if info(i,2)-(2*sd) >= 0 % check non negativity
            lb = vertcat(lb, info(i,2)-(2*sd));
        else
            lb = vertcat(lb, info(i,2));
        end
        ub = vertcat(ub, info(i,2)+(2*sd));
    end
    base(1) = [];
    lb(1) = [];
    ub(1) = [];
    


end