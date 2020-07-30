%% Compute standard error of the estimate to judge goodness of DRAM model
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         sigma_est = sqrt(SSE/num of pairs of scores)
%                SSE = (acutal - predicted)^2

function sigma_est = error_est(observed, predicted)
n = length(observed);
sse = sum((observed-predicted).^2);

sigma_est = sqrt(sse/n);
end