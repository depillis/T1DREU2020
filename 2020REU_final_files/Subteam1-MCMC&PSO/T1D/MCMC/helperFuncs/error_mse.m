%% Compute mean square error to judge goodness of DRAM model
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         mean square error

function mse = error_mse(observed, predicted)% Sum-of-squares function
n = length(observed); 
ss = sum((observed-predicted).^2);

mse = ss/n;



end