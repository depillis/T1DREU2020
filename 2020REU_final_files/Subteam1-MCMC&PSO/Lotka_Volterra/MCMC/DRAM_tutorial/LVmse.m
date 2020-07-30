%% Mean square error
% Authors:      M. Watanabe
%
% Date:         July 2020
%
% Descr:        Script to calculate mean square error

function [mse_prey mse_pred] = LVmse(params, data)% Sum-of-squares function

time   = data(:,1); 
ydata  = data(:,2:end);

ymodel = lotkaVolterrafun(time,params,ydata);
%ss = sum(sum((ymodel - ydata).^2));
ss = sum((ymodel-ydata).^2);
mse_prey = ss(:,1)/length(data);
mse_pred = ss(:,2)/length(data);
end
