%% Lotka-Volterra Sum of Squares
% N. Tania (May 30, 2018), edited by Christina Catlett (June 3, 2020)
%
% Sum of squares of model predictions, observed data

function ss = lotkaVolterrass(params, data)
% Sum-of-squares function

time   = [0:(length(data(:,1)) - 1)]; 
ydata  = data(:,2:end);

ymodel = lotkaVolterrafun(time,params,ydata);
ss = sum(sum((ymodel - ydata).^2));