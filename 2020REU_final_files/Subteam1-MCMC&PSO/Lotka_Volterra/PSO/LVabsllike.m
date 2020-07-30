%% Lotka-Volterra Abs Log Likelihood Func
%  Author:      Christina Catlett
%  Date:        June 9, 2020     
%  Desc:        log10 Gaussian likelihood function for LV model; used as
%               objective func

function absllike = LVabsllike(params, data)

% Separate time series (first col) and collected data
time   = 0:(length(data(:,1)) - 1); 
ydata  = data(:,2:end);

% Solve ODE for current params
modelOutput = LVfun(time, params, ydata);

% Calculate log10 normal prob of residuals
normH = log10(normpdf(modelOutput(:,1) - ydata(:,1), 0, sqrt(params(5)))); 
normL = log10(normpdf(modelOutput(:,2) - ydata(:,2), 0, sqrt(params(6))));

% Sum over all datapoints for final llike
llike = sum(normH + normL);

% Take abs val -> to minimize in PSO
absllike = abs(llike);
end