%% Lotka-Volterra Sum of Squares Obj Func
%  Author:      Christina Catlett
%  Date:        June 9, 2020     
%  Desc:        SS objective function for LV model
%               (used in fitLV.m)

function ss = LVss(params, data)

% Separate time series (first col) and collected data
time   = 0:(length(data(:,1)) - 1); 
ydata  = data(:,2:end);

% Solve ODE for current params
modelOutput = LVfun(time, params, ydata);

% Calculate difference squared for columns
Hres = (modelOutput(:,1) - ydata(:,1)).*(modelOutput(:,1) - ydata(:,1));
Lres = (modelOutput(:,2) - ydata(:,2)).*(modelOutput(:,2) - ydata(:,2));


% Sum rows, columns
ss = sum(Lres+Hres);
end