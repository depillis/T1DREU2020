function [mse, params] = fitInitialParams(params, data)
% Extract the initial guess for your parameters (must be in vector form)
for i = 1:length(params)
    initGuess(i) = params{i}{2}; % value of parameter
    lb(i) = params{i}{3}; % lower bound of uniform prior
    ub(i) = params{i}{4}; % upper bound of uniform prior
end

%Run a quick fminsearch to get initial parameter and s2 guesses for DRAM
if length(data(:,1)) > length(params)
    fun = @(params) lotkaVolterrass(params, data);
    opt = optimset('Display', 'off');
    [theta, ss0] = fmincon(fun, initGuess, [], [], [], [], lb, ub, [], opt);
    mse = ss0/(length(data(:,1))-length(params));
    
    %Now replace params starting values with fminsearch output
    for i = 1:length(params)
        params{i}{2} = theta(i);
    end
else
    mse = 1; %THIS MAY NEED TO BE CHANGED! Normal guess for mse breaks down! Need more data pts than parameters
end