%Function that can be used to calculate change in parameter values across
%iterations


function [] = ParamChange(final_params, InitialParamGuess, numParams)
param_change = final_params - InitialParamGuess;
param_change_pct = (param_change ./ InitialParamGuess) * 100;

for y = 1:numParams
    y
    param_change_pct(y)
end
end

