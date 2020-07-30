%% Making bounds for initial conditions (Mathews & Li Data)
%  Author:       Christina Catlett
%  Date:         June 19, 2020
%  Desc:         Create upper and lower bounds for initial conditions,
%                percentage of varibility allowed; etxracts initial
%                parameters from struct to array/matrix

function [lb, ub] = makeICbounds(struct, percent)
% Array of ICs from struct
ICarr = formatArray(struct);
% Make lb, ub
lb = (1-percent).*ICarr;
ub = (1+percent).*ICarr;