% This function creates the parameter structure needed to initiate the DRAM
% algorithm. Priors (uniform bounds) can be changed as needed to ensure
% converging chains.

function params = paramStruct

% All parameters are constrained to be positive. The initial
% concentrations are also unknown and are treated as extra parameters.

% ParamName, starting value, uniform prior bounds, or normally distributed
%priors with Inf indicator

params = {
    {'mumax', 0.5,  0}
    {'rhoa',  0.03, 0}
    {'rhoz',  0.1,  0}
    {'k',     10,   0}
    {'alpha', 0.02, 0}
    {'th',    1.14, 0, Inf, 1.14, 0.2}  % N(1.14, 0.2^2){th>0} prior
% initial values for the model states
    {'A0', 0.77, 0, Inf, 0.77, 2 }
    {'Z0', 1.3,  0, Inf, 1.3,  2 }
    {'P0', 10,   0, Inf, 10,   2 }
    };
end