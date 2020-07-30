%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Author: Bjorn Gustavsson
% (https://www.mathworks.com/matlabcentral
%   /answers/464916-how-to-generate-random-sampling
%   -from-a-probability-distribution-function-generated
%   -by-kernel-density?s_tid=answers_rc1-2_p2_MLT), edited by Christina
%   Catlett & Maya Watanabe
% Date: 5/26/2020
% Description: Random Sampling according to estimated PDF

function randVals = randVals(x, y, n)
cdfKD = cumtrapz(x, y);
cdfKD = (cdfKD - cdfKD(1)) / (cdfKD(end)-cdfKD(1)); % normalizing the cdf to be between zero and 1 at the end of your intervall
Q = rand(n); 
randVals = interp1(cdfKD,x,Q,'pchip-or-linear');
end