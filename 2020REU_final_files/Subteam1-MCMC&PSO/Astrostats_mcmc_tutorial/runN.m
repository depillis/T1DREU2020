%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: Maya Watanabe & Christina Catlett
% Date: 5/25/2020
% Description: Plotting potential linear fittings after MCMC

function randlines = runN(N,a0, a0den, a1, a1den, ysig, ysigden, obsdata_x, obsdata_y)
% Sampling random values from posterior PDFs
a0rands = randVals(a0den,a0, N);
a1rands = randVals(a1den, a1, N);
ysigrands = randVals(ysigden,ysig, N);

randlines = zeros(N, 3);
% Selecting values of parameters from sample
for i=1:N
randlines(i, 1) = a0rands(randi(numel(a0rands)));
randlines(i, 2) = a1rands(randi(numel(a1rands)));
randlines(i, 3) = ysigrands(randi(numel(ysigrands)));
end 

% plot sampled solutions (lines) to the linear model
%   overtop the simulated data (obsdata)
scatter(obsdata_x, obsdata_y)
hold on
xaxis=0:15; % x axis range
for j=1:N
    y=randlines(j,1) + randlines(j,2)*xaxis + randlines(j,3);
    hold on
    plot(xaxis,y)
end

end