%% Log Prior Function
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         Informative logprior function for mcmcrun (mcmcstat
%                library). Hard-coded as a lognormal prior, but this can be
%                updated to suit the model's needs. Default is a
%                uninformative prior.

function logprior = logprior(params,parammu,paramsig)
logprior = 0;
if ~isempty(parammu) && ~isempty(paramsig)
% for i = 1:length(params)
%    prior = normpdf(params(i), parammu(i), paramsig(i)); % Informative log priors for params: 
%    logprior = logprior + log(prior);
% end

% for i = 1:length(params)
%    prior = norpf(params(i), parammu(i), paramsig(i)); % Informative log priors for params: 
%    logprior = logprior + log(prior);
% end
lp = lognorpf(params, parammu, paramsig);
logprior = sum(lp);
% DEFAULT
% Uninformative prior: uniform between 1e5 and 1
else
    for i = 1:length(params)
   prior = normpdf(params(i), parammu(i), paramsig(i));
   logprior = logprior + log(prior);
    end
end
    
end