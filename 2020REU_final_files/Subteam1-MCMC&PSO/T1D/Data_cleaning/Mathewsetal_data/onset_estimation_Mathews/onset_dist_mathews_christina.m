%% Function to find onset of diabetes distributions for Mathews et al data
function [mu, sig, mulb, muub, siglb, sigub] = onset_dist_mathews_christina(type)

% Load data
generate_onset_times;

% Define shifted data
progressive = shifted_progressive;
acute = shifted_acute;

% Select data based on func input
switch type
    case 'prog'
        X = progressive;
    case 'acute'
        X = acute;
end

% Determine columns with all NaN (no data)
colstodelete = removeCols(X);

threshold = zeros(1,1); % Initialize onset time vector
for n = 2:2:length(X(1,:)) % For each mouse...
    
 % Skip column if all NaN
if ~ismember(n, colstodelete)
          
  % Find second consecutive reading over 240
  numReadings = 25;
    for j = 1:numReadings
        currentReading = X(j,n);
        % Check next reading if applicable
        if currentReading >= 240 && j ~= 25
            nextReading = X(j+1,n);
            if nextReading >= 240
               firstTime = j+1;
               break;
            end
        end
    end
end
  
    threshold = [threshold X(firstTime, n-1)/7];
end

% Remove 0 inserted for initialization 
threshold(1) = [];
onset_dist = fitdist(threshold', 'lognormal');

mu = onset_dist.mu;
sig = onset_dist.sigma;

ci = paramci(onset_dist);
mulb = ci(1,1);
muub = ci(2,1);
siglb = ci(1,2);
sigub = ci(2,2);
end
