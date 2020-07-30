%% Function to find onset of diabetes distributions for Mathews et al data
%  Author:       C. Catlett, M. Watanabe, using code from Shtylla et al
%                2019 
%  Date:         June 2020
%  Desc:         Collect onset of diabetes times and fit distributions.
%                Used in DRAM algorithm to allow Mathews et al to be in 
%                absolute time.

function [onset_distX, onset_distY] = onset_time_dist_mathews

% Load shifted data (Shtylla et al 2019)
generate_onset_times;

% Define shifted data
X = shifted_progressive;
Y = shifted_acute;

% Determine columns with all NaN (no data)
colstodeleteX = removeCols(X);
colstodeleteY = removeCols(Y);

thresholdX = zeros(1,1); % Initialize onset time vector
thresholdY = zeros(1,1); 

% progressive onset time dist
for n = 2:2:length(X(1,:)) % For each mouse...
    
 % Skip column if all NaN
if ~ismember(n, colstodeleteX)
          
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
  
    thresholdX = [thresholdX X(firstTime, n-1)];
end

% acute onset dist
for n = 2:2:length(Y(1,:)) % For each mouse...
    
 % Skip column if all NaN
if ~ismember(n, colstodeleteY)
          
  % Find second consecutive reading over 240
  numReadings = 25;
    for j = 1:numReadings
        currentReading = Y(j,n);
        % Check next reading if applicable
        if currentReading >= 240 && j ~= 25
            nextReading = Y(j+1,n);
            if nextReading >= 240
               firstTime = j+1;
               break;
            end
        end
    end
end
  
    thresholdY = [thresholdY Y(firstTime, n-1)];
end

% Remove 0 inserted for initialization 
thresholdX(1) = [];
thresholdY(1) = [];

onset_distX = fitdist(thresholdX', 'lognormal');
onset_distY = fitdist(thresholdY', 'lognormal');

% PLOTTING
% % Historgram of onset times IN DAYS
% figure (1); clf
% histogram(thresholdX)
% xlabel('days')
% title('Histogram of diabetes onset times: Progressive')
% figure(2); clf
% histogram(thresholdY)
% xlabel('days')
% title('Histogram of diabetes onset times: Acute')
% 
% % % Plot fitted log-normal
% figure (3); clf
% onset_distX = fitdist(thresholdX', 'lognormal');
% x = 0:1:240;
% y = pdf(onset_distX,x');
% plot(x,y)
% xlabel('days')
% title('Distribution of diabetes onset times: Progressive')
% 
% 
% figure (4); clf
% onset_distY = fitdist(thresholdY', 'lognormal');
% xx = 0:1:240;
% yy = pdf(onset_distY, xx');
% plot(xx, yy)
% xlabel('days')
% title('Distribution of diabetes onset times: Acute')
end


