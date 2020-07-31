%% Shift and average Li et al mice
% determine diabetes onset time
function [threshold, onsetdist] = onset_dist_Li(mousetype)

threshold = zeros(1,1); % vector to hold diabetes onset times

if mousetype == 1 % acute
for n = 1:11 %loop over mice 3 to 11
    n_str = int2str(n);
    file = strcat('dat',n_str,'.csv'); %create file name
    X = readtable(file);
    X = X{:,:};
    X = flip(X);
    %X = X'; %set up table as wide   
  % FIND SECOND CONSECUTIVE READING OVER 240
  numReadings = length(X);
    for j = 1:numReadings
        currentReading = X(j,2);
        if currentReading >= 240
            if j ~= numReadings
                nextReading = X(j+1,2);
                if nextReading >= 240
                    firstTime = j+1;
                    break;
                end
            end 
        end
    end
    

        if n ~= 1 && n ~= 5
            min_time = X(firstTime, 2); %minimum time available
            max_time = X(end, 2); %maximum time available
            vq = interp1(X(:,1), X(:,2), [min_time:1:max_time], 'pchip', 'extrap'); %interpolate
   
            threshold = [threshold X(firstTime, 1)]; % collect onset times
        end
end

    threshold(1) = [];
    onsetdist = fitdist(threshold','lognormal');

%     figure(1); clf
%     histogram(threshold)
%     xlabel('weeks')
%     title('Histogram of diabetes onset times')
% 
%     figure(2); clf
%     x = 0:1:60;
%     y = pdf(onsetdist, x');
%     plot(x,y)
%     xlabel('weeks')
%     title('Distribution of diabetes onset times')

else
    for n = 1:11 %loop over mice 3 to 11
    n_str = int2str(n);
    file = strcat('dat',n_str,'.csv'); %create file name
    X = readtable(file);
    X = X{:,:};
    X = flip(X);
    %X = X'; %set up table as wide
    
          
  % FIND SECOND CONSECUTIVE READING OVER 240
  numReadings = length(X);
    for j = 1:numReadings
        currentReading = X(j,2);
        if currentReading >= 240
            if j ~= numReadings
                nextReading = X(j+1,2);
                if nextReading >= 240
                    firstTime = j+1;
                    break;
                end
            end 
        end
    end
    

        if n ~= 1 && n ~= 5
            min_time = X(firstTime, 2); %minimum time available
            max_time = X(end, 2); %maximum time available
            vq = interp1(X(:,1), X(:,2), [min_time:1:max_time], 'pchip', 'extrap'); %interpolate
   
            threshold = [threshold X(firstTime, 1)]; % collect onset times
        end
end

%     threshold(1) = [];
%     onsetdist = fitdist(threshold','lognormal');
% 
%     figure(1); clf
%     histogram(threshold)
%     xlabel('weeks')
%     title('Histogram of diabetes onset times')
% 
%     figure(2); clf
%     x = 0:1:60;
%     y = pdf(onsetdist, x');
%     plot(x,y)
%     xlabel('weeks')
%     title('Distribution of diabetes onset times')
end
    

end
