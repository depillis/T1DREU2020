%% Shift and average Li data sets (consolidation of mouse data)
%  Author:       Edited by C. Catlett & M. Watanabe, from D. Shenker 2020
%  Date:         June 2020
%  Desc:         Shift Li data to onset time point and average data. Used
%                to create a single data set
function mean_data = lietal_meandat(mousetype)
counts = zeros(1, 39); %vector to hold number of readings for each time point
vq_sum = zeros(1,39); %vector to hold sum of readings

if mousetype == 1 % acute mouse: all except 1 and 5
for i = 1:11 %loop over all mice 
    i_str = int2str(i);
    file = strcat('dat',i_str,'.csv'); %create file name
    X = readtable(file);
    X = X{:,:};
    X = flip(X);
    X = X'; %set up table as wide
    [~, col] = size(X);
    X_new = zeros(2, col);
    X_new(:, 2:col + 1) = X;
    X_new(:,1) = [1;100]; %add initial condition 
    min_time = round(X_new(1,2)); %minimum time available
    max_time = round(X_new(1,end)); %maximum time available
    for t = min_time:max_time
        counts(1,t) = counts(1,t) + 1; %add to counts
    end
    
    if i ~= 1 && i ~= 5
    vq = interp1(X_new(1,:), X_new(2,:), [min_time:1:max_time], 'pchip', 'extrap'); %interpolate
    vq_update(11:min_time - 1) = zeros(1, min_time - 11); %set zeros for early values
    vq_update(max_time + 1:39) = zeros(1, 39 - max_time); %set zeros for later values
    vq_update(min_time:max_time) = vq; 
    vq_sum = vq_sum + vq_update; %add to sum of glucose levels
    end
end

mean_data = vq_sum(1,11:33)./counts(1,11:33); %find mean by dividing sum by counts for each entry
% figure(1)
% p1 = plot([11:1:33], mean_data, '-o'); %plot the mean data

else
    for i = 1:11 %loop over all mice 
    i_str = int2str(i);
    file = strcat('dat',i_str,'.csv'); %create file name
    X = readtable(file);
    X = X{:,:};
    X = flip(X);
    X = X'; %set up table as wide
    [~, col] = size(X);
    X_new = zeros(2, col);
    X_new(:, 2:col + 1) = X;
    X_new(:,1) = [1;100]; %add initial condition 
    min_time = round(X_new(1,2)); %minimum time available
    max_time = round(X_new(1,end)); %maximum time available
    for t = min_time:max_time
        counts(1,t) = counts(1,t) + 1; %add to counts
    end
    
    vq = interp1(X_new(1,:), X_new(2,:), [min_time:1:max_time], 'pchip', 'extrap'); %interpolate
    vq_update(11:min_time - 1) = zeros(1, min_time - 11); %set zeros for early values
    vq_update(max_time + 1:39) = zeros(1, 39 - max_time); %set zeros for later values
    vq_update(min_time:max_time) = vq; 
    vq_sum = vq_sum + vq_update; %add to sum of glucose levels
    
    end

    mean_data = vq_sum(1,11:33)./counts(1,11:33); %find mean by dividing sum by counts for each entry
% figure(1)
% p1 = plot([11:1:33], mean_data, '-o'); %plot the mean data

end
end

