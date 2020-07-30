%% Function to find onset of diabetes distributions for Li et al data
function [meandat,onset_dist] = onset_time_dist_full
%TO GET INITIAL CONDITIONS, RUN ODE UNTIL GLUCOSE REACHES 250. AT THIS
%POINT, WE HAVE OUR INITIAL STATE CONDITIONS.
wave = 0;
allParam_vals;
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
T = 280;
tspan = [1 T];
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
sol = ode15s(@(t, y) T1D_ODE(t, y, f1n, f2n, wave), tspan, InitialState, options); %Use ODE solver
for t = 1:T %Iterate over times
    curr_sol = deval(sol, t); %get solution for that time
    curr_sol(6);
    if curr_sol(6) >= 250 %Check if glucose exceeds threshold
        initial_state = curr_sol; %Set initial conditions to that value
        break; %break out of loop
    end
end

newdat = zeros(2,1);

for i = 1:11 

    i_str = int2str(i);
    file = strcat('dat',i_str,'.csv'); %create file name
    X = readtable(file);
    X = X{:,:};
    X = flip(X);
    X = X'; %set up table as wide
    
    %FIND FIRST READING OVER 250
    numReadings = size(X,2);
    for j = 1:numReadings
        currentReading = X(2,j);
        if currentReading >= 250
            firstTime = j;
            break;
        end
    end
    
    topReading = X(2,firstTime); %Reading once reach 250
    bottomReading = X(2,firstTime - 1); %Reading before reaching 250
    topTime = X(1,firstTime); %Time when over 250
    bottomTime = X(1,firstTime - 1); %Time before being over 250
    slope = (topReading - bottomReading) / (topTime - bottomTime); %slop between before and after 250
    delta_t = (250 - bottomReading) / slope; %Change in time to reach 250
    
    X_new = zeros(2, numReadings - firstTime + 2); %Create new X vector
    X_new(:, 2:end) = X(:, firstTime:end); %Shift everything over a spot
    X_new(:, 1) = [bottomTime + delta_t; 250]; %Add initial time point
    
    %INTERPOLATE VALUES TO FILL IN X_NEW
    min_time = round(X_new(1,1)); %minimum time available
    max_time = round(X_new(1,end)); %maximum time available
    vq = interp1(X_new(1,:), X_new(2,:), [min_time:1:max_time], 'pchip', 'extrap'); %interpolate
    if i ~= 1 || i ~=5 
        newdat = horzcat(newdat, [min_time:1:max_time;vq]);
    end
       
    % collect time at 250
    threshold(i) = bottomTime+delta_t;    
end

newdat = sort(newdat');
newdat = newdat(2:end, :);

meandat = zeros(2,1);
for j = newdat(1,1):newdat(end,1)
indices = find(newdat(:,1) == j);
result = mean(newdat(indices, 2));
meandat = horzcat(meandat, [j; result]);
end
meandat=meandat(:,2:end);
%meandat = meandat(:, 8:end);
% Remove mice 1 and 5 (outliers)
threshold(1)=[];
threshold(5)=[];

% fit a distribution
% histogram(threshold) % plot onset times
onset_dist=fitdist(threshold', 'normal');
% 
%plot distribution
% x = 0:1:39;
% y = pdf(onset_dist,x');
% plot(x,y)
end


