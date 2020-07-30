%DOES PARAMETER ESTIMATION FOR 10 MICE FROM LI ET AL DATASET. FIRST SHIFTS
%MICE OVER TO TIME ZERO WHERE THE GLUCOSE VALUE AT TIME ZERO IS FIRST AFTER
%T1D ONSET HAS BEGUN


%TO GET INITIAL CONDITIONS, RUN ODE UNTIL GLUCOSE REACHES 250. AT THIS
%POINT, WE HAVE OUR INITIAL STATE CONDITIONS.
wave = 0;
fig_2_3_4_5_Parameters;
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
T = 280;
tspan = [1 T];
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, truepar), tspan, InitialState, options); %Use ODE solver
for t = 1:T %Iterate over times
    curr_sol = deval(sol, t); %get solution for that time
    curr_sol(6)
    if curr_sol(6) >= 250 %Check if glucose exceeds threshold
        initial_state = curr_sol; %Set initial conditions to that value
        break; %break out of loop
    end
end



InitialState = initial_state; %Update initial state values

InitialParamGuess = [10^-5; 0.0334; sqrt(2000); .72; .194; aEmday/5000; aEmday/5000]; %Set initial parameter guesses

%--------RUN LOOP FOR EVERY MOUSE ---------------------------------------%
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
    vq
    
    
    
end


Take MEans
Run again

