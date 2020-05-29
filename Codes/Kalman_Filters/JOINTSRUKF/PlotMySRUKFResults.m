% LdeP April 13, 2020 - Updated to reflect this is being used with 
% Square Root UKF code.
% LdeP April 6, 2020 - Script to plot results of UKF run - just to keep
% code cleaner
% LdeP April 7, 2020 - Updated to work with ukf_ALBERSmain.m (square root
% UKF)
%
% Must first run ukf_ALBERSmain.m

%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%

%LdeP Plot states and estimates
figure('units','normalized','outerposition',[0 0 0.7 0.7]) %LdeP make figure window larger
% Plot TRUE solution
%p = plot(timeVector,xTrueGlucosePerDL,'-r');
p = plot(timeSteps,xTrueGlucosePerDL,'-b');
p.LineWidth = 3;
% Plot UKF solution
hold on;
xtime = [1:dt:Nsteps]; % LdeP Time steps for UKF states
p = plot(xtime,xV(3,:),'-r',xtime,zV(1,:),'-om');
p(1).LineWidth=3; %LdeP increase line width of UKF guess
p(2).LineWidth=2; %LdeP increase line width of Measurement
p(2).MarkerSize=12;
ppar1 = plot(xtime,xV(7,:),'-c',xtime,zV(2,:),'-oc');
ppar2 = plot(xtime,xV(8,:),'-g',xtime,zV(3,:),'-og');
ppar3 = plot(xtime,xV(9,:),'-k',xtime,zV(4,:),'-ok');
ppar1(1).LineWidth=3; %LdeP increase line width of UKF guess
ppar1(2).LineWidth=2; %LdeP increase line width of Measurement
ppar1(2).MarkerSize=12;
ppar2(1).LineWidth=3; %LdeP increase line width of UKF guess
ppar2(2).LineWidth=2; %LdeP increase line width of Measurement
ppar2(2).MarkerSize=12;
ppar3(1).LineWidth=3; %LdeP increase line width of UKF guess
ppar3(2).LineWidth=2; %LdeP increase line width of Measurement
ppar3(2).MarkerSize=12;
set(gca, 'LineWidth',2); %LdeP Make axes and boxes thicker
set(gca, 'FontSize', 20);
legend('True','SR-UKF state estimate','State Measured','SR-UKF E','E meas','SR-UKF Vi','Vi meas', 'SR-UKF ti','ti meas');
xlabel('Time [min]', 'FontSize', 20);
ylabel('Glucose [mg/dL]', 'FontSize', 20);
%LdeP comparing outputs of different measurment noise models
title('Albers Type 2 Diabetes Model - SR-UKF fit to generated data: Additive Meas Noise');


%LdeP Plot error
figure('units','normalized','outerposition',[0 0 0.7 0.7]) %LdeP make figure window larger
p=plot(xtime, e(1,:), '-o',xtime,e(2,:),'-.', xtime,e(3,:),'-', xtime,e(4,:),'*-');
set(gca, 'LineWidth',2); %LdeP Make axes and boxes thicker
set(gca, 'FontSize', 20);
p(1).LineWidth=2; p(2).LineWidth=2; p(3).LineWidth=2; p(4).LineWidth=2;
legend('State Error','E error', 'Vi Error', 'ti Error');
xlabel('Time [min]', 'FontSize', 20);
ylabel('Residual (or innovation)', 'FontSize', 20);
title('Albers Type 2 Diabetes Model - Residuals from SR-UKF: Additive Meas Noise');

%LdeP Plot error and estimates on same graph
figure('units','normalized','outerposition',[0 0 0.7 0.7]) %LdeP make figure window larger
% Plot TRUE solution
%p = plot(timeVector,xTrueGlucosePerDL,'-r');
p = plot(timeSteps,xTrueGlucosePerDL,'-b');
p.LineWidth = 3;
% Plot UKF solution and error
hold on;
p = plot(xtime,xV(3,:),'-c',xtime,zV(1,:),'-om', xtime, e(1,:), '-oy');
p(1).LineWidth=3; %LdeP increase line width of UKF guess
p(2).LineWidth=2; %LdeP increase line width of Measurement
p(2).MarkerSize=12;
p(3).LineWidth=2; %LdeP increase line width of Measurement
p(3).MarkerSize=12;
set(gca, 'LineWidth',2); %LdeP Make axes and boxes thicker
set(gca, 'FontSize', 20);
legend('True','SR-UKF estimate','Measured', 'Meas Err');
xlabel('Time [min]', 'FontSize', 20);
ylabel('Glucose [mg/dL]', 'FontSize', 20);
title('Albers Type 2 Diabetes Model - SR-UKF fit to generated data: Additive Meas Noise and Error');
