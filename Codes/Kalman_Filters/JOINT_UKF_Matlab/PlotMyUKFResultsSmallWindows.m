% LdeP April 6, 2020 - Script to plot results of UKF run - just to keep
% code cleaner
%
% Must first run AlbersUKF_ParState.m

%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%

%LdeP Plot states and estimates
figure();
% Plot TRUE solution
%p = plot(timeVector,xTrueGlucosePerDL,'-r');
p = plot(timeSteps,xTrueGlucosePerDL,'-b');
p.LineWidth = 3;
% Plot UKF solution
hold on;
p = plot(Trandstepvec,xCorrectedUKF(:,3),'-r',Trandstepvec,yMeas(:,1),'-om');
p(1).LineWidth=3; %LdeP increase line width of UKF guess
p(2).LineWidth=2; %LdeP increase line width of Measurement
p(2).MarkerSize=12;
ppar1 = plot(Trandstepvec,xCorrectedUKF(:,7),'-r',Trandstepvec,yMeas(:,2),'-or');
ppar2 = plot(Trandstepvec,xCorrectedUKF(:,8),'-g',Trandstepvec,yMeas(:,3),'-og');
ppar3 = plot(Trandstepvec,xCorrectedUKF(:,9),'-k',Trandstepvec,yMeas(:,4),'-ok');
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
legend('True','UKF state estimate','State Measured','UKF E','E meas','UKF Vi','Vi meas', 'UKF ti','ti meas');
xlabel('Time [min]', 'FontSize', 20);
ylabel('Glucose [mg/dL]', 'FontSize', 20);
%LdeP comparing outputs of different measurment noise models
title('Albers Type 2 Diabetes Model - UKF fit to generated data: Additive Meas Noise');


%LdeP Plot error
figure();
p=plot(Trandstepvec, e(:,1), '-o',Trandstepvec,e(:,2),'-.', Trandstepvec,e(:,3),'-', Trandstepvec,e(:,4),'*-');
set(gca, 'LineWidth',2); %LdeP Make axes and boxes thicker
set(gca, 'FontSize', 20);
p(1).LineWidth=2; p(2).LineWidth=2; p(3).LineWidth=2; p(4).LineWidth=2;
legend('State Error','E error', 'Vi Error', 'ti Error');
xlabel('Time [min]', 'FontSize', 20);
ylabel('Residual (or innovation)', 'FontSize', 20);
title('Albers Type 2 Diabetes Model - Residuals: Additive Meas Noise');

%LdeP Plot error and estimates on same graph
figure();
% Plot TRUE solution
%p = plot(timeVector,xTrueGlucosePerDL,'-r');
p = plot(timeSteps,xTrueGlucosePerDL,'-r');
p.LineWidth = 3;
% Plot UKF solution and error
hold on;
p = plot(Trandstepvec,xCorrectedUKF(:,3),'-b',Trandstepvec,yMeas(:,1),'-om', Trandstepvec, e(:,1), '-ok');
p(1).LineWidth=3; %LdeP increase line width of UKF guess
p(2).LineWidth=2; %LdeP increase line width of Measurement
p(2).MarkerSize=12;
p(3).LineWidth=2; %LdeP increase line width of Measurement
p(3).MarkerSize=12;
set(gca, 'LineWidth',2); %LdeP Make axes and boxes thicker
set(gca, 'FontSize', 20);
legend('True','UKF estimate','Measured', 'Meas Err');
xlabel('Time [min]', 'FontSize', 20);
ylabel('Glucose [mg/dL]', 'FontSize', 20);
title('Albers Type 2 Diabetes Model - UKF fit to generated data: Additive Meas Noise and Error');
