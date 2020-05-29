% Authors: LdeP
% Date: March 14, 2020
% Summary: UKF implementation of Albers study 2018

% Date: March 31, 2020
% Updates: LdeP Working to test alternately spaced sampling times

% Date: April 1, 2020 LdeP
% Using the "uneven data sample" general structure to allow for an even data sample.

% Updates: LdeP March 2020 - working to fix the large filter shift
% March 25, 2020 - filter shift fixed:
% LdeP: The problem is the scaling of the Glucose (3rd state)
% When we call the AlbersODE function, we need to rescale to
% compute total glucose, and then we need to scale back again
% to do the fits with glucose per dL so that we match the "data."
% Update: LdeP March 30, 2020 - experimenting to see what happens when we
% use the UKF implementation parameters given by Albers 2017 paper in
% supplement 1 (Albers2017S1Appendix.pdf).

% Including to get I_g line 31
% AlbersParamVals;

% LdeP Read in the Albers dual UKF algorithm parameters
% Rv is the Assumed process noise covariance matrix
% Rn is the Assumed measurement noise covariance matrix 
AlbersParamsForUKFAlgorithmFrom2017AppendixS1; 

% Initial conditions taken from Albers code
% I_p = 200
% I_i = 200
% G = 12000
% h_1 = 0.1
% h_2 = 0.2
% h_3 = 0.1
initialStateGuessODE = [200; 200; 12000; 0.1; 0.2; 0.1]; %LdeP This is the same set of initial values used in the Albers code.
%initialStateGuess = [200; 200; 3000; 0.1; 0.2; 0.1];
%timeFinal = 100; % t = 9 days worth of minutes = 12960 minutes
timeFinal=500; %LdeP Try a different length of time. Albers went to over 1200.


% LdeP To get randomly spaced sample times, change T below
T = 1; % [s] Filter sample time
timeVector = 0:T:timeFinal; %LdeP Fiddling with the time vector

% UNEVEN DATA SAMPLE - Uncommented in the Uneven data sample code.
% LdeP Trying out a time vector with random sampling times
% Downloaded randspace from MathWorks
% Note that you may get fewer sample points than you specify if your
% MaxTimeStep is too large compared to timeFinal - randspace will just cut
% off the entries that get too large. So you need to be careful that
% dimensions work out when you use randspace. There is likely a better way
% to generate a monotonically increasing list of numbers with random
% intervals between them.
% NoSamples=100;
% MaxTimeStep=8;
% Trandstepvec=ceil(randspace(0,NoSamples,timeFinal,[1,MaxTimeStep]));

% EVEN DATA SAMPLE - allow Trandstepvec simply to equal the timeVector
Trandstepvec = timeVector+1; %LdeP Add 1 so the timestep vector starts at 1
NumSamples = numel(Trandstepvec); %LdeP

% I = [[1:timeFinal]' I_g*ones(timeFinal, 1)];
% [] space holder for special options
[timeSteps,xTrue]=ode45(@AlbersODE,timeVector,initialStateGuessODE);

% System has 10L of glucose = 100dL
% We compute total mg of glucose in system
% but plot mg dL^(-1)
% So we have to divide glucose by 100 before plotting
xTrueGlucosePerDL = xTrue(:,3)./100; %LdeP Save glucose measure in ml/dL


% Moved construction of filter to input different initial state guesses
% This allows UKF and yMeas to start at same initial values
initialStateGuessUKF = [200; 200; 120; 0.1; 0.2; 0.1]; %LdeP State vector with glucose scaled to mg/dL

% ORIGINAL MEASUREMENT NOISE (looks additive, called as non-additive) - uncomment below
% Construct the filter
% LdeP: Note: With AlbersNoiseFcn_orig, the noise function is called as though it is not
% additive, but the form of the noise is yk = xk(3)+vk (which is
% additive). Results using this form look much like the "additive noise" results
% ukf = unscentedKalmanFilter(...
%     @AlbersStateFcn,... % State transition function
%     @AlbersNoiseFcn_orig,... % Measurement function
%     initialStateGuess,...
%     'HasAdditiveMeasurementNoise',false);

% NON-ADDITIVE MEASUREMENT NOISE - uncomment below
%LdeP experimenting with additive an non-additive noise
%LdeP The way we have written it, non-aditive measurement noise will produce a
%much more smooth curve.
% ukf = unscentedKalmanFilter(...
%     @AlbersStateFcn,... % State transition function
%     @AlbersNoiseFcn_NonAdd,... % Measurement function
%     initialStateGuessUKF,...
%     'HasAdditiveMeasurementNoise',false);

% ADDITIVE MEASUREMENT NOISE - uncomment below
%LdeP The way we have written it, additive measurement noise will produce a
%much more jagged curve.
ukf = unscentedKalmanFilter(...
    @AlbersStateFcn,... % State transition function
    @AlbersNoiseFcn_Additive,... % Measurement function
    initialStateGuessUKF,...
    'HasAdditiveMeasurementNoise',true);

% LdeP: Making R large generates more noisy measurements - and estimates will
% favor following the model rather than the measured data.
% LdeP: smaller R -> closer estimate will match measured data
% R = 0.2; % Variance of the measurement noise v[k]
% R = 0.02;
% R = 1.0;
% R = .1;
% sqrtR=1.5; %LdeP Standard Deviation of measurement noise
% sqrtR = 5; %LdeP Supposing that measurments of blood glucose will deviate by + or -5 per dL.

% LdeP March 30, 2020 - see what happens when we use Albers UKF parameters:
% Rn - process noise covariance
% Rv - measurement noise covariance
sqrtR = sqrt(Rn); % LdeP Rn = measurment noise covariance from Albers2017
ukf.MeasurementNoise = sqrtR^2; %LdeP Variance of measurement noise.

% LdeP: Making Process Noise too large won't allow sufficient smoothing
% The algorithm will think the "noisy" measurements are correct.
% small noise values -> large gap between measured and UKF estimate

% LdeP: March 30, 2020 - experimenting with using the process noise used in
% Albers2017, Rv (a 6x6 matrix). When we add parameters, we probably make
% the parameter process noise zero. Note, however, that the process noise
% from Albers is really large, so that will essentially make the UKF follow
% the measured data, and not the ODE model. The smaller the process noise,
% the closer to the smoothed ODE model we will be.
% But larger noise values ->  Our UKF estimate follows the "noisy" measure and won't smooth as much.
% ukf.ProcessNoise = diag([0.02 0.1 0.04 0.2 0.5 0.01]);
ukf.ProcessNoise = sqrt(Rv)/1000; %LdeP Making this smaller to adhere to the ODE model more than to the random data


% rng(1); % Fix the random number generator for reproducible results
rng('default'); % LdeP Matlab said use this
% xTrue(:,1) only pulls out third vector in matrix
% yTrue = xTrue(:,3); % Extract glucose, G
yTrue = xTrueGlucosePerDL; % Glucose is the only state we currently measure, G 

% LdeP Now generate some sample measurements
% yMeas = yTrue + some noise
% Noise = StandardDeviation*(random value between -1 and 1)
% sqrt(R): Standard deviation of noise
% LdeP If we set sqrtR (above) to be standard deviation of measurement noise,
% then sqrtR^2 is the variance of measurement noise, and we add that to our
% measure yMeas:
% NumberOfSamplePoints = size(yTrue); %LdeP just compute number of sample points once
% yMeas = yTrue + (sqrtR*randn(NumberOfSamplePoints)); %This follows the Matlab help example

%LdeP Trying randomly spaced samples
% NumberOfSamplePoints = NumSamples; %LdeP Prespecified number of sample points
% LdeP Trandstepvec is the vector of randomly spaced sample times
% LdeP 
yMeas = yTrue(Trandstepvec) + (sqrtR*randn(NumSamples,1)); %This follows the Matlab help example

%LdeP Create vector of time step jumps
DelT = Trandstepvec(2:end)-Trandstepvec(1:end-1);
DelT(NumSamples)=(mean(DelT(1:NumSamples-1))); %LdeP Fake final time step to make DelT correct dimension
for k=1:numel(yMeas)
    % Let k denote the current time.
    % Residuals (or innovations): Measured output - Predicted output
    % ukf.State is x[k|k-1] at this point
    e(k) = yMeas(k) - AlbersMeasFcn(ukf.State);
    
    %LdeP: Try predicting first, and then correcting
    % Note: It doesn't seem to matter in this example
    % whether we predict then correct or correct then predict.
    % Predict the states at time step k. 
    %[xPredictedUKF(k,:), PPredictedUKF(k,:,:)]=predict(ukf);
    % LdeP - send the of the time step (for data sampling) through predict
%     DelT = Trandstepvec(k+1)-Trandstepvec(k);
    [xPredictedUKF(k,:), PPredictedUKF(k,:,:)]=predict(ukf,DelT(k));

    
    % Incorporate measurements at time k into state estimates using
    % the "correct" command. This updates the State and StateCovariance
    % properties of the filter to contain x[k|k] and P[k|k]. These values
    % are also produced as the output of the "correct" command.
    [xCorrectedUKF(k,:), PCorrected(k,:,:)] = correct(ukf,yMeas(k));
    
    % Predict the states at next time step, k+1. This updates the State and
    % StateCovariance properties of the filter to contain x[k+1|k] and
    % P[k+1|k]. These will be utilized by the filter at the next time step.
    % predict(ukf);
    %[xPredictedUKF(k,:), PPredictedUKF(k,:,:)]=predict(ukf); %LdeP
    
end

%LdeP Update plotting 
figure();
% p = plot(timeVector,xTrueGlucosePerDL,'-r', timeVector,xCorrectedUKF(:,3), '-b',...
%     timeVector,yMeas(:), '-*m');
p = plot(timeVector,xTrueGlucosePerDL,'-r');
p.LineWidth = 3;
hold on;
p = plot(Trandstepvec,xCorrectedUKF(:,3),'-b',Trandstepvec,yMeas(:),'-om');
p(1).LineWidth=3; %LdeP increase line width of UKF guess
p(2).LineWidth=2; %LdeP increase line width of Measurement
p(2).MarkerSize=12;
set(gca, 'LineWidth',2); %LdeP Make axes and boxes thicker
set(gca, 'FontSize', 20);
legend('True','UKF estimate','Measured');
xlabel('Time [min]', 'FontSize', 20);
ylabel('Glucose [mg/dL]', 'FontSize', 20);
%LdeP comparing outputs of different measurment noise models
%title('Albers Type 2 Diabetes Model - UKF fit to generated data');
%title('Albers Type 2 Diabetes Model - UKF fit to generated data: Original Meas Noise Model');
% title('Albers Type 2 Diabetes Model - UKF fit to generated data: Non-Additive Meas Noise');
title('Albers Type 2 Diabetes Model - UKF fit to generated data: Additive Meas Noise');


figure();
% plot(timeVector, e, '-o'); %LdeP replace with sample size e
plot(Trandstepvec, e, '-o');
set(gca, 'LineWidth',2); %LdeP Make axes and boxes thicker
set(gca, 'FontSize', 20);
xlabel('Time [min]', 'FontSize', 20);
ylabel('Residual (or innovation)', 'FontSize', 20);
%title('Albers Type 2 Diabetes Model - Residuals: Original Meas Noise Model');
% title('Albers Type 2 Diabetes Model - Residuals: Non-Additive Meas Noise');
title('Albers Type 2 Diabetes Model - Residuals: Additive Meas Noise');

%LdeP Plot error and estimates on same graph
figure();
p = plot(timeVector,xTrueGlucosePerDL,'-r');
p.LineWidth = 3;
hold on;
p = plot(Trandstepvec,xCorrectedUKF(:,3),'-b',Trandstepvec,yMeas(:),'-om', Trandstepvec, e, '-ok');
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
%LdeP comparing outputs of different measurment noise models
%title('Albers Type 2 Diabetes Model - UKF fit to generated data');
%title('Albers Type 2 Diabetes Model - UKF fit to generated data: Original Meas Noise Model');
% title('Albers Type 2 Diabetes Model - UKF fit to generated data: Non-Additive Meas Noise');
title('Albers Type 2 Diabetes Model - UKF fit to generated data: Additive Meas Noise and Error');