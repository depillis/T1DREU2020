%% Simulate predator-prey data
% Authors:      C. Catlett
%
% Date:         July 2020
%
% Descr:        Script to simulate hare and lynx data. The Lotka-Volterra
%               ODE model is evaluated and then 10% Guassian noise is added
%               to each data point.
               
function simData = makeSimData(params)
% Load datasets
load('HaresLynxData.mat')    % Loads as 'Lotka_Volterra_Data'
data = Lotka_Volterra_Data;
% Select time, data from file
time  = 0:(length(data(:,1)) - 1); % Change time from dates to years from 0
ydata = data(:,2:end);
% Solve model
y0 = ydata(1,:);
[tMod,yMod] = ode45(@(t,y) lotkaVolterrasys(t, y, params), time, y0);
% Determine noise (10%  Gaussian)
noiseSigmaH = .1*yMod(:,1);
noiseSigmaL = .1*yMod(:,2);
noiseH = noiseSigmaH' .* randn(1, length(time));
noiseL = noiseSigmaL' .* randn(1, length(time));
% Add noise to model predictions
wNoiseH = yMod(:,1) + noiseH';
wNoiseL = yMod(:,2) + noiseL';
% Create, save simulated data
simData = [tMod wNoiseH wNoiseL];
save('makesimData.mat', 'simData');
