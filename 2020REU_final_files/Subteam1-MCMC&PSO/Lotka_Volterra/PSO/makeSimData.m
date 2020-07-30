%% Simulate data for Lotka-Volterra
%  Author:      Christina Catlett
%  Date:        July 3, 2020
%  Desc:        Starting at y0 given in data, simulate data by solving
%               model for given params, adding specified percent of Gaussian noise
%               (expressed in decimal, 0-1)

function simData = makeSimData(params, percent)
% Load data
load('HaresLynxData.mat') % Loads as 'Lotka_Volterra_Data'
data = Lotka_Volterra_Data;

% Select time, data from file
time  = 0:(length(data(:,1)) - 1); % Change time from dates to years from 0
ydata = data(:,2:end);

% Solve model
y0 = ydata(1,:);
[tMod,yMod] = ode45(@(t,y) LVsys(t, y, params), time, y0);

% Determine noise (percent%  Gaussian)
noiseSigmaH = percent*yMod(:,1);
noiseSigmaL = percent*yMod(:,2);
noiseH = noiseSigmaH' .* randn(1, length(time));
noiseL = noiseSigmaL' .* randn(1, length(time));

% Add noise to model predictions
wNoiseH = yMod(:,1) + noiseH';
wNoiseL = yMod(:,2) + noiseL';

% Create, save simulated data
simData = [tMod wNoiseH wNoiseL];
save('makesimData.mat', 'simData');


