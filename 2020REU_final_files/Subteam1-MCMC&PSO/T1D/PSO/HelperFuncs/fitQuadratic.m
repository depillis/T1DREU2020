%% Toy example: fitting a quadratic function with PSO
%  Author:       Christina Catlett
%  Date:         June 9, 2020
%  Desc:         Main function to fit coefficients a,b,c to simple
%                quadratic model (f(x) = ax^2 + bx + c) for simulated data.
%                Requires Global Optimization Toolbox.

clear all
rng default 

%%%%%%%%%%%%%%%%
%% Generate data
%%%%%%%%%%%%%%%%

nDat = 25; % Number of datapoints
x = sort(unifrnd(0, 1:nDat)); % Randomly chosen input vals
sNoise = 5 * range(x); % Standard deviation for noise of simulated data: 500% of range

% Define true params
a = -9;
b = 6;
c = 4;

% Run model (data to fit quadFun)
modelOut = quadFun([x; x]', [a, b, c]);

% Add noise
noise = normrnd(0, sNoise, nDat);
modelOut = modelOut' + noise(:, 1); 
simData = [x' modelOut]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Format particleswarm options, input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npars = 3; % Number of parameters being optimized

% Set lower, upper bounds for a, b
lb = [-10; -10; -10];
ub = [ 10;  10;  10];

fn = @(params) quadSS(simData, params); % Objective func -> sum of squares

options = optimoptions('particleswarm','MinNeighborsFraction',1); % Global search 
options = optimoptions(options,'PlotFcn',@pswplotbestf); % Plot iterations
options.SwarmSize = 10; % Swarm size of 10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run particleswarm, plot iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[out,fval,exitflag,output] = particleswarm(fn, npars, lb, ub, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot best fit over data
%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
plot(simData(:,1), simData(:,2),'ob');
hold on;
plot(simData(:,1), quadFun(simData, out), '-k');
title('Best fit to observed data');
xlabel('x');
ylabel('y');

%%%%%%%%%%%%%%%%%%
%% Model functions
%%%%%%%%%%%%%%%%%%

% Quadratic function
function model = quadFun(data, params)
x = data(:,1);
for i=1:length(x)
   model(i) = params(1)*x(i)^2 + params(2)*x(i) + params(3);
end
end

% Sum of squares
function SS = quadSS(data, params)
SS = sum((data(:,2) - quadFun(data, params)').^2);
end
