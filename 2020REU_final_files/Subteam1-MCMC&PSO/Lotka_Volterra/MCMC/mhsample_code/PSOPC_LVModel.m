clear all
close all

% Select years 1908-1935
idx = 60:91;
mytime = (1:length(idx))';
load('HaresLynxData.mat')
mydata(:,1) = Lotka_Volterra_Data(idx,2);
mydata(:,2) = Lotka_Volterra_Data(idx,3);


objfun = @(x) least_squares(x,mydata, mytime);

%% Optimisation using PSOPC with the same upper/lower bounds used in
%% simplex method
[k least_squares] = PSOPC(objfun, 4, [0.4 0 0.4 0], [0.8 0.8 0.8 0.8], 200);

%% Plot model with estimated parameters
%y0(1) = 21.5; y0(2) = 3.4;
y0(1) = mydata(1,1); y0(2) = mydata(1,2);
[t,y] = ode45(@Lotka_Volterra_Model,mytime,y0,[],k);

subplot(2,1,1)
hold on
plot(mydata(:,1),'O');
plot(y(:,1),'--b')

subplot(2,1,2)
hold on
plot(mydata(:,2),'rO');
plot(y(:,2),'--r')

