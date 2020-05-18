%This function loads and plots the data from the algaedata example
%and related ODE model.

%Clear workspace
clear model data params options

%Load the data structure
load algaedata.mat

%Plot the data first distributions of xdata
figure(1); clf
for i =1:3
  subplot(2,3,i)
  plot(data.xdata(:,1),data.xdata(:,i+1),'-k');
  title(data.xlabels(i+1)); xlim([1,120])
end

%Plot the state variable time courses
subplot(2,1,2)
plot(data.ydata(:,1),data.ydata(:,2:end),'o-');
title('model state variable observations');
legend(data.ylabels(2:end),'Location','best');
xlabel('days');