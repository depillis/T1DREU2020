load('activedata.mat')
load('mouse6PSOfit.csv')
load('liMouse6.mat');

% Select time series, data of interest
t = dat6.VarName1.*7;
y = dat6.VarName2;

plot(t,y,'+k', 'LineWidth', 2)
hold on;
plot(avgdataPSOfit(:,1),avgdataPSOfit(:,2),'-g', 'LineWidth', 2)
plot(mouse6PSOfit(:,1),mouse6PSOfit(:,2), '-m', 'LineWidth', 2)
legend('Observed','Active parameter subset', 'All parameters', 'Location', 'best')
set(gca,'FontSize',18);
ylabel('glucose (ug/dl)');
xlabel('day from birth');
