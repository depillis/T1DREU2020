%Script to produce figure that shows performance of Joint UKF, Dual UKF,
%and PSO on fit-then-average approach for averaged data

clear all;

li = load('avglietal.mat');
li_mat = li.avglietal2;
X = li_mat{:,:};

X = X'; %set up table as wide
X = X(:, 1:11);

fit_first = 'Fit_Then_Average_Joint.csv';
Joint_fit_sol = readtable(fit_first);
Joint_fit_sol = Joint_fit_sol{:,:};

fit_first = 'Fit_Then_Average_Dual.csv';
Dual_fit_sol = readtable(fit_first);
Dual_fit_sol = Dual_fit_sol{:,:};

fit_first = 'Fit_Then_Average_PSO.mat';
PSO_fit_sol = load(fit_first);
PSO_fit_sol = PSO_fit_sol.fit_then_avg;

fig = figure(1);
p1 = plot(7 * X(1,:), X(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
p1.MarkerSize = 3;
p2 = plot(Joint_fit_sol(:,1), Joint_fit_sol(:,2), 'b', 'Linewidth', 2);
p3 = plot(Dual_fit_sol(:,1), Dual_fit_sol(:,2), 'r', 'Linewidth', 2);
p4 = plot(PSO_fit_sol(:,1), PSO_fit_sol(:,2), 'g', 'Linewidth', 2);
ylabel('glucose(mg/dl)', 'Fontsize', 15);
xlabel('days', 'Fontsize', 15);
legend('raw', 'Joint UKF', 'Dual UKF', 'PSO');
