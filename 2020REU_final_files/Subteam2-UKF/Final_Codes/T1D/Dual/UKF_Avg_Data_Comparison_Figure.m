%Script to produce figure that compares averaging techniques across Dual
%UKF, Joint UKF, and PSO

clear all;

li = load('avglietal.mat');
li_mat = li.avglietal2;
X = li_mat{:,:};

X = X'; %set up table as wide
X = X(:, 1:11);

average_first = 'Average_Then_Fit_Joint.csv';
Joint_avg_sol = readtable(average_first);
Joint_avg_sol = Joint_avg_sol{:,:};

fit_first = 'Fit_Then_Average_Joint.csv';
Joint_fit_sol = readtable(fit_first);
Joint_fit_sol = Joint_fit_sol{:,:};

average_first = 'Average_Then_Fit_Dual.csv';
Dual_avg_sol = readtable(average_first);
Dual_avg_sol = Dual_avg_sol{:,:};

fit_first = 'Fit_Then_Average_Dual.csv';
Dual_fit_sol = readtable(fit_first);
Dual_fit_sol = Dual_fit_sol{:,:};

fit_first = 'Fit_Then_Average_PSO.mat';
PSO_fit_sol = load(fit_first);
PSO_fit_sol = PSO_fit_sol.fit_then_avg;

average_first = 'Average_Then_Fit_PSO.csv';
PSO_avg_sol = readtable(average_first);
PSO_avg_sol = PSO_avg_sol{:,:};

fig = figure(1);
subplot(2,2,1);
p1 = plot(7 * X(1,:), X(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
p1.MarkerSize = 3;
p5 = plot(Joint_fit_sol(:,1), Joint_fit_sol(:,2), 'b', 'Linewidth', 2);
p2 = plot(Joint_avg_sol(:,1), Joint_avg_sol(:,2), 'b--', 'Linewidth', 2); hold off;
legend('Averaged Acute', 'Fit First', 'Avg First', 'Location', 'northwest');
title('Joint UKF');


subplot(2,2,2);
p1 = plot(7 * X(1,:), X(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
p1.MarkerSize = 3;
p5 = plot(Dual_fit_sol(:,1), Dual_fit_sol(:,2), 'r', 'Linewidth', 2);
p2 = plot(Dual_avg_sol(:,1), Dual_avg_sol(:,2), 'r--', 'Linewidth', 2); hold off;
legend('Averaged Acute', 'Fit First', 'Avg First', 'Location', 'northwest');
title('Dual UKF');

subplot(2,2,3);
p1 = plot(7 * X(1,:), X(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
p1.MarkerSize = 3;
p5 = plot(PSO_fit_sol(:,1), PSO_fit_sol(:,2), 'g', 'Linewidth', 2);
p2 = plot(PSO_avg_sol(:,1), PSO_avg_sol(:,2), 'g--', 'Linewidth', 2); hold off;
legend('Averaged Acute', 'Fit First', 'Avg First', 'Location', 'northwest');
title('PSO');


subplot(2,2,4);
p1 = plot(7 * X(1,:), X(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
p1.MarkerSize = 3;
p4 = plot(PSO_fit_sol(:,1), PSO_fit_sol(:,2), 'g', 'Linewidth', 2);
p3 = plot(Dual_fit_sol(:,1), Dual_fit_sol(:,2), 'r', 'Linewidth', 2);
p2 = plot(Joint_fit_sol(:,1), Joint_fit_sol(:,2), 'b', 'Linewidth', 2);
legend('Averaged Acute', 'Fit First PSO', 'Fit First Dual UKF', 'Fit First Joint UKF', 'Location', 'northwest');





han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han, 'glucose (mg/dl)', 'fontsize', 15);
xlabel(han, 'days since birth', 'fontsize', 15);
yh = get(gca,'ylabel'); % handle to the label object
p = get(yh,'position'); % get the current position property
p(1) = 1.5*p(1);
set(yh,'position',p); 

figure(5)
p1 = plot(7 * X(1,:), X(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
p1.MarkerSize = 3;
p4 = plot(PSO_fit_sol(:,1), PSO_fit_sol(:,2), 'g', 'Linewidth', 2);
p3 = plot(Dual_fit_sol(:,1), Dual_fit_sol(:,2), 'r', 'Linewidth', 2);
p2 = plot(Joint_fit_sol(:,1), Joint_fit_sol(:,2), 'b', 'Linewidth', 2);
legend('Averaged Acute', 'Fit First PSO', 'Fit First Dual UKF', 'Fit First Joint UKF', 'Location', 'northwest');
xlabel('days since birth', 'fontsize', 15);
ylabel('glucose (mg/dl)', 'fontsize', 15);

