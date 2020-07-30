function [] = Joint_UKF_Real_Data_Figures(xhat, x, T, truepar, y, param_final)

figure(5)
p1 = plot(1:T, xhat(1,:), 'k', 'Linewidth', 2); hold on;
p2 = plot(1:T, xhat(2,:), 'r', 'Linewidth', 2); hold off;
ylabel('Populations', 'fontsize', 15);
xlabel('Year', 'fontsize', 15);
legend('PreyPredicted', 'PredatorPredicted');
title('State Estimation Using Real D');


%PLOT REAL DATA
figure(6)
p1 = plot(1:T, x(1,:), 'k', 'Linewidth', 2); hold on;
p2 = plot(1:T, x(2,:), 'r', 'Linewidth', 2); hold off;
ylabel('Populations', 'fontsize', 15);
xlabel('Year', 'fontsize', 15);
legend('PreyReal', 'PredatorReal');
title('Raw Predator Prey Data');

%PLOT REAL DATA AND PREDICTIONS ON TOP OF EACH OTHER
figure(9)
subplot(2,1,1);
p1 = plot(1845:1935, x(1,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
p1.MarkerSize = 3;
p2 = plot(1845:1935, xhat(1,:), 'b', 'Linewidth', 2); hold off;
ylabel('population (thousands)', 'fontsize', 15);
xlabel('year', 'fontsize', 15);
legend('Raw', 'Predicted');
title('Prey');

subplot(2,1,2);
p1 = plot(1845:1935, x(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold on;
p1.MarkerSize = 3;
p2 = plot(1845:1935, xhat(2,:), 'r', 'Linewidth', 2); hold off;
ylabel('population (thousands)', 'fontsize', 15);
xlabel('year', 'fontsize', 15);
legend('Raw', 'Predicted');
title('Predator');

%PLOT PARAMETER ESTIMATES
alpha_vec = truepar(1) * ones(1, T);
gamma_vec = truepar(2) * ones(1, T);
beta_vec = truepar(3) * ones(1, T);
delta_vec = truepar(4) * ones(1, T);
figure (11)
subplot(2,2,1);
p1 = plot(1845:1935, alpha_vec, 'k+', 'Linewidth', 2); hold on;
p2 = plot(1845:1935, xhat(3,:), 'b', 'Linewidth', 2); hold off;
ylabel('\alpha', 'fontsize', 15);
xlabel('year', 'fontsize', 15);
legend('Initial', 'Predicted');
title('\alpha Parameter Prediction');

subplot(2,2,2);
p1 = plot(1845:1935, gamma_vec, 'k+', 'Linewidth', 2); hold on;
p2 = plot(1845:1935, xhat(4,:), 'b', 'Linewidth', 2); hold off;
ylabel('\gamma', 'fontsize', 15);
xlabel('year', 'fontsize', 15);
legend('Initial', 'Predicted');
title('\gamma Parameter Prediction');

subplot(2,2,3);
p1 = plot(1845:1935, beta_vec, 'k+', 'Linewidth', 2); hold on;
p2 = plot(1845:1935, xhat(5,:), 'b', 'Linewidth', 2); hold off;
ylabel('\beta', 'fontsize', 15);
xlabel('year', 'fontsize', 15);
legend('Initial', 'Predicted');
title('\beta Parameter Prediction');

subplot(2,2,4);
p1 = plot(1845:1935, delta_vec, 'k+', 'Linewidth', 2); hold on;
p2 = plot(1845:1935, xhat(6,:), 'b', 'Linewidth', 2); hold off;
%ylim([0 0.1]);
ylabel('\delta', 'fontsize', 15);
xlabel('year', 'fontsize', 15);
legend('Initial', 'Predicted');
title('\delta Parameter Prediction');

%PLOT THE ERROR
x_error = x(1:2,:) - xhat(1:2,:); %use matrix subtraction to get error
error_norm = vecnorm(x_error); %take columnwise norm
figure(20) %Plot norm of the error over time
p1 = plot(1:T, error_norm(1,:), 'b', 'Linewidth', 2);
ylabel('Error Norm', 'fontsize', 15);
xlabel('Time', 'fontsize', 15);
title('Norm of Error Over Time');

%PLOT ODE WITH FINAL PARAMETERS OVER DATA
tspan = [0 91];
x0 = [y(1,1); y(2,1); 0.62526; 0.6607; 0.1896; 0.0468];
sol = ode45(@(t, y) Lotka_Volterra_Model(t, y, param_final), tspan, x0); %Use ODE solver
x_with_final = zeros(6, 91);
for j = 1:91
    x_with_final(:,j) = deval(sol, j);
end

figure(100)
subplot(2,1,1);
p1 = plot(1845:1935, x_with_final(1,:), 'b', 'Linewidth', 2); hold on;
p2 = plot(1845:1935, x(1,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold off;
p2.MarkerSize = 3;
ylabel('population (thousands)', 'fontsize',15);
xlabel('year', 'fontsize', 15);
legend('ODE Fit', 'Raw');
title('Prey');

subplot(2,1,2);
p1 = plot(1845:1935, x_with_final(2,:), 'r', 'Linewidth', 2); hold on;
p2 = plot(1845:1935, x(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold off;
p2.MarkerSize = 3;
ylabel('population (thousands)', 'fontsize',15);
xlabel('year', 'fontsize', 15);
legend('ODE Fit', 'Raw');
title('Predator');

tspan = [0 91];
params0 = [0.62526 0.6607 0.1896 0.0468];
x0 = [y(1,1); y(2,1); 0.62526; 0.6607; 0.1896; 0.0468];
sol = ode45(@(t, y) Lotka_Volterra_Model(t, y, params0), tspan, x0); %Use ODE solver
x_start = zeros(6, 91);
for j = 1:91
    x_start(:,j) = deval(sol, j);
end

figure(102)
p1 = plot(1:T, x_start(1,:), 'b', 'Linewidth', 2); hold on;
p2 = plot(1:T, x(1,:), 'b+', 'Linewidth', 2); hold off;
ylabel('Prey Population', 'fontsize', 15);
xlabel('Year', 'fontsize', 15);
title('Prey Simulated Using Initial Paramaters');

figure(103)
p1 = plot(1:T, x_start(2,:), 'b', 'Linewidth', 2); hold on;
p2 = plot(1:T, x(2,:), 'b+', 'Linewidth', 2); hold off;
ylabel('Predator Population', 'fontsize', 15);
xlabel('Year', 'fontsize', 15);
title('Predator Simulated Using Initial Paramaters');

%CALCULATE MODEL ERROR WITH FINAL PARAMETERS
prey_final_error = x(1,:) - x_with_final(1,:);
prey_final_error = prey_final_error.^2;
prey_final_error = sum(prey_final_error);
MSE_prey_final = prey_final_error /91;
MSE_prey_final

predator_final_error = x(2,:) - x_with_final(2,:);
predator_final_error = predator_final_error.^2;
predator_final_error = sum(predator_final_error);
MSE_predator_final = predator_final_error /91;
MSE_predator_final

%CALCULATE MODEL ERROR WITH STARTING PARAMETERS
prey_start_error = x(1,:) - x_start(1,:);
prey_start_error = prey_start_error.^2;
prey_start_error = sum(prey_start_error);
MSE_prey_start = prey_start_error /91;
MSE_prey_start

predator_start_error = x(2,:) - x_start(2,:);
predator_start_error = predator_start_error.^2;
predator_start_error = sum(predator_start_error);
MSE_predator_start = predator_start_error / 91;
MSE_predator_start



