function [] = T1D_DualUKF_Figures(numReadings,X_final, xhat, Nruns, params_alltime, all_param_estimates, X, wave)

fig_2_3_4_5_Parameters;
T = numReadings;

 %-----PLOT GLUCOSE ESTIMATES--------%
    figure(2)
    p1 = plot(1:T, X_final(2,:), 'r+', 'Linewidth', 2); hold on;
    p2 = plot(1:T, xhat(6,:), 'r', 'Linewidth', 2); hold off;
    ylabel('Glucose', 'fontsize', 15);
    xlabel('Week', 'fontsize', 15);
    legend('G generated', 'G predicted');
    
    %------PLOT PARAMETER ESTIMATES------%
    figure (3)
    p2 = plot(1: T * Nruns, params_alltime(24,:), 'b', 'Linewidth', 2); 
    ylabel('D_ss', 'fontsize', 15);
    xlabel('Day', 'fontsize', 15);


    figure (4)
    p2 = plot(1: T * Nruns, params_alltime(9,:), 'b', 'Linewidth', 2); 
    ylabel('alpha_B', 'fontsize', 15);
    xlabel('Day', 'fontsize', 15);
  

    figure (5)
    p2 = plot(1: T * Nruns, params_alltime(18,:), 'b', 'Linewidth', 2); 
    ylabel('GI', 'fontsize', 15);
    xlabel('Day', 'fontsize', 15);

    figure (6)
    p2 = plot(1: T * Nruns, params_alltime(15,:), 'b', 'Linewidth', 2);
    ylabel('SI', 'fontsize', 15);
    xlabel('Day', 'fontsize', 15);

    figure (7)
    p2 = plot(1: T * Nruns, params_alltime(19,:), 'b', 'Linewidth', 2); 
    ylabel('Qpanc', 'fontsize', 15);
    xlabel('Day', 'fontsize', 15);

    figure (8)
    p2 = plot(1: T * Nruns, params_alltime(35,:), 'b', 'Linewidth', 2);
    ylabel('mues_e', 'fontsize', 15);
    xlabel('Day', 'fontsize', 15);


    figure (9)
    p2 = plot(1: T * Nruns, params_alltime(34,:), 'b', 'Linewidth', 2); 
    ylabel('mues_r', 'fontsize', 15);
    xlabel('Day', 'fontsize', 15);

    
X(1,:) = X(1,:) * 7;
    
T = 350;
final_params = all_param_estimates(:,5);
%final_params = all_param_estimates(:,3);
%final_params = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8); 
 %                    f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
  %                   1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
   %                  0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21];
%final_params(1) = 10^5 - 300
%final_params = [10^5 - 300; 0.04; sqrt(18000); .8; .210; aEmday/5000; aEmday/5000];
%final_params(1:6) = all_param_estimates(1:6,1);
%final_params = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8); 
%                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
%                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
%                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21];
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
final_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, final_params), [0 350], InitialState,options); %Simulate ODE's up until first time point
x_fake = zeros(12, 350);
for j = 1:350
    x_fake(:,j) = deval(final_sol, j);
end
figure(205)
p1 = plot(1:T, x_fake(6,:), 'r', 'Linewidth', 2); hold on;
p2 = plot(X(1,:), X(2,:), 'r+', 'Linewidth', 2); hold off;
ylabel('Glucose Set 5', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);
x_predicted = deval(final_sol, X(1,:));
error = X(2,:) - x_predicted(6,:);
error_norm = norm(error);
error_norm


%PLOT THE STATE ESTIMATES WITH WAVE

figure(301)
sgtitle('Parameter Estimates + Wave');
subplot(3,4,1);
plot(1:T, x_fake(1,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Resting Macrophages', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,2);
plot(1:T, x_fake(2,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Active Mcrophages', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,3);
plot(1:T, x_fake(3,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Healthy Beta Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,4);
plot(1:T, x_fake(4,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Apoptotic Beta Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,5);
plot(1:T, x_fake(5,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Necrotic Beta Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,6);
plot(1:T, x_fake(7,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Insulin', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,7);
plot(1:T, x_fake(8,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('D Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,8);
plot(1:T, x_fake(9,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('tD Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,9);
plot(1:T, x_fake(10,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Effector T Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,10);
plot(1:T, x_fake(11,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Regulatory T Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,11);
plot(1:T, x_fake(12,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Memory T Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);



%PLOT STATE ESTIMATES IF NO WAVE
final_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, 0, final_params), [0 350], InitialState,options); %Simulate ODE's up until first time point
x_fake_no_wave = zeros(12, 350);
for j = 1:350
    x_fake_no_wave(:,j) = deval(final_sol, j);
end


figure(302)
sgtitle('Parameter Estimates + NO Wave');
subplot(3,4,1);
plot(1:T, x_fake_no_wave(1,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Resting Macrophages', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,2);
plot(1:T, x_fake_no_wave(2,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Active Mcrophages', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,3);
plot(1:T, x_fake_no_wave(3,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Healthy Beta Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,4);
plot(1:T, x_fake_no_wave(4,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Apoptotic Beta Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,5);
plot(1:T, x_fake_no_wave(5,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Necrotic Beta Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,6);
plot(1:T, x_fake_no_wave(7,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Insulin', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,7);
plot(1:T, x_fake_no_wave(8,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('D Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,8);
plot(1:T, x_fake_no_wave(9,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('tD Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,9);
plot(1:T, x_fake_no_wave(10,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Effector T Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,10);
plot(1:T, x_fake_no_wave(11,:), 'r', 'Linewidth', 2);
set(gca, 'YScale', 'log')
ylabel('Regulatory T Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,11);
plot(1:T, x_fake_no_wave(12,:), 'r', 'Linewidth', 2); 
set(gca, 'YScale', 'log')
ylabel('Memory T Cells', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

subplot(3,4,12);
plot(1:T, x_fake_no_wave(6,:), 'r', 'Linewidth', 2);
ylabel('Glucose', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);



T = 280;
final_params = all_param_estimates(:,2);
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
final_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, final_params), [0 280], InitialState,options); %Simulate ODE's up until first time point
x_fake = zeros(12, 280);
for j = 1:280
    x_fake(:,j) = deval(final_sol, j);
end
figure(202)
p1 = plot(1:T, x_fake(6,:), 'r', 'Linewidth', 2); hold on;
p2 = plot(X(1,:), X(2,:), 'r+', 'Linewidth', 2); hold off;
ylabel('Glucose Set 2', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);
x_predicted = deval(final_sol, X(1,:));
error = X(2,:) - x_predicted(6,:);
error_norm = norm(error);
error_norm



T = 280;
final_params = all_param_estimates(:,1);
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
final_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, final_params), [0 280], InitialState,options); %Simulate ODE's up until first time point
x_fake = zeros(12, 280);
for j = 1:280
    x_fake(:,j) = deval(final_sol, j);
end
figure(201)
p1 = plot(1:T, x_fake(6,:), 'r', 'Linewidth', 2); hold on;
p2 = plot(X(1,:), X(2,:), 'r+', 'Linewidth', 2); hold off;
ylabel('Glucose Set 1', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);
x_predicted = deval(final_sol, X(1,:));
error = X(2,:) - x_predicted(6,:);
error_norm = norm(error);
error_norm





end

