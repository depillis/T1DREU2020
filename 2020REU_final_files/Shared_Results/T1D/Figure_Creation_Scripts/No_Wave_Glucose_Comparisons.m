clear all;
file = 'all_values_IMPORTANT.csv';
vals = readtable(file);
vals = vals{:,:};
T = 350;
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
fig_2_3_4_5_Parameters;
wave = wave_basal;
m = 6;
m_str = int2str(m);
mouse_dat_file = strcat('dat', m_str,'.csv');
mouse_dat = readtable(mouse_dat_file);
mouse_dat = mouse_dat{:,:};
mouse_dat = flip(mouse_dat);
mouse_dat = mouse_dat';
mouse_dat(1,:) = mouse_dat(1,:) * 7;
num_glucose = size(mouse_dat, 2);
mouse_params = vals(:,m);
sim_params = mouse_params(2:42);
options=odeset('RelTol',1e-12,'AbsTol',1e-12);

no_wave_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, 0, sim_params), [0 350], InitialState,options); %Simulate ODE's up until first time point
x_sim_no_wave = zeros(12, 350);
for j = 1:350
    x_sim_no_wave(:,j) = deval(no_wave_sol, j);
end

plot(1:T, x_sim_no_wave(6,:), 'r', 'Linewidth', 2); hold on;
ylabel('Glucose', 'fontsize', 15);
xlabel('Day', 'fontsize', 15);

file = 'mcmc_nowave_glucose.csv';
mcmc_glucose = readtable(file);
mcmc_glucose = mcmc_glucose{:,:};
mcmc_glucose = mcmc_glucose';
plot(0:T, mcmc_glucose(2,:), 'm', 'Linewidth', 2);

file = 'mouse6_joint_nowave.csv';
joint_glucose = readtable(file);
joint_glucose = joint_glucose{:,:};
joint_glucose = joint_glucose';
plot(joint_glucose(1,:), joint_glucose(2,:), 'b', 'Linewidth', 2);

file = 'PSO_m6_nowave.mat';
PSO_sol = load(file);
pso_glucose = PSO_sol.PSOnowave;
pso_glucose = pso_glucose';
plot(pso_glucose(1,:), pso_glucose(2,:), 'g', 'Linewidth', 2);

baseline_params = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8);
                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21; 0.01];
                 
base_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, 0, baseline_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
x_base_nowave = zeros(12, 350);
for j = 1:350
    x_base_nowave(:,j) = deval(base_sol, j);
end


plot(1:T, x_base_nowave(6,:), 'k--', 'Linewidth', 2); hold off;
xlabel('days since birth', 'fontsize', 15);
ylabel('glucose (mg/dl)', 'fontsize', 15);
legend('Dual UKF', 'DRAM MCMC', 'Joint UKF', 'PSO', 'Baseline');





