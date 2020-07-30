%CREATES THE FIGURE OF ALL THE VARIOUS ALGORITHMS FITS ON LI MOUSE 6 DATA

%PLOT ALL FITS ON TOP OF EACH OTHER
file = 'dat6.csv'; %create file name
mouse6 = readtable(file);
mouse6 = mouse6{:,:};
mouse6 = flip(mouse6);
mouse6 = mouse6';

file = 'mouse6PSOfit.csv';
PSO_sol = readtable(file);
PSO_sol = PSO_sol{:,:};

file = 'jul10_mouse_run1(noIC)_acute_NOD_waveOn_lietal_predmodData.csv';
MCMC_sol = readtable(file);
MCMC_sol = MCMC_sol{:,:};

file = 'mouse6_DualUKFSol.csv';
Dual_sol = readtable(file);
Dual_sol = Dual_sol{:,:};

file = 'mouse6_joint.csv';
Joint_sol = readtable(file);
Joint_sol = Joint_sol{:,:};

figure(1);
p1 = plot(7 * mouse6(1,:), mouse6(2,:), 'k+', 'Linewidth', 2); hold on;
p2 = plot(PSO_sol(:,1), PSO_sol(:,2), 'g', 'Linewidth', 2);
p3 = plot(MCMC_sol(:,1), MCMC_sol(:,2), 'm', 'Linewidth', 2);
p4 = plot(Dual_sol(:,1), Dual_sol(:,2), 'r', 'Linewidth', 2);
p5 = plot(Joint_sol(:,1), Joint_sol(:,2), 'b', 'Linewidth', 2);
xlabel('Day', 'fontsize', 15);
ylabel('Glucose', 'fontsize', 15);
legend('Raw', 'PSO', 'MCMC', 'Dual UKF', 'Joint UKF');
ylim([0 600]);

%CALCULATE ALL MSE VALUES
file = 'jul10_mouse_run1(noIC)_acute_NOD_waveOn_lietal_errorestData.csv';
MCMC_mse = readtable(file);
MCMC_mse = MCMC_mse{:,:};
MCMC_mse = MCMC_mse';
MCMC_final = (mouse6(2,:) - MCMC_mse(2,:)).^2;
MCMC_final  = sum(MCMC_final);
MCMC_final = MCMC_final / 10;
MCMC_final

file = 'mouse6PSOmse.csv';
PSO_mse = readtable(file);
PSO_mse = PSO_mse{:,:};
PSO_mse = PSO_mse';
PSO_final = (mouse6(2,:) - PSO_mse(2,:)).^2;
PSO_final  = sum(PSO_final);
PSO_final = PSO_final / 10;
PSO_final



