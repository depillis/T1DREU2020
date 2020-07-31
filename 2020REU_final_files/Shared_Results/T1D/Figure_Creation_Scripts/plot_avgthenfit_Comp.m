%% Plot glucose predictions from different parameterizations
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         Function to plot the glucose predictions from PSO, MCMC,
%                UKFs parameterizations. PSO/MCMC run on averaged data.
%                UKFs perform parameterization multiple times and average
%                final parameter values


addpath('averagethenfit_comparison_data')
% load file with all mean predictions
% joint, dual, pso, dram
avg = readmatrix('joint_dual_pso_dram_avgPred.csv');
% load averaged li data
li = readmatrix('avg_lietal_2.csv');
li = sortrows(li);
li = li(1:end-1, :);
li = [li(:,1)*7, li(:,2)];

figure(1); clf
p5 = scatter(li(:,1), li(:,2), 'k','filled', 'o', 'LineWidth', 2); % averaged
hold on
p3 = plot(avg(:,5), avg(:,6), 'LineWidth', 2, 'color', 'g'); % pso
hold on
p4 = plot(avg(:,7), avg(:,8), 'LineWidth', 2, 'color', 'm'); % dram mcmc
hold on
p2 = plot(avg(:,3), avg(:,4), 'LineWidth', 2, 'color', 'r'); % dual ukf
hold on
p1 = plot(avg(:,1), avg(:,2), 'LineWidth', 2, 'color', 'b'); % joint ukf

% labeling
set(gca, 'FontSize', 15)
xlabel('day from birth')
ylabel('glucose (mg/dl)')
legend('Averaged acute', 'PSO', 'MCMC', 'Dual UKF', 'Joint UKF', 'Location', 'northwest');

% CALCULATE MSEs
% joint, dual, pso, dram
mses = readmatrix('joint_dual_pso_dram_avgMSEdata.csv');
mse_joint = error_mse(li(:,2), mses(:, 2)); % joint
mse_dual = error_mse(li(:,2), mses(:, 4)); % dual
mse_pso = error_mse(li(:,2), mses(:, 6)); % pso
mse_dram = error_mse(li(:,2), mses(:, 8)); % dram

final_mse = vertcat(mse_joint, mse_dual, mse_pso, mse_dram);
writematrix(final_mse, 'final_mse_avgComp.csv');




