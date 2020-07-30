%% Compare MH and DRAM parameter means
% Authors:      M. Watanabe
%
% Date:         July 2020
%
% Descr:        Script to compare the means of the Lotka-Volterra
%               parameters that result from the MH and DRAM
%               parameterizations. This script generates bar charts of the
%               means of the 4 parameters with error bars created from 1
%               standard deviation.
                
%               Note that this code uses mean and standard deviations
%               specifically from the workspaces "final_mh.mat" and
%               "final_dram.mat". To plot means from different runs, use
%               must specify values.
%
% Directions:

% plotting mean parameter values with standard deviation bars
% keep values at same magnitude together


% MH                                  DRAM
%                 mean         std          mean         std
% -----------------------------------|----------------------------------
%      alpha    0.77857     0.13427  |     0.78674     0.13898
%       beta   0.030054   0.0054025  |     0.030346   0.0055804
%      gamma    0.41156    0.087064  |     0.40665    0.090777
%      delta   0.009945   0.0023043  |     0.0098302   0.0023888
% -----------------------------------|----------------------------------
% 

x2 = [0.77857, 0.78674
      0.41156, 0.40665];
err2 = [0.13427, 0.13898
        0.087064, 0.090777];
       
figure(2)
b2 = bar(x2);
set(gca, 'XTickLabel', {'\alpha', '\gamma'}, 'FontSize', 25);
hold on
ngroups = size(x2, 1);
nbars = size(x2, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, x2(:,i), err2(:,i), '.');
end

% Labelling
set(gca, 'FontSize', 15)
xlabel('Parameter')
ylabel('parameter value')
%title('Comparison of Parameter Means')
legend('Metropolis', 'DRAM')
b2(1).FaceColor = 'k';
b2(2).FaceColor = 'w';
%--------------------------------------------------------------------------

x3 = [0.030054, 0.030346   
      0.009945, 0.0098302];
err3 = [0.0054025, 0.0055804
        0.0023043, 0.0023888];
figure(3)
b3 = bar(x3);
set(gca, 'XTickLabel', {'\beta', '\delta'});
hold on
ngroups = size(x3, 1);
nbars = size(x3, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, x3(:,i), err3(:,i), '.');
end
hold off
set(gca, 'FontSize', 15)
xlabel('Parameter')
ylabel('parameter value')
%title('Comparison of Parameter Means')
legend('MH', 'DRAM')
b3(1).FaceColor = 'k';
b3(2).FaceColor = 'w';



