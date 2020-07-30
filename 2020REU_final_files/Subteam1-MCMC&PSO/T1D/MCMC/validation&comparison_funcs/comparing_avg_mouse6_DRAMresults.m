%% Create bar charts to compare means and standard deviations of parameter means
%  Author:       M. Watanabe
%  Date:         July 2020
%  Desc:         Compare means and standard deviations of parameter means
%                using bar charts. Error bars are created from 1 standard 
%                deviation. Cluster results by scale/magnitude of mean.

% DATA: taken from chainstats

% AVERAGE MOUSE RESULTS              | % MOUSE 6 RESULTS
%                  mean         std  |    mean         std      
% ---------------------------------- | ---------------------------------
%         e1 9.9199e-05  5.6117e-05  | 5.9625e-07  4.8352e-07
%         e2 9.2374e-05  5.9001e-05  | 1.0391e-07   5.029e-08
%    delta_B   0.017692    0.001341  | 0.020498   0.0057686
%         SI     0.6345      0.5506  | 1.2582     0.50705   
%         GI     30.125        6.02  | 25.949      7.0361  
%     mues_r  0.0013304  0.00095325  | 3.3264e-06  2.5351e-06
%     mues_e  0.0010473   0.0010099  | 8.8867e-05  5.2077e-05
%  alpha_eta    0.04688    0.029704  | 0.11996    0.057717
%   beta_eta     24.723      2.0772  | 23.507      2.8037  
%  eta_basal   0.012852   0.0020148  | 0.015418   0.0031396
% ---------------------------------- | ----------------------------------
%                                    |
% rejected: 0.433700000000000        |   % rejected: 0.541000000000000


% 10e-06 to 10e-03: e2, e1, mues_r, mues_e - these are off by too large of
% a magnitude so it is very hard to visualize the values together

% 10e-02 to 10e-0.1: delta_B, alpha_eta, eta_basal
x2 = [0.017692,0.020498
      0.04688, 0.11996 
      0.012852, 0.015418];
err2 = [0.001341,0.0057686
        0.029704, 0.057717
        0.0020148, 0.0031396];

       
figure(2)
b2 = bar(x2);
set(gca, 'XTickLabel', {'\delta_B', '\alpha_\eta', '\eta'}, 'FontSize', 50);
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
ylabel('Value')
title('Comparison of Parameter Means')
legend('Average Data', 'Mouse 6')
b2(1).FaceColor = 'k';
b2(2).FaceColor = 'w';
%--------------------------------------------------------------------------
% 10s: GI, beta_eta
x3 = [30.125, 25.949  
      24.723, 23.507];
err3 = [6.02, 7.0361 
        2.0772, 2.8037];
figure(3)
b3 = bar(x3);
set(gca, 'XTickLabel', {'GI', '\beta_{\eta}'});
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
ylabel('Value')
title('Comparison of Parameter Means')
legend('Average Data', 'Mouse 6')
b3(1).FaceColor = 'k';
b3(2).FaceColor = 'w';



