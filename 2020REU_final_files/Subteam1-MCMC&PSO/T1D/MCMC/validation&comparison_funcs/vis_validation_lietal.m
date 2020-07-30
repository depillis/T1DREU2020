%% Visual Validation (using Li et al data)
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         Compare predicted model (using DRAM parameters) with the
%                Li et al data by plotting both.

function vis_validation_lietal(workspace, modeltype, disease, parVer, wave)
addpath('../T1D_data_files/Lietal_data');
addpath('../mcmcstat');
addpath('../past_run_workspaces');

% LOAD DATA
% load previous workspace you want to use to compare fits
load(workspace);

% load Li data & shift time to days not weeks
dat1 = readmatrix('dat1.csv');
dat1 =[dat1(:,1)*7 dat1(:,2)];
dat1 = sortrows(dat1);

dat2 = readmatrix('dat2.csv');
dat2 =[dat2(:,1)*7 dat2(:,2)];
dat2 = sortrows(dat2);

dat3 = readmatrix('dat3.csv');
dat3 =[dat3(:,1)*7 dat3(:,2)];
dat3 = sortrows(dat3);

dat4 = readmatrix('dat4.csv');
dat4 =[dat4(:,1)*7 dat4(:,2)];
dat4 = sortrows(dat4);

dat5 = readmatrix('dat5.csv');
dat5 =[dat5(:,1)*7 dat5(:,2)];
dat5 = sortrows(dat5);

dat6 = readmatrix('dat6.csv');
dat6 =[dat6(:,1)*7 dat6(:,2)];
dat6 = sortrows(dat6);

dat7 = readmatrix('dat7.csv');
dat7 =[dat7(:,1)*7 dat7(:,2)];
dat7 = sortrows(dat7);

dat8 = readmatrix('dat8.csv');
dat8 =[dat8(:,1)*7 dat8(:,2)];
dat8 = sortrows(dat8);

dat9 = readmatrix('dat9.csv');
dat9 =[dat9(:,1)*7 dat9(:,2)];
dat9 = sortrows(dat9);

dat10 = readmatrix('dat10.csv');
dat10 =[dat10(:,1)*7 dat10(:,2)];
dat10 = sortrows(dat10);

dat11 = readmatrix('dat11.csv');
dat11 =[dat11(:,1)*7 dat11(:,2)];
dat11 = sortrows(dat11);


% PLOT MEAN PARAM VALUE PREDICTION
meanParams = results.theta;
[t, y] = T1Dfun_meanpred(modeltype, parVer, wave, meanParams, data);
[t2, y2] = T1Dfun_r(modeltype, parVer, wave, meanParams, data);

figure(1); clf
% plot predicted model
plot(t,y(:,6), 'LineWidth',2);
hold on
scatter(data(:,1), data(:,2), 'LineWidth', 2);
hold off
%plot Li et al data
plot(dat1(:,1), dat1(:,2), '--');
hold on
plot(dat2(:,1), dat2(:,2), '-*');
hold on
plot(dat3(:,1), dat3(:,2), '-*');
hold on
plot(dat4(:,1), dat4(:,2), '-*');
hold on
plot(dat5(:,1), dat5(:,2), '--');
hold on
plot(dat6(:,1), dat6(:,2), '-*');
hold on
plot(dat7(:,1), dat7(:,2), '-*');
hold on
plot(dat8(:,1), dat8(:,2), '-*');
hold on
plot(dat9(:,1), dat9(:,2), '-*');
hold on
plot(dat10(:,1), dat10(:,2), '-*');
hold on
plot(dat11(:,1), dat11(:,2), '-*');
hold on

% Axes & legend
set(gca, 'FontSize', 20);
xlabel('days from birth');
ylabel('glucose (mg/dl)');

if disease == 0
    legend('Model prediction', 'Li data, progressive', 'Location', 'southeast');
    %title('Glucose model prediction w/mean parameter values: progressive','FontSize', 15);

else
    legend('Model prediction', 'Li data, acute', 'Location', 'southeast');
    %title('Glucose model prediction w/mean parameter values: acute','FontSize', 15);
    
end
% figure(2); clf
% % plot truncated predicted model
% plot(t2,y2, 'LineWidth',2);
% hold on
% scatter(data(:,1), data(:,2), 'LineWidth',2);
% hold off
% 
% % plot Li et al data
% % plot(dat1(:,1), dat1(:,2), '-*');
% hold on
% plot(dat2(:,1), dat2(:,2), '-*');
% hold on
% plot(dat3(:,1), dat3(:,2), '-*');
% hold on
% plot(dat4(:,1), dat4(:,2), '-*');
% hold on
% % plot(dat5(:,1), dat5(:,2), '-*');
% hold on
% plot(dat6(:,1), dat6(:,2), '-*');
% hold on
% plot(dat7(:,1), dat7(:,2), '-*');
% hold on
% plot(dat8(:,1), dat8(:,2), '-*');
% hold on
% plot(dat9(:,1), dat9(:,2), '-*');
% hold on
% plot(dat10(:,1), dat10(:,2), '-*');
% hold on
% plot(dat11(:,1), dat11(:,2), '-*');
% 
% % Axes & legend
% set(gca, 'FontSize', 15)
% xlabel('days');
% ylabel('glucose (mg/dL)');
% 
% if disease == 0
%     legend('Model prediction', 'Li data, progressive', 'Location', 'southeast');
%     title('Glucose model prediction w/mean parameter values: progressive','FontSize', 15);
% 
% else
%     legend('Model prediction', 'Li data, acute', 'Location', 'southeast');
%     title('Glucose model prediction w/mean parameter values: acute','FontSize', 15);
%     
% end

% % PLOT MCMCSTAT MODPRED
% figure(2); clf
% plot(dat1(:,1), dat1(:,2), '-*');
% hold on
% plot(dat2(:,1), dat2(:,2), '-*');
% hold on
% plot(dat3(:,1), dat3(:,2), '-*');
% hold on
% plot(dat4(:,1), dat4(:,2), '-*');
% hold on
% plot(dat5(:,1), dat5(:,2), '-*');
% hold on
% plot(dat6(:,1), dat6(:,2), '-*');
% hold on
% plot(dat7(:,1), dat7(:,2), '-*');
% hold on
% plot(dat8(:,1), dat8(:,2), '-*');
% hold on
% plot(dat9(:,1), dat9(:,2), '-*');
% hold on
% plot(dat10(:,1), dat10(:,2), '-*');
% hold on
% plot(dat11(:,1), dat11(:,2), '-*');
% 
% hold on
% 
% % plot mcmcstat predplot
% mcmcpredplot(out); % altered mcmcpredplot to create a 4x3 pane figure
% hold on
% scatter(data(:,1), data(:,2)); % our data is only from Glucose measurements so we only plot data over 1 of the densities
% 
% xlabel('days');
% ylabel('glucose (mg/dL)');
% 
% if disease == 0
%     title('Glucose model prediction: progressive');
% 
% else
%     title('Glucose model prediction: acute');
% end
end
