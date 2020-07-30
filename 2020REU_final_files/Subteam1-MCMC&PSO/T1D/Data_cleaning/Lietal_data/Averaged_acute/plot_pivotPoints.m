% Script to plot pivot points of 
% load data
addpath('../T1D_data_files/Lietal_data')
dat2 = readmatrix('dat2.csv');
dat2 =[dat2(:,1)*7 dat2(:,2)];
dat2 = sortrows(dat2);

dat3 = readmatrix('dat3.csv');
dat3 =[dat3(:,1)*7 dat3(:,2)];
dat3 = sortrows(dat3);

dat4 = readmatrix('dat4.csv');
dat4 =[dat4(:,1)*7 dat4(:,2)];
dat4 = sortrows(dat4);

% dat5 = readmatrix('dat5.csv');
% dat5 =[dat5(:,1)*7 dat5(:,2)];
% dat5 = sortrows(dat5);

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

% plot

plot(dat2(:,1), dat2(:,2), '--o', 'color', 'k', 'LineWidth', 2);
hold on
scatter(231.091, 206.3, 'r', 'o', 'filled')
hold on
plot(dat3(:,1), dat3(:,2), '--o','color', 'k', 'LineWidth', 2);
hold on
scatter(189.471, 220.4, 'r', 'o', 'filled')
hold on
plot(dat4(:,1), dat4(:,2), '--o','color', 'k', 'LineWidth', 2);
hold on
scatter(189.357, 234.51, 'r', 'o', 'filled')
hold on
plot(dat6(:,1), dat6(:,2), '--o','color', 'k', 'LineWidth', 2);
hold on
scatter(175.161, 181.61, 'r', 'o', 'filled')
hold on
plot(dat7(:,1), dat7(:,2), '--o','color', 'k', 'LineWidth', 2);
hold on
scatter(175.126, 222.17, 'r', 'o', 'filled')
hold on
plot(dat8(:,1), dat8(:,2), '--o','color', 'k', 'LineWidth', 2);
hold on
scatter(160.867, 234.51, 'r', 'o', 'filled')
hold on
plot(dat9(:,1), dat9(:,2), '--o','color', 'k', 'LineWidth', 2);
hold on
scatter(160.839, 264.48, 'r', 'o', 'filled')
hold on
plot(dat10(:,1), dat10(:,2), '--o','color', 'k', 'LineWidth', 2);
hold on
scatter(147.28, 259.19, 'r', 'o', 'filled')
hold on
plot(dat11(:,1), dat11(:,2), '--o','color', 'k', 'LineWidth', 2);
hold on
scatter(132.748, 208.06, 'r', 'o', 'filled')
hold on
set(gca, 'FontSize', 15)
xlabel('days from birth')
ylabel('glucose (mg/dl)')
legend('Acute mouse data', 'Pivot Points', 'Location', 'northwest')