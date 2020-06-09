% Solve the single compartment model and plot the results


function fig_2_3_4_5_Solver
close all
clear all

%Get parameter values
fig_2_3_4_5_Parameters;

%Initial Values for DCs and tDCs
D0 = 0;
tD0 = 0;

%Initial values
EndTime =1000;% 300;
Tspan = [0 EndTime]; % time in days
IC = [4.77*10^5 0 0 0 300 100 10 D0 tD0 0 0 0];% start at the healthy rest state, using Topp healthy rest state for beta cells, glucose, insulin


%Wave is off
wave = 0;

%Solve ODE for wild type
[Tnowave, Ynowave] = ode15s(@(t,y)fig_2_3_4_5_ODE(t,y,f1,f2,wave),Tspan,IC); % Solve ODE

%Solve ODE for NOD
[Tnnowave, Ynnowave] = ode15s(@(t,y)fig_2_3_4_5_ODE(t,y,f1n,f2n,wave),Tspan,IC); % Solve ODE


%Wave is on
wave = wave_basal;%0.75;

%Solve ODE for wild type
[Twave, Ywave] = ode15s(@(t,y)fig_2_3_4_5_ODE(t,y,f1,f2,wave),Tspan,IC); % Solve ODE

%Solve ODE for NOD

options=odeset('Events',@sickEvents);
 [Tnwave, Ynwave,Te,Ye,ie] =...
     ode15s(@(t,y)fig_2_3_4_5_ODE(t,y,f1n,f2n,wave),Tspan,IC,options); % Solve ODE


Time_sick=Te./7
glucose = Ynwave(:,6); 



%Plot the BIG results with and without the wave

%%%%%%%%%%%%%%%Below we replot using weeks%%%%%%%%%%%%%%%
%Convert days into weeks
 EndTime=EndTime./7;
 Twave=Twave./7;
 Tnwave=Tnwave./7;
 Tnnowave=Tnnowave./7;
 Tnowave=Tnowave./7;
%Plot the BIG results with and without the wave using weeks
figure;
subplot(2,2,1)
semilogy(Tnnowave, Ynnowave(:,5),'-' , 'LineWidth', 2.5,'color', [ 0,.5,.1]); hold on;
semilogy(Tnnowave, Ynnowave(:,6),':' , 'LineWidth', 3  , 'color', [ 1,.1,.1]);
semilogy(Tnnowave, Ynnowave(:,7),'-.', 'LineWidth', 3  , 'color', [ 0,.8,1]);
semilogy(Tnnowave, ones(size(Tnnowave)) * 250, 'LineWidth', 1, 'color', 'k');
titlef = ['A) NOD mouse, no apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',20); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',20);
legend('B (mg)','G (mg/dl)','I (\muU)')
axis([0,EndTime,1,10^3])

subplot(2,2,2)
semilogy(Tnowave, Ynowave(:,5),'-' , 'LineWidth', 2.5, 'color', [ 0,.5,0.1]); hold on;
semilogy(Tnowave, Ynowave(:,6),':' , 'LineWidth', 3  , 'color', [ 1,.1,.1]);
semilogy(Tnowave, Ynowave(:,7),'-.', 'LineWidth', 3  , 'color', [ 0,.8,1]);
semilogy(Tnowave, ones(size(Tnowave)) * 250, 'LineWidth', 1, 'color', 'k');
titlef = ['B) Balb/c mouse, no apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('B','G','I')
axis([0,EndTime,1,10^3])


subplot(2,2,3)
semilogy(Tnwave, Ynwave(:,5),'-' , 'LineWidth', 2.5, 'color', [ 0,.5,.1]); hold on;
semilogy(Tnwave, Ynwave(:,6),':' , 'LineWidth', 3  , 'color', [ 1,.1,.1]);
semilogy(Tnwave, Ynwave(:,7),'-.', 'LineWidth', 3  , 'color', [ 0,.8,1]);
semilogy(Tnwave, ones(size(Tnwave)) * 250, 'LineWidth', 1, 'color', 'k');
%semilogy(Time_sick,250,'*r');
%titlef = '\beta-cell, Glucose, and Insulin';
titlef = ['C) NOD mouse, apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('B','G','I')
axis([0,EndTime,1,10^3])
% 

subplot(2,2,4)
semilogy(Twave, Ywave(:,5),'-' , 'LineWidth', 2.5, 'color', [ 0,.5,.1]); hold on;
semilogy(Twave, Ywave(:,6),':' , 'LineWidth', 3  , 'color', [ 1,.1,.1]);
semilogy(Twave, Ywave(:,7),'-.', 'LineWidth', 3  , 'color', [ 0,.8,1]);
semilogy(Twave, ones(size(Twave)) * 250, 'LineWidth', 1, 'color', 'k');
%titlef = '\beta-cell, Glucose, and Insulin';
titlef = ['D) Balb/c mouse, apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('B','G','I')
axis([0,EndTime,1,10^3])


figure;
%Plot the immune cells with and without the wave
subplot(2,2,1)
semilogy(Tnnowave, Ynnowave(:,10),'-' , 'LineWidth', 2.5, 'color', [ 0.9,.7, 0.1]); hold on;
semilogy(Tnnowave, Ynnowave(:,11),':' , 'LineWidth', 3  , 'color', [.3, 0.6, 0.2]);
semilogy(Tnnowave, Ynnowave(:,12),'-.', 'LineWidth', 3  , 'color', [ .4,.4, .4]);
%titlef = 'Specific Immune Cells';
titlef = ['A) NOD mouse, no apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
legend('E (cells/ml)','R (cells/ml)','Em (cells/ml)')
axis([0,EndTime,1,10^10])

subplot(2,2,2)
semilogy(Tnowave, Ynowave(:,10),'-' , 'LineWidth', 2.5, 'color', [ 0.9,.7, 0.1]); hold on;
semilogy(Tnowave, Ynowave(:,11),':' , 'LineWidth', 3  , 'color', [.3, 0.6, 0.2]);
semilogy(Tnowave, Ynowave(:,12),'-.', 'LineWidth', 3  , 'color', [ 0.4,.4, 0.4]);
%titlef = 'Specific Immune Cells';
titlef = ['B) Balb/c mouse, no apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('E','R','Em')
axis([0,EndTime,1,10^10])

subplot(2,2,3)
semilogy(Tnwave, Ynwave(:,10),'-' , 'LineWidth', 2.5, 'color', [ 0.9,.7, 0.1]); hold on;
semilogy(Tnwave, Ynwave(:,11),':' , 'LineWidth', 3  , 'color', [.3, 0.6, 0.2]);
semilogy(Tnwave, Ynwave(:,12),'-.', 'LineWidth', 3  , 'color', [0.4,0.4,0.4]);
%titlef = 'Specific Immune Cells';
titlef = ['C) NOD mouse, apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('E','R','Em')
axis([0,EndTime,1,10^10])



subplot(2,2,4)
semilogy(Twave, Ywave(:,10),'-' , 'LineWidth', 2.5, 'color', [ 0.9,.7, 0.1]); hold on;
semilogy(Twave, Ywave(:,11),':' , 'LineWidth', 3  , 'color', [.3, 0.6, 0.2]);
semilogy(Twave, Ywave(:,12),'-.', 'LineWidth', 3  , 'color', [ 0.4,.4, 0.4]);
%titlef = 'Specific Immune Cells';
titlef = ['D) Balb/c mouse, apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('E','R','Em')
axis([0,EndTime,1,10^10])


%Plot the dead beta cells and dendritic cells
figure;
subplot(2,2,1)
semilogy(Tnnowave, Ynnowave(:,3),'-' , 'LineWidth', 2.5, 'color', [.3, 0.6, 0.2]);hold on;
semilogy(Tnnowave, Ynnowave(:,4),'-' , 'LineWidth', 3  , 'color', [.9, 0.7,.1]);
semilogy(Tnnowave, Ynnowave(:,9),'--', 'LineWidth', 2.5, 'color', [ .3, 0.6, 0.2]);
semilogy(Tnnowave, Ynnowave(:,8),'-.', 'LineWidth', 3  , 'color', [.9, 0.7, 0.1]);
%titlef = 'Dendritic Cell Dynamics';
titlef = ['A) NOD mouse, no apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
legend('Ba (cells/ml)','Bn (cells/ml)','tD (cells/ml)','D (cells/ml)');
axis([0,EndTime,1,10^10]);

subplot(2,2,2)
semilogy(Tnowave, Ynowave(:,3),'-' , 'LineWidth', 2.5, 'color', [.3, 0.6, 0.2]); hold on;
semilogy(Tnowave, Ynowave(:,4),'-' , 'LineWidth', 3, 'color', [.9, 0.7,.1]);
semilogy(Tnowave, Ynowave(:,9),'-.', 'LineWidth', 3  , 'color', [.3, 0.6, 0.2]);
semilogy(Tnowave, Ynowave(:,8),'--', 'LineWidth', 2.5, 'color', [.9, 0.7,.1]);
%titlef = 'Dendritic Cell Dynamics';
titlef = ['B) Balb/c mouse, no apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('Ba','Bn','D','tD');
axis([0,EndTime,1,10^10]);

subplot(2,2,3)
semilogy(Tnwave, Ynwave(:,3),'-' , 'LineWidth', 2.5, 'color', [.3, 0.6, 0.2]);hold on;
semilogy(Tnwave, Ynwave(:,4),'-' , 'LineWidth', 3  , 'color', [.9, 0.7,.1]);
semilogy(Tnwave, Ynwave(:,9),'-.', 'LineWidth', 3  , 'color', [.3, 0.6, 0.2]); 
semilogy(Tnwave, Ynwave(:,8),'--', 'LineWidth', 2.5, 'color',[.9, 0.7,.1]);
%titlef = 'Dendritic Cell Dynamics';
titlef = ['C) NOD mouse, apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('Ba','Bn','D','tD');
axis([0,EndTime,1,10^10]);

subplot(2,2,4)
semilogy(Twave, Ywave(:,3),'-' , 'LineWidth', 2.5, 'color', [.3, 0.6, 0.2]); hold on;
semilogy(Twave, Ywave(:,4),'-' , 'LineWidth', 3, 'color', [.9, 0.7,.1]);
semilogy(Twave, Ywave(:,9),'-.', 'LineWidth', 3  , 'color', [.3, 0.6, 0.2]);
semilogy(Twave, Ywave(:,8),'--', 'LineWidth', 2.5, 'color', [.9, 0.7,.1]);
%titlef = 'Dendritic Cell Dynamics';
titlef = ['D) Balb/c mouse, apoptotic wave'];
title(titlef,'FontSize',20);
xlabel('Time (weeks)','FontSize',22); 
yalabel = 'Species population';
yalabel = [yalabel sprintf('\n') '(Log scale)'];
ylabel(yalabel,'FontSize',22);
%legend('Ba','Bn','D','tD');
axis([0,EndTime,1,10^10]);

% 
%Extract and plot data
for i = 1: 11
    %load data
dat = csvread(['Lietal_data/dat',num2str(i),'.csv']) ;
    % find when mouse gets sick 
    
week(i)= min(dat(dat(:,2)>250,1));

figure(4);
plot(dat(:,1), dat(:,2),'*-');
hold on; 
plot(Tnwave, Ynwave(:,6),':' , 'LineWidth', 3  , 'color', 'r');
hold on
plot(Tnwave, ones(size(Tnwave)) * 250, 'LineWidth', 1, 'color', 'k');
clear dat
titlef = 'Comparison of Model Glucose Levels for NOD with apoptotic wave';
title(titlef,'FontSize',20);
xlabel('Time (weeks)', 'FontSize',22); 
yalabel = 'Glucose Levels';
ylabel(yalabel,'FontSize',22);
%legend('B','G','I')
axis([0,40,90,600])


end

