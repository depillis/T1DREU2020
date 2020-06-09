%Author: B. Shtylla
%Last Updated: 05/07/2020

%This code loads the results from the PSO or any fitting subroutine and
%the data that the model was fitted against and then plots them together. 


%load the data points
data_points = csvread('data/synthetic_data_points_Nsample15.csv');
%load the data times
data_time = csvread('data/synthetic_data_time_Nsample15.csv');

%Load the fitted parameters and assign them to the paramter vector and
%script
P_fitted = csvread('pso_param_set1_Nsample15.csv');

global P;
P=P_fitted;
pso_parameters;

%Declare the fixed initial conditions for ODE solver
pso_init_vec;
      
      
 % Solve the ODE with the fitted paramters
  
 [fitted_t, fitted_y]=ode15s(@IGFRODE,[0 max(data_time)],init_vec);
    
 %Plot the results
 %We need to plot the following indexed-terms from the solution matrix      
 ipAKT=22;
 ipGSK3 =28;
 ipMAPK=19;
 ipmTOR=35;
 ipP70S6K=37;
 ipTSC2=31;
 
 %Assign the variables for the the stimulated result plots
 y=fitted_y;
 t=fitted_t;
 yp=data_points;
 tp=data_time;
 
 %Generate the plot
figure(1)
subplot(2,3,1)
plot(t,y(:,ipAKT),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipAKT),'*b','Linewidth',3)
hold off;
title('p-AKT:stimulated');
xlabel('Time t');
ylabel('Solution y');
legend('p-AKT','in silico p-AKT')

subplot(2,3,2)
plot(t,y(:,ipGSK3),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipGSK3),'*b','Linewidth',3)
hold off;
title('p-GSK3: stimulated');
xlabel('Time t');
ylabel('Solution y');
legend('p-GSK3','in silico p-GSK3')

subplot(2,3,3)
plot(t,y(:,ipMAPK),'-r','Linewidth',3) 
hold on
plot(tp,yp(:,ipMAPK),'*b','Linewidth',3) 
hold off;
title('p-MAPK:stimulated');
xlabel('Time t');
ylabel('Solution y');
legend('p-MAPK','in silico p-MAPK')

subplot(2,3,4)
plot(t,y(:,ipmTOR),'-r','Linewidth',3) 
hold on;
plot(tp,yp(:,ipmTOR),'*b','Linewidth',3)
hold off;
hold off;
title('p-mTOR:stimulated');
xlabel('Time t');
ylabel('Solution y');
legend('p-mTOR','in silico p-mTOR')

subplot(2,3,5)
plot(t,y(:,ipP70S6K),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipP70S6K),'*b','Linewidth',3)
hold off;
title('p-p7SK06:stimulated');
xlabel('Time t');
ylabel('Solution y');
legend('p-p7SK06','in silico p-p7SK06')

subplot(2,3,6)
plot(t,y(:,ipTSC2),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipTSC2),'*b','Linewidth',3)
hold off;
title('p-TSC2:stimulated');
xlabel('Time t');
ylabel('Solution y');
legend('p-TSC2','in silico p-TSC2')