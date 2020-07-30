%FUNCTION TO CALCULATE THE LEAST SQUARES ERROR OF A MOUSE FIT TO RAW DATA

function [error] = determineError(X,final_params, wave)

fig_2_3_4_5_Parameters; %load parameters to get f1n and f2n
X(1,:) = X(1,:) * 7; %rescale time to days
    
T = 350; %will simulate through 350 days


InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
final_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, final_params), [0 350], InitialState,options); %Simulate ODE
x_predicted = deval(final_sol, X(1,:)); %get solutions to ODE at given time point
%calculate error of prediction as comapred to raw data
error = X(2,:) - x_predicted(6,:);
error = error.^2;
error = sum(error);
