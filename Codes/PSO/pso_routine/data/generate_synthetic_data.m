function y_sample= generate_synthetic_data
%N_sample is the number of sample points used to sample the ODE
clear all;
close all;

N_sample = 15; %Pick the total number of time points
N=65; %Number of reactions
t0= 0;
scale_factor=1;
tfinal = scale_factor*60;
meas_error= 0.005;

global stimulated;

t_interpolate= linspace(0,tfinal,N_sample);

%Solve the ODE, with stimulated conditions of IGF-1
%and sample the points with added error

stimulated=1;

%Enter the model parameters

reduced_model_parameters;

[t,y] = ode15s(@IGFRODE,[0 tfinal],init_vec);

[t_sample,y_sample] = ode15s(@IGFRODE,t_interpolate, init_vec);

y_sample = y_sample + meas_error*randn(size(y_sample));

y_sample=abs(y_sample);%Take the absolute value to avoid negative numbers

%save the output
csvwrite('synthetic_data_points_Nsample15.csv',y_sample);
csvwrite('synthetic_data_time_Nsample15.csv',t_interpolate);


end