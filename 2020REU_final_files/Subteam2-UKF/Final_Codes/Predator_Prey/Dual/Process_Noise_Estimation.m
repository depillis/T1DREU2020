%FUNCTION TO ESTIMATE MEASUREMENT NOISE

HLData = load('HaresLynxData.mat');  %Load dataset
rawData = HLData.Lotka_Volterra_Data;
x(1:2,:) = rawData(:, 2:3)'; %The real data
Nx = 2;
T = 91;
x_fake = zeros(Nx,T);
x_fake(:,1) = [x(1,1); x(2,1)];
tspan = [0 T];
sol = ode45(@(t, y) Lotka_Volterra_Model(t, y, truepar), tspan, x_fake(:,1)); %Use ODE solver
for t=2:T %Generate data for each time point
  x_fake(:,t) = deval(sol, t); %Get ODE solution for time t
end

prey_noise = x(1,:) - x_fake(1,:); %Difference between real and ODE
predator_noise = x(2,:) - x_fake(2,:); 
cov(prey_noise, predator_noise) %Print covariance matrix

figure(1)
p1 = histogram(prey_noise);
xlabel('Prey Process Error', 'Fontsize', 15);
ylabel('Frequency', 'Fontsize', 15);

figure(2)
p1 = histogram(predator_noise);
xlabel('Predator Process Error', 'Fontsize', 15);
ylabel('Frequency', 'Fontsize', 15);
mean(prey_noise)
mean(predator_noise)






