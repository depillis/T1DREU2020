function dx= Lotka_Volterra_Model(t,x,k)
% Lotka-Volterra predator-prey model.
dx = zeros(6,1);
alpha = k(1); 
beta = k(3); 
gamma = k(2); 
delta = k(4);

A = [alpha - beta*x(2), 0; 0, ...
    -gamma + delta*x(1)];
dx(1:2) = A*x(1:2);
dx(3:6) = [0;0;0;0];
end