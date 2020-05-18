function dy = Lotka_Volterra_Model(t,y,k)
% Lotka-Volterra predator-prey model.
dy = zeros(2,1);
alpha = k(1); beta = k(2); gamma = k(3); delta = k(4);
A = [alpha - beta*y(2), 0; 0, ...
    -gamma + delta*y(1)];
dy = A*y;
end