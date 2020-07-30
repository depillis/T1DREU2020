function newx =Lordyn(x,u,v,t,par)
%For the fourth-order Runge-Kutta below I used  delta = .01
%.167 below comes from delta/6
%x1-x3 are the 3 state variables of the Lorenz system
%par holds the true parameter values needed for the Runge-Kutta method

if ~isempty(x)
par= x(4:6,:);
k1 = Lorenz(t, x(1:3,:), par);
k2 = Lorenz(t+.05, x(1:3,:) +.005*k1, par);
k3 = Lorenz(t+.05, x(1:3,:) +.005*k2, par);
k4 = Lorenz(t+.1, x(1:3,:) +.01*k3, par);
newx(1:3,:) = x(1:3,:) + (.0017)*(k1 + 2*k2 + 2*k3 + k4);

newx(4:9,:) = x(4:9,:);
end    

if ~isempty(v)
   newx(1:9,:) = newx(1:9,:) + v(1:9,:);
end