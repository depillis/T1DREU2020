function newx =LordynO(x,u,v,t,par)
%For the fourth-order Runge-Kutta below I used  delta = .01
%.167 below comes from delta/6
%x1-x3 are the 3 state variables of the Lorenz system
%x5-x6 are the model parameters and
%x7-x9 are the measurement errors associated with y1-y3 at time t

if ~isempty(x)
par= x(4:6,:);
k1 = Lorenz(x(1:3,:), par);
k2 = Lorenz(x(1:3,:) +.005*k1, par);
k3 = Lorenz(x(1:3,:) +.005*k2, par);
k4 = Lorenz(x(1:3,:) +.01*k3, par);
newx(1:3,:) = x(1:3,:) + (.0017)*(k1 + 2*k2 + 2*k3 + k4);

newx(4:6,:) = x(4:6,:);
end    

if ~isempty(v)
   newx(1:6,:) = newx(1:6,:) + v(1:6,:);
end