function val=least_squares(k, mydata, mytime)

y0(1) = 21.5; y0(2) = 3.4;
modelfun=@(t,x)Lotka_Volterra_Model(t,x,k);
[t ycalc]=ode45(modelfun,mytime,y0);
resid = (ycalc-mydata).*(ycalc-mydata);
val = sum(sum(resid));
end