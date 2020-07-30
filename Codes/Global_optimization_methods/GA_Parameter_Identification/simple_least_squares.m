function val = simple_least_squares(k)
a=k(1); b=k(2);
mytime=(0:6)';
mydata(:,1) = [-1 2 11 26 47 74 107];
mydata(:,2) = [ 1 3 09 19 33 51 73];
y0(1) = -1; y0(2) = 1;
modelfun=@(t,x)simple_model(t,x,a,b);
[t ycalc]=ode45(modelfun,mytime,y0);
resid = (ycalc-mydata).*(ycalc-mydata);
val = sum(sum(resid));