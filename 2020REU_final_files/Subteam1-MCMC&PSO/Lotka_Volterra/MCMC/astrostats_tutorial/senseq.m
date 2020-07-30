%--------------------------------------------------------------------------
% N. Tania (May 30, 2018)
%
% Computes the sensitivity matrix dy/dpars and 
% the relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, 
%
% senseq measures sensitivity relative to parameter
% senseq_ic measures sensitivity relative to initial conditions
%
% The code below are modified from
%   https://web.math.ncsu.edu/cdg/supplemental-material/
%--------------------------------------------------------------------------

function [sens, sensrel, flagg, y, sol] = senseq(pars,x,Init)

%[sol,y,flagg] = model_sol(pars,x,Init);
obspops = [Init(:,1) Init(:,2)]; % Collected data for population size, 2x90 matrix
pops_init = [Init(1,1) Init(1,2)]';
[y,sol] = ode45(@(t,y) Lotka_Volterra_Model(t,y,pars),x, pops_init);
y=0;

%Depending on form of y, the sensitivity matrix could be dymodel/dpars or 
%dr/dpars where r = ymodel - ydata

% flagg is to catch integration failures
if y == 0
[sens, sensrel] = diffjac(pars,@myfun,y,x,Init);
else
 sens = 0; senrel=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function y = myfun(pars,x,Init)

%[sol,y,flagg]  = model_sol(pars,x,Init);
obspops = [Init(:,1) Init(:,2)]; % Collected data for population size, 2x90 matrix
pops_init = [Init(1,1) Init(1,2)]';
[y,sol] = ode45(@(t,y) Lotka_Volterra_Model(t,y,pars),x, pops_init);
ye=0;
    
function [jac, jacrel] = diffjac(x,f,f0,time,Init)
% Compute a forward difference dense Jacobian f'(x), return lu factors.
%
% uses dirder
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
%
% inputs:
%          x, f = point and function
%          f0   = f(x), preevaluated

n = length(x);
jac = zeros(length(f0), n);
jacrel = zeros(length(f0), n);
for j = 1:n
    zz = zeros(n,1);
    zz(j) = 1;
    [jac(:,j), jacrel(:,j)] = dirder(x,zz,f,f0,time,Init);
end

function [z, zrel] = dirder(x,w,f,f0,time,Init)
% Compute a finite difference directional derivative.
% Approximate f'(x) w
% 
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function z = dirder(x,w,f,f0)
%
% inputs:
%           x, w = point and direction
%           f = function
%           f0 = f(x), in nonlinear iterations
%                f(x) has usually been computed
%                before the call to dirder

% Hardwired difference increment.
global DIFF_INC

epsnew = DIFF_INC;

n = length(x);

% scale the step
if norm(w) == 0
    z = zeros(n,1);
return
end

% Now scale the difference increment.

xs=(mtimes(x, w))/norm(w);
if xs ~= 0.d0
    epsnew=epsnew*max(abs(xs),1.d0)*sign(xs);
end
epsnew=epsnew/norm(w);

% del and f1 could share the same space if storage
% is more important than clarity.
del = x+epsnew.*w;
f1  = feval(f,del,time,Init);
z   = (f1 - f0)/epsnew;

% %Relative sensitivity
zrel = z;
gh = find(w == 1);
zrel = zrel*x(gh)./f0;
