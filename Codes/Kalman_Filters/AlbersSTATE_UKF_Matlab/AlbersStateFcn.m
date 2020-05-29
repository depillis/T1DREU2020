function [xnew] = AlbersStateFcn(x)
% Author: LdP and CLe
% LdeP Update March 2020: Converting to using MATLAB ODE solver (instead of Euler
% step)
% Date: November 19, 2019
% Summary:  Euler integration of continuous-time
%           dynamics x'=f(x) w/ sample time dt

% Inputs:   xk = states, x[k]
% Outputs:  xk1 = propagated states, x[k+1]

%LdeP comment: The following is equivalent to the simplest Euler-style
%integration. Instead, we should make use of a built-in ODE solver.
 dt = 0.1;
%%% t = 0;
%%% x = x + AlbersODE(t, x)*dt;

% %LdeP Convert scale from glucose/dL back to total glucose (multiply by
% 100)
x(3)=100.*x(3);

% LdeP comment
% x is the current state (at "time" k-1). Use ODE45 to get x at the next
% "time" k.
initialStateGuess = x; %Current state vector at time k

timeVector = 0:dt:1; % Assuming a time-invariant system, the step from k-1 to k should be independent of t, so we should always be a able to go from 0 to 1.
[timeSteps,xout]=ode45(@AlbersODE,timeVector,initialStateGuess);
%[timeSteps,xout]=ode23s(@AlbersODE,timeVector,initialStateGuess); %LdeP Stiff?


%LdeP Convert scale back to glucose per dL (divide by 100)
% 
xout(:,3)=xout(:,3)./100;
xnew = xout(end,:)'; % LdeP Just return the final time solution


end