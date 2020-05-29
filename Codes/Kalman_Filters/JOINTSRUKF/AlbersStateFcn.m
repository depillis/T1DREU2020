function [xnew] = AlbersStateFcn(x,DelT)
% Author: LdP and CLe
% LdeP March 31, 2020 - modify to allow for sample time input (instead of
% assuming that we just step from k to k+1, we go from k to k+DelT).
%
% LdeP Update April 2, 2020 - Allow for 3 appended parameters 
% These don't need to go through the ODE since they will just be propagated as they are
% 
% LdeP Update March 2020: Converting to using MATLAB ODE solver (instead of Euler
% step)
% Date: November 19, 2019
% Summary:  Euler integration of continuous-time
%           dynamics x'=f(x) w/ sample time dt

% Inputs:   xk = states, x[k] - LdeP Update: There will be 9 entries in x (6 States 3 Params)
% Outputs:  xk1 = propagated states, x[k+1]

%LdeP comment: The following is equivalent to the simplest Euler-style
%integration. Instead, we should make use of a built-in ODE solver.
%dt = 0.1;
%%% t = 0;
%%% x = x + AlbersODE(t, x)*dt;

% %LdeP Convert scale from glucose/dL back to total glucose (multiply by 100)
x(3)=100.*x(3);

% LdeP comment
% x is the current state (at "time" k-1). Use ODE45 to get x at the next
% "time" k.
initialStateGuess = x; %Current state vector at time k, use all 9 states (last 3 are params)

%timeVector = 0:dt:1; % Assuming a time-invariant system, the step from k-1 to k should be independent of t, so we should always be a able to go from 0 to 1.
%LdeP modified time vector to allow for user-defined step size from previous k
%Since we assume time invariance, we can assume we move from time 0 to time dt.

% CALL ODE SOLVER TO GET SOLUTION FOR x at next time step
%timeVector=0:dt:DelT;
%[timeSteps,xout]=ode45(@AlbersODE,timeVector,initialStateGuess);
%[timeSteps,xout]=ode23s(@AlbersODE,timeVector,initialStateGuess); %LdeP Stiff?
% LdeP just send through initial and final time.
% [timeSteps,xout]=ode45(@AlbersODE,[0 DelT],initialStateGuess);
[timeSteps,xout]=ode45(@AlbersODEwParams,[0 DelT],initialStateGuess);


%LdeP Convert scale back to glucose per dL (divide by 100)
% 
xout(:,3)=xout(:,3)./100;
xnew = xout(end,:)'; % LdeP Just return the final time solution

%%LdeP Append the parameters that were sent in.
%xnew = [xnew x(7:9)];

end
