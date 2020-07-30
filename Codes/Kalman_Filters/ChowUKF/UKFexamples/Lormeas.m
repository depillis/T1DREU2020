function y = Lormeas(x,u,n,t,par)
%Par is not used but kept here so don't have to modify the generic UKF
%algorithm

y(1:3,:) = x(1:3,:) + x(7:9,:);
