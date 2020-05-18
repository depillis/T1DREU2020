function y = LormeasO(x,u,n,t,par)
%Par is not used but kept here so don't have to modify the generic UKF
%algorithm

if ~isempty(n)
y(1:3,:) = x(1:3,:) + n(1:3,:);

else
y(1:3,:) =x(1:3,:);
end