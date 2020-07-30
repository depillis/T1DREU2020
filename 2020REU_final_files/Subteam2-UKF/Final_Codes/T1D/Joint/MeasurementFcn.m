function y = MeasurementFcn(x,u,n,t,par)
%Par is not used but kept here so don't have to modify the generic UKF
%algorithm

%This describes the measurement function (often denoted as H) for the UKF

%Changed so only glucose

if ~isempty(n)
y(1,:) = x(6,:) + n(6,:);

else
y(1,:) = x(6,:);
end

end