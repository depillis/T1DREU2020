function y = MeasurementFcn(x,u,n,t,par)
%Par is not used but kept here so don't have to modify the generic UKF
%algorithm

%This describes the measurement function (often denoted as H) for the UKF 

if ~isempty(n)
y(1:2,:) = x(1:2,:) + n(1:2,:);
%y(1:6,:) = x(1:6,:) + n(1:6,:);
else
y(1:2,:) = x(1:2,:);
%y(1:6,:) = x(1:6,:);
end

end

