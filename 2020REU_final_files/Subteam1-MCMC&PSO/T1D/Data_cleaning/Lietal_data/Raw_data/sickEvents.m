%% Extract post-UKF parameter means/vars
%  Author:       Shtylla et al
%  Date:         2019
%  Desc:         Function to determine when a mouse has achieved diabetes

function [value, isterminal,direction]=sickEvents(t,y)
    
    value = y(6)-250; % Detect glucose=250
    isterminal = 0; % Do not stop integration
    direction = 1; % The zero can be approached from left to right only (positive direction)

end