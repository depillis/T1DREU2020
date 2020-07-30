%% Function to determine progressive or acute mouse type
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         Determine mouse type, not fully working. Not used in
%                current DRAM execution
%function [acute progressive] = det_mouse_type
% for n = 1:11 %loop over mice 
%     n_str = int2str(n);
%     file = strcat('dat',n_str,'.csv'); %create file name
%     X = readtable(file);
%     X = X{:,:};
%     X = flip(X);
%     %X = X'; %set up table as wide
%     
%     count=0;
%     numReadings=length(X);
%     for j = 1:numReadings % check glucose reading in each day
%       if isnan(X(j,2)) || X(j,2) < 200 %if it is less than 200 add to the acute counter
%           count=count+1;
%       end
%     end
%     if count==numReadings % if we documented 25 readings where the data was <200 add it to the acute mice
%          acute=vertcat(acute, glucose(i,:));
%     else
%          progressive=vertcat(progressive, glucose(i,:));
%     end
%     
% 
% end
% 
% acute=zeros(1,25);
% progressive=zeros(1,25);
% 
% end
