%% Function to determine progressive vs acute mice

function [acute, progressive] = data_category
% read in data file
raw=readmatrix('acute+progressive.csv');

% gather data
time=raw(3, 4:end); % day 0:25 = day -24:0 in Mathews paper
glucose=raw(4:end,4:end); 

% loop to find diabetes onset
% define acute: <200 mg/dL glucose prior to onset (day 25)
% define progressive: >= 200 mg/dL glucose prior to onset

acute=zeros(1,25);
progressive=zeros(1,25);

for i = 1:489 % for each mouse
    count=0;
    for j = 1:23 % check glucose reading in each day
      if isnan(glucose(i,j)) || glucose(i,j) < 200 %if it is less than 200 add to the acute counter
          count=count+1;
      end
    end
    if count==23 % if we documented 25 readings where the data was <200 add it to the acute mice
         acute=vertcat(acute, glucose(i,:));
    else
         progressive=vertcat(progressive, glucose(i,:));
    end
end
         
acute=acute(2:end,:);
progressive=progressive(2:end,:);

% plot
% figure(1); clf
% subplot(2,1,1)
% for i = 1:132
%   hold on
%   plot(time, acute(i,:))
% end
% subplot(2,1,2)
% for j = 1:357
%     hold on
%     plot(time, progressive(i,:))
% end
end