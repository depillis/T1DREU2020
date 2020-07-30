% save data into file 

for i = 1: 9
dat = [eval(['Data00', num2str(i),'(:,1)'])...
    eval(['Data00', num2str(i),'(:,2)'])];
csvwrite(['dat',num2str(i),'.csv'], dat) 
end

for i = 10:11
dat = [eval(['Data0', num2str(i),'(:,1)'])...
    eval(['Data0', num2str(i),'(:,2)'])];
csvwrite(['dat',num2str(i),'.csv'], dat) 
end


% find days when mouse gets sick 

for i = 1: 11
    %load data
dat = csvread(['dat',num2str(i),'.csv']) ;
    % find when mouse gets sick 
day(i)= min(dat(dat(:,2)>250,1));
clear dat
end
