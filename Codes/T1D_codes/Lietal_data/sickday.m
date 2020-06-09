% find days when mouse gets sick 

for i = 1: 11
    %load data
dat = csvread(['dat',num2str(i),'.csv']) ;
    % find when mouse gets sick 
    
week(i)= min(dat(dat(:,2)>250,1));

figure(1);
plot(dat(:,1), dat(:,2),'*-');
xlabel('Weeks', 'FontSize', 30)
ylabel('Glucose Level', 'FontSize', 30)
hold on; 

clear dat


end

mean(week)

figure(2); 
hist(week)
xlabel('Weeks', 'FontSize', 30)
ylabel('Glucose Level', 'FontSize', 30)