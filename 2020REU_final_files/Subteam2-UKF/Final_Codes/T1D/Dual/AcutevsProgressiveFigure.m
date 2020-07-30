%CODE TO PRODUCE FIGURE DEMONSTRATING DIFFERENCE BETWEEN LI ACUTE AND
%PROGRESSIVE MICE

figure(1)
for i = 1:11
    i_str = int2str(i);
    file = strcat('dat',i_str,'.csv'); %create file name
    X = readtable(file);
    X = X{:,:};
    X = flip(X);
    X = X'; %set up table as wide
    if i == 1 || i == 5 
        p = plot(X(1,:), X(2,:), 'r', 'Linewidth', 2); hold on;
        if i == 5
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    else
        p = plot(X(1,:), X(2,:), 'b', 'Linewidth', 2); hold on;
        if p == 2
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
end
hold off;
ylabel('Glucose', 'fontsize', 15);
xlabel('Week', 'fontsize', 15);
legend('Progressive', 'Acute');
