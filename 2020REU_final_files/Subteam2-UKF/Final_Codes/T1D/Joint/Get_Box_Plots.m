  
%Script to produce box plots of the acute mice key parameter values


clear all;
file = 'all_values_IMPORTANT_multiple_eta.csv';
vals = readtable(file);
vals = vals{:,2:end};


vals = vals(2:42,:);
numParams = 41;
dist_properties = zeros(41,2);
%param_important = [5 6 10 15 18 34 39 40 41]
param_important = [5 6 15 18 34 40];
param_names = ["e_1"; "e_2"; "SI"; "GI"; "\mu"; "\alpha_{\eta}"];
fig = figure(1);
for i = 1:6
    param_num = param_important(i);
    param_set = vals(param_num,[1 2 3 4 5 6 7 8 9 10 11]);
    
    subplot(2,3,i);
    boxplot(param_set);
    set(gca,'XTickLabel',{' '})
    title(param_names(i), 'fontsize', 15);
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han, 'parameter value', 'fontsize', 15);
yh = get(gca,'ylabel'); % handle to the label object
p = get(yh,'position'); % get the current position property
p(1) = 1.5*p(1);
set(yh,'position',p);