%Function to produce histograms of parameters seen in write up and also to
%fit normal distributions to the parameter values.

clear all;
file = 'all_values_IMPORTANT.csv';
vals = readtable(file);
vals = vals{:,:};


vals = vals(2:42,:);
numParams = 41;
dist_properties = zeros(41,2);
param_important = [5 6 10 15 18 34 39 40 41]
param_names = ["e_1"; "e_2"; "\delta_B"; "S_I"; "G_I"; "\mu"; "\alpha_{\eta}"; "\beta_{\eta}"; "\eta"];
fig = figure(1);
for i = 1:9
    param_num = param_important(i);
    param_set = vals(param_num,[2 3 4 6 7 8 9 10 11]);
    
    subplot(3,3,i);
    histogram(param_set, 10);
    title(param_names(i, :));
    
    pd = fitdist(param_set', 'Normal');
    dist_properties(i,1) = pd.mu;
    dist_properties(i,2) = pd.sigma;
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han, 'frequency', 'fontsize', 15);
xlabel(han, 'parameter value', 'fontsize', 15);
yh = get(gca,'ylabel'); % handle to the label object
p = get(yh,'position'); % get the current position property
p(1) = 1.5*p(1);
set(yh,'position',p);
dist_properties
