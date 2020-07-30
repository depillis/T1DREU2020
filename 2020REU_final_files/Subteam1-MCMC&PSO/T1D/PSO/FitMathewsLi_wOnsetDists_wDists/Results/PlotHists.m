%% Plotting histograms and finding mean, median, mode params
%  Author:       C. Catlett
%  Date:         July 1, 2020
%  Desc:         Load results, make histograms for each 

load('paramsdist.mat');
load('paramsdist2.mat');
load('testDists.mat'); % loads as 'ans'
testDists = ans;

all = [paramsDist; paramsDist2; testDists];
modelst = mode(all(:,1));
medlst = median(all(:,1));
meanlst = mean(all(:,1));

% Find mean, median, mode of all params
for n=2:60
    modelst = [modelst mode(all(:,n))];
    medlst = [medlst median(all(:,n))];
    meanlst = [meanlst mean(all(:,n))];
end

% eFAST sensitive params
inds = [3 10 38 39 27 28 26];
for k=inds
histogram(all(:,k));
pause
end