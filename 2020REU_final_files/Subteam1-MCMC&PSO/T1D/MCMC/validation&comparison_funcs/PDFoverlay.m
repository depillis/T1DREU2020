%% Plot parameter PDFs from multiple T1D parameterization scenarios 
%  Author:       M. Watanabe
%  Date:         July 2020
%  Desc:         Script to plot overlaid parameter PDFs from 
%                parameterizations using from Mice 3, 6, 11 and 
%                averaged data


% load past work spaces
clear all
addpath('../past_run_workspaces');

% Average 
load('jul10_avg_run1(noIC)_acute_NOD_waveOn_lietal');
chain_avg = chain; % red

% Mouse 6
load('jul10_mouse6_run1(noIC)_acute_NOD_waveOn_lietal');
chain6 = chain; % blue

% Mouse 3
load('jul13_mouse3_acute_NOD_waveOn_lietal');
chain3 = chain; % green

% Mouse 11
load('jul13_mouse11_acute_NOD_waveOn_lietal');
chain11 = chain; % magenta

titles{1} = 'e1';
titles{2} = 'e2';
titles{3} = '\delta_B';
titles{4} = 'SI';
titles{5} = 'GI';
titles{6} = '\mu_R';
titles{7} = '\mu_E';
titles{8} = '\alpha_\eta';
titles{9} = '\beta_\eta';
titles{10} = '\eta';


figure(1)
for i = 1:10
h=subplot(4,3,i);
    
    % Average
    [y,x]=density(chain_avg(:,i),[],2);
    plot(x,y,'-k', 'color', [0.46 0.72 0.58],'LineWidth', 2)
    hold on
    
    % Mouse 3
    [y,x]=density(chain3(:,i),[],2);
    plot(x,y,'-k', 'color', [0.31 0.58 0.83],'LineWidth', 2);
    hold on
    
    % Mouse 6
    [y,x]=density(chain6(:,i),[],2);
    plot(x,y,'-k', 'color', [0.87, 0.31, 0.31],'LineWidth', 2);
    hold on
    
    % Mouse 11
    [y,x]=density(chain11(:,i),[],2);
    plot(x,y,'-k', 'color', [0.62 0.45 0.71],'LineWidth', 2);
     
    set(gca, 'FontSize', 15);
    title(titles{i});
end

han=axes(figure(1),'visible','off');
han.XLabel.Visible='on';
xlabel('parameter value')
set(gca, 'FontSize', 20)
hold off
