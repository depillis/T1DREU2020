%Script to produce figures of all various fits. Alterations can be made
%within this file to produce the various figures in the write up. We also
%print out MSE values throughout this file.


clear all;
file = 'all_values_IMPORTANT.csv';
vals = readtable(file);
vals = vals{:,:};
T = 350;
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
fig_2_3_4_5_Parameters;
wave = wave_basal;
fig = figure(1);
count = 1;
for m = 1:11
    if m ~= 1 && m ~= 5
    m_str = int2str(m);
    mouse_dat_file = strcat('dat', m_str,'.csv');
    mouse_dat = readtable(mouse_dat_file);
    mouse_dat = mouse_dat{:,:};
    mouse_dat = flip(mouse_dat);
    mouse_dat = mouse_dat';
    mouse_dat(1,:) = mouse_dat(1,:) * 7;
    num_glucose = size(mouse_dat, 2);
    mouse_params = vals(:,m);
    least_squares = mouse_params(1);
    MSE = least_squares / num_glucose;
    m
    MSE
    sim_params = mouse_params(2:42);
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    sim_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, sim_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
    x_sim = zeros(12, 350);
    for j = 1:350
        x_sim(:,j) = deval(sim_sol, j);
    end
    subplot(3,3,count);
    p1 = plot(1:T, x_sim(6,:), 'r', 'Linewidth', 2); hold on;
    p2 = plot(mouse_dat(1,:), mouse_dat(2,:), 'ko', 'MarkerFaceColor', 'k', 'Linewidth', 2); hold off;
    p2.MarkerSize = 3;
    sub_title = strcat('Mouse ', m_str);
    title(sub_title);
    count = count + 1;
    end
    
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han, 'glucose(mg/dl)', 'fontsize', 15);
xlabel(han, 'day', 'fontsize', 15);
yh = get(gca,'ylabel'); % handle to the label object
p = get(yh,'position'); % get the current position property
p(1) = 1.5*p(1);
set(yh,'position',p);  