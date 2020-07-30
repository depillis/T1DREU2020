file = 'all_values_important_multiple_eta.csv'; %contains all params for all mice
vals = readtable(file); %
vals = vals{:,2:12};
T = 350;
T1D_ODE_Parameters;
wave = wave_basal;
    
%Produce figure for acute mice
fig = figure(1);
count = 1;
for m = 1:11
    %If you wish to run this on mice 9-11 make sure you set eta to 0.018
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
    %least_squares = mouse_params(1);
    %MSE = least_squares / num_glucose;
    %m
    %MSE
    sim_params = mouse_params(2:42);
    InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
    InitialState = [InitialState; sim_params];
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    if m < 9
        sim_sol = ode15s(@(t, y) ODE_allparams(t, y, f1n, f2n, wave, sim_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
    else
        sim_sol = ode15s(@(t, y) ODE_allparams2(t, y, f1n, f2n, wave, sim_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
    end
    x_sim = zeros(12, 350);
    for j = 1:350
        temp = deval(sim_sol, j);
        x_sim(:,j) = temp(1:12);
    end
    subplot(3,3,count);
    p1 = plot(1:T, x_sim(6,:), 'r', 'Linewidth', 2); hold on;
    p2 = plot(mouse_dat(1,:), mouse_dat(2,:), 'ko', 'Linewidth', 2); hold off;
    p2.MarkerSize = 3;
    p2.MarkerFaceColor = 'k';
    %ylabel('Glucose', 'fontsize', 15);
    %xlabel('Day', 'fontsize', 15);
    sub_title = strcat('Mouse ', m_str);
    title(sub_title);
    count = count + 1;
    end
    
    han = axes(fig, 'visible', 'off');
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabel(han, 'glucose(mg/dl)', 'fontsize', 15);
    xlabel(han, 'day', 'fontsize', 15);
    yh = get(gca, 'ylabel');
    p = get(yh, 'position');
    p(1) = 1.5*p(1);
    set(yh, 'position', p);
end

%Produce figure for progressive mice
fig2 = figure(2);
count = 1;
for m = 1:11
    if m == 1 || m == 5
    m_str = int2str(m);
    mouse_dat_file = strcat('dat', m_str,'.csv');
    mouse_dat = readtable(mouse_dat_file);
    mouse_dat = mouse_dat{:,:};
    mouse_dat = flip(mouse_dat);
    mouse_dat = mouse_dat';
    mouse_dat(1,:) = mouse_dat(1,:) * 7;
    num_glucose = size(mouse_dat, 2);
    mouse_params = vals(:,m);
    %least_squares = mouse_params(1);
    %MSE = least_squares / num_glucose;
    %m
    %MSE
    sim_params = mouse_params(2:42);
    InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
    InitialState = [InitialState; sim_params];
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    if m < 9
        sim_sol = ode15s(@(t, y) ODE_allparams(t, y, f1n, f2n, wave, sim_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
    else
        sim_sol = ode15s(@(t, y) ODE_allparams2(t, y, f1n, f2n, wave, sim_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
    end
    x_sim = zeros(12, 350);
    for j = 1:350
        temp = deval(sim_sol, j);
        x_sim(:,j) = temp(1:12);
    end
    subplot(1,2,count);
    p1 = plot(1:T, x_sim(6,:), 'r', 'Linewidth', 2); hold on;
    p2 = plot(mouse_dat(1,:), mouse_dat(2,:), 'ko', 'Linewidth', 2); hold off;
    p2.MarkerSize = 3;
    p2.MarkerFaceColor = 'k';
    sub_title = strcat('Mouse ', m_str);
    title(sub_title);
    count = count + 1;
    end
    
    han2 = axes(fig2, 'visible', 'off');
    han2.Title.Visible = 'on';
    han2.XLabel.Visible = 'on';
    han2.YLabel.Visible = 'on';
    ylabel(han2, 'glucose(mg/dl)', 'fontsize', 15);
    xlabel(han2, 'day', 'fontsize', 15);
    yh2 = get(gca, 'ylabel');
    p2 = get(yh, 'position');
    p2(1) = 1.5*p2(1);
    set(yh2, 'position', p2);
    
    
end