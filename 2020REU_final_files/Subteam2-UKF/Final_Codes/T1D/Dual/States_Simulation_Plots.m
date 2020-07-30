%Script to produce figures of all the various state simulations with and
%without the apoptotic wave


clear all;
file = 'all_values_IMPORTANT.csv';
vals = readtable(file);
vals = vals{:,:};
T = 350;
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
fig_2_3_4_5_Parameters;
wave = wave_basal;

count = 1;
for m = 1:11
    if m == 6
        m_str = int2str(m);
        mouse_dat_file = strcat('dat', m_str,'.csv');
        mouse_dat = readtable(mouse_dat_file);
        mouse_dat = mouse_dat{:,:};
        mouse_dat = flip(mouse_dat);
        mouse_dat = mouse_dat';
        mouse_dat(1,:) = mouse_dat(1,:) * 7;
        num_glucose = size(mouse_dat, 2);
        mouse_params = vals(:,m);
        sim_params = mouse_params(2:42);
        options=odeset('RelTol',1e-12,'AbsTol',1e-12);
        sim_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, sim_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
        x_sim_wave = zeros(12, 350);
        for j = 1:350
            x_sim_wave(:,j) = deval(sim_sol, j);
        end
        
        
        baseline_params = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8);
                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21; 0.018];
        base_wave_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, baseline_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
        x_base_wave = zeros(12, 350);
        for j = 1:350
            x_base_wave(:,j) = deval(base_wave_sol, j);
        end
        
        fig = figure(2);
        sgtitle('Dual UKF');
            
        %subplot(3,4,4);
        subplot(2,2,1);
        p1 = plot(1:T, x_sim_wave(4,:), 'r', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(4,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Apoptotic beta cells', 'fontsize', 15);
        ylabel('cells/ml', 'fontsize', 15);
        
        

        %subplot(3,4,6);
        subplot(2,2,2);
        p1 = plot(1:T, x_sim_wave(7,:), 'r', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(7,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Insulin', 'fontsize', 15);
        ylabel('\muU', 'fontsize', 15);
        
       

        %subplot(3,4,9);
        subplot(2,2,3);
        p1 = plot(1:T, x_sim_wave(10,:), 'r', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(10,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Effector T cells', 'fontsize', 15);
        ylabel('cells/ml', 'fontsize', 15);
        
        

        %subplot(3,4,10);
        subplot(2,2,4);
        p1 = plot(1:T, x_sim_wave(11,:), 'r', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(11,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Regulatory T cells', 'fontsize', 15);
        ylabel('cells\ml', 'fontsize', 15);
        legend('fitted', 'baseline', 'Location', 'southeast');
      
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        xlabel(han, 'days from birth', 'fontsize', 20);
        
    
        

      
        
        %{
        baseline_params = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8);
                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21; 0.01];

        base_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, 0, baseline_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
        x_base_wave = zeros(12, 350);
        for j = 1:350
            x_base_wave(:,j) = deval(base_sol, j);
        end
        figure(15)
        sgtitle('Baseline Estimates + No Wave');
        %{
        subplot(3,4,1);
        plot(1:T, x_base_wave(1,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Resting Macrophages', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,2);
        plot(1:T, x_base_wave(2,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Active Mcrophages', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,3);
        plot(1:T, x_base_wave(3,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Healthy Beta Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,4);
        plot(1:T, x_base_wave(4,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Apoptotic Beta Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,5);
        plot(1:T, x_base_wave(5,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Necrotic Beta Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,6);
        plot(1:T, x_base_wave(7,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Insulin', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,7);
        plot(1:T, x_base_wave(8,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('D Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,8);
        plot(1:T, x_base_wave(9,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('tD Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,9);
        plot(1:T, x_base_wave(10,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Effector T Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,10);
        plot(1:T, x_base_wave(11,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Regulatory T Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}

        %{
        subplot(3,4,11);
        plot(1:T, x_base_wave(12,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Memory T Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        %}
        
        %subplot(3,4,12);
        plot(1:T, x_base_wave(6,:), 'r', 'Linewidth', 2);
        ylabel('Glucose', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        
        
        
        
        no_wave_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, 0, sim_params), [0 350], InitialState,options); %Simulate ODE's up until first time point
        x_sim_no_wave = zeros(12, 350);
        for j = 1:350
            x_sim_no_wave(:,j) = deval(no_wave_sol, j);
        end


        figure(3)
        sgtitle('Dual Parameter Estimates + NO Wave');
        %{
        subplot(3,4,1);
        plot(1:T, x_sim_no_wave(1,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Resting Macrophages', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,2);
        plot(1:T, x_sim_no_wave(2,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Active Mcrophages', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,3);
        plot(1:T, x_sim_no_wave(3,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Healthy Beta Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,4);
        plot(1:T, x_sim_no_wave(4,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Apoptotic Beta Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,5);
        plot(1:T, x_sim_no_wave(5,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Necrotic Beta Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,6);
        plot(1:T, x_sim_no_wave(7,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Insulin', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,7);
        plot(1:T, x_sim_no_wave(8,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('D Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,8);
        plot(1:T, x_sim_no_wave(9,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('tD Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,9);
        plot(1:T, x_sim_no_wave(10,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Effector T Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,10);
        plot(1:T, x_sim_no_wave(11,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Regulatory T Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,11);
        plot(1:T, x_sim_no_wave(12,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Memory T Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(3,4,12);
        %}
        plot(1:T, x_sim_no_wave(6,:), 'r', 'Linewidth', 2);
        ylabel('Glucose', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);
        
        
        baseline_params = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8);
                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21; 0.018];
        base_wave_sol = ode15s(@(t, y) fig_2_3_4_5_ODE(t, y, f1n, f2n, wave, baseline_params), [0 T], InitialState,options); %Simulate ODE's up until first time point
        x_base_wave = zeros(12, 350);
        for j = 1:350
            x_base_wave(:,j) = deval(base_wave_sol, j);
        end
        figure(20)
        sgtitle('Baseline Parameters + Wave');
        subplot(2,2,1);
        plot(1:T, x_base_wave(4,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Apoptotic Beta Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(2,2,2);
        plot(1:T, x_base_wave(7,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Insulin', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(2,2,3);
        plot(1:T, x_base_wave(10,:), 'r', 'Linewidth', 2); 
        set(gca, 'YScale', 'log')
        ylabel('Effector T Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        subplot(2,2,4);
        plot(1:T, x_base_wave(11,:), 'r', 'Linewidth', 2);
        set(gca, 'YScale', 'log')
        ylabel('Regulatory T Cells', 'fontsize', 15);
        xlabel('Day', 'fontsize', 15);

        
        %}
        
        
        
        
    end  
end