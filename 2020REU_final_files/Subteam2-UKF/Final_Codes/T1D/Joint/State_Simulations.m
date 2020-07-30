%Script to produce figures of all the various state simulations with and
%without the apoptotic wave


clear all;
file = 'all_values_IMPORTANT_multiple_eta.csv';
vals = readtable(file);
vals = vals{:,2:end};
T = 350;
InitialState = [4.77*10^5; 0; 0; 0; 300; 100; 10; 0; 0; 0; 0; 0];
T1D_ODE_Parameters;
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
        sim_sol = ode15s(@(t, y) ODE_allparams(t, y, f1n, f2n, wave, sim_params), [0 T], [InitialState; sim_params] ,options); %Simulate ODE's up until first time point
        x_sim_wave = zeros(12, 350);
        for j = 1:350
            temp = deval(sim_sol, j)
            x_sim_wave(:,j) = temp(1:12);
        end
        
        
        baseline_params = [5*10^4; 0.4; 0.09; 0.1; 1*10^(-8); 1*10^(-8);
                     f2s*DCtoM*scale_factor*10^(-5); f2s*tDCtoM*scale_factor*10^(-5); 0.0334; 1/60; 2.59*10^5; 90; 864; 
                     1.44; .72; 43.2; 432; sqrt(20000); .194; 3; 0.1; 0.51; 0.1; 10^5; 0.487e-5; 0.487e-5; .1199; 370; 12; 
                     0.01; 10^(-3); 10^(-3); 0.01; aEmday/5000; aEmday/5000; 2.12e5; 0.5; 1; 36; 0.11; 21];
        base_wave_sol = ode15s(@(t, y) ODE_allparams2(t, y, f1n, f2n, wave, baseline_params), [0 T], [InitialState; baseline_params],options); %Simulate ODE's up until first time point
        x_base_wave = zeros(12, 350);
        for j = 1:350
            temp2 = deval(base_wave_sol, j)
            x_base_wave(:,j) = temp2(1:12);
        end
        
        fig = figure(m);
        %sgtitle('Joint UKF');
            
        %subplot(3,4,4);
        subplot(2,2,1);
        p1 = plot(1:T, x_sim_wave(4,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(4,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Apoptotic beta cells', 'fontsize', 15);
        ylabel('cells/ml', 'fontsize', 15);
        
        

        %subplot(3,4,6);
        subplot(2,2,2);
        p1 = plot(1:T, x_sim_wave(7,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(7,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Insulin', 'fontsize', 15);
        ylabel('\muU', 'fontsize', 15);
        
       

        %subplot(3,4,9);
        subplot(2,2,3);
        p1 = plot(1:T, x_sim_wave(10,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(10,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Effector T cells', 'fontsize', 15);
        ylabel('cells/ml', 'fontsize', 15);
        
        

        %subplot(3,4,10);
        subplot(2,2,4);
        p1 = plot(1:T, x_sim_wave(11,:), 'b', 'Linewidth', 2); hold on;
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
        

      
        
        base_nowave_sol = ode15s(@(t, y) ODE_allparams2(t, y, f1n, f2n, 0, baseline_params), [0 T], [InitialState; baseline_params],options); %Simulate ODE's up until first time point
        x_base_nowave = zeros(12, 350);
        for j = 1:350
            temp = deval(base_nowave_sol, j);
            x_base_nowave(:,j) = temp(1:12);
        end
        
        no_wave_sol = ode15s(@(t, y) ODE_allparams(t, y, f1n, f2n, 0, sim_params), [0 350], [InitialState; sim_params], options); %Simulate ODE's up until first time point
        x_sim_nowave = zeros(12, 350);
        for j = 1:350
            temp = deval(no_wave_sol, j);
            x_sim_nowave(:,j) = temp(1:12);
        end
        fig2 = figure(m+11)
        %sgtitle('Baseline Estimates + No Wave');
        
                %subplot(3,4,4);
        subplot(2,2,1);
        p1 = plot(1:T, x_sim_nowave(4,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(4,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Apoptotic beta cells', 'fontsize', 15);
        ylabel('cells/ml', 'fontsize', 15);
        
        

        %subplot(3,4,6);
        subplot(2,2,2);
        p1 = plot(1:T, x_sim_nowave(7,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(7,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Insulin', 'fontsize', 15);
        ylabel('\muU', 'fontsize', 15);
        
       

        %subplot(3,4,9);
        subplot(2,2,4);
        p1 = plot(1:T, x_sim_nowave(10,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(10,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Effector T cells', 'fontsize', 15);
        ylabel('cells/ml', 'fontsize', 15);
        legend('fitted', 'baseline', 'Location', 'southeast');
      
        

        %subplot(3,4,10);
        subplot(2,2,3);
        p1 = plot(1:T, x_sim_nowave(11,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(11,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Regulatory T cells', 'fontsize', 15);
        ylabel('cells\ml', 'fontsize', 15);
        
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        xlabel(han, 'days from birth', 'fontsize', 20);
        
        
        
        
       
        fig3 = figure(m+22)
        %sgtitle('Joint Parameter Estimates + NO Wave');
        
        subplot(3,4,1);
        p1 = plot(1:T, x_sim_wave(1,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(1,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Resting Macrophages', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,2);
        p1 = plot(1:T, x_sim_wave(2,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(2,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Active Macrophages', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,3);
        p1 = plot(1:T, x_sim_wave(3,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(3,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Healthy Beta cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,4);
        p1 = plot(1:T, x_sim_wave(4,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(4,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Apoptotic Beta cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,5);
        p1 = plot(1:T, x_sim_wave(5,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(5,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Necrotic Beta cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,6);
        p1 = plot(1:T, x_sim_wave(7,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(7,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Insulin', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
       
        subplot(3,4,7);
        p1 = plot(1:T, x_sim_wave(8,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(8,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('D Cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,8);
        p1 = plot(1:T, x_sim_wave(9,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(9,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('tD cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,9);
        p1 = plot(1:T, x_sim_wave(10,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(10,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Effector T cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,10);
        p1 = plot(1:T, x_sim_wave(11,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(11,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Regulatory T cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,11);
        p1 = plot(1:T, x_sim_wave(12,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(12,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Memory T cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,12);
        p1 = plot(1:T, x_sim_wave(6,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_wave(6,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Glucose', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        legend('fitted', 'baseline', 'Location', 'southeast');
        
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        xlabel(han, 'days from birth', 'fontsize', 20);
        
        
       fig3 = figure(m+33)
        %sgtitle('Joint Parameter Estimates + NO Wave');
        
        subplot(3,4,1);
        p1 = plot(1:T, x_sim_nowave(1,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(1,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Resting Macrophages', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,2);
        p1 = plot(1:T, x_sim_nowave(2,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(2,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Active Macrophages', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,3);
        p1 = plot(1:T, x_sim_nowave(3,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(3,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Healthy Beta cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,4);
        p1 = plot(1:T, x_sim_nowave(4,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(4,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Apoptotic Beta cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,5);
        p1 = plot(1:T, x_sim_nowave(5,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(5,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Necrotic Beta cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,6);
        p1 = plot(1:T, x_sim_nowave(7,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(7,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Insulin', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
       
        subplot(3,4,7);
        p1 = plot(1:T, x_sim_nowave(8,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(8,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('D Cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,8);
        p1 = plot(1:T, x_sim_nowave(9,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(9,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('tD cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,9);
        p1 = plot(1:T, x_sim_nowave(10,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(10,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Effector T cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,10);
        p1 = plot(1:T, x_sim_nowave(11,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(11,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Regulatory T cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,11);
        p1 = plot(1:T, x_sim_nowave(12,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(12,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Memory T cells', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        
        subplot(3,4,12);
        p1 = plot(1:T, x_sim_nowave(6,:), 'b', 'Linewidth', 2); hold on;
        p2 = plot(1:T, x_base_nowave(6,:), 'k--', 'Linewidth', 2); hold off;
        set(gca, 'YScale', 'log')
        title('Glucose', 'fontsize', 15)
        ylabel('cells\ml', 'fontsize', 15);
        legend('fitted', 'baseline', 'Location', 'southeast');
        
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        xlabel(han, 'days from birth', 'fontsize', 20);
        
        
        
        
    end  
end