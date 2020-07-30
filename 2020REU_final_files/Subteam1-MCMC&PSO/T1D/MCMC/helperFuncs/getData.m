%% Determine Data Set
%  Author:       M. Watanabe
%  Date:         July 2020
%  Desc:         Function to create a data set for parameterization routine
%                Options are Li et al., Mathews et al., and a simulated
%                data set.

function data = getData(dataset, disease)

if dataset == 0 % MATHEWS ET AL DATA
    % add folder where Mathews et al data is stored
    addpath('T1D_data_files/Mathewsetal_data/Mathews_Data_Digitized_with_Headers');
    addpath('T1D_data_files/Mathewsetal_data/onset_estimation_Mathews');

    % Read in MATHEWS et al mean data
    data_all = readmatrix('acute+progressive_seperated_1_raw.csv');

    % progressive glucose readings
    data_prog = data_all(:,3);
    % acute glucose readings
    data_acute = data_all(:, 6);

    % get onset of diabetes time distributions X = progressive, Y = acute
    [onsetdistX, onsetdistY] = onset_time_dist_mathews;

    % SHIFTED DATA
    if disease == 0
    % progressive
    mu=onsetdistX.mu;
    sigma=onsetdistX.sigma;
    onset=round(exp(mu+sigma^2/2)); % convert lognormal to normal onset time
    data = [[onset-24:1:onset]', data_prog]; % shifted data
    distCI=paramci(onsetdistX); % range of mu and sigma

    % PLOT DATA
    figure(1); clf
    plot(data(:,1), data(:,2), '-o')
    hold on
    set(gca, 'FontSize', 15)
    xlabel('days');
    ylabel('glucose (mg/dl)');
    title('Mean glucose measurements for progressive mice');

    elseif disease == 1
    % acute
    mu=onsetdistY.mu;
    sigma=onsetdistY.sigma;
    onset=round(exp(mu+sigma^2/2)); % convert lognormal to normal onset time
    data=[[onset-24:1:onset]', data_acute]; % shifted data
    distCI=paramci(onsetdistY); % range of mu and sigma

    % PLOT DATA
    figure(1); clf
    plot(data(:,1), data(:,2), '-o')
    hold on
    set(gca, 'FontSize', 15)
    xlabel('days');
    ylabel('glucose (mg/dl)');
    title('Mean glucose measurements for acute mice');
    end

elseif dataset == 1 % Load Li et al data
   
    addpath('T1D_data_files/Lietal_data');
    data = readmatrix('avg_lietal_2.csv');
    data = sortrows(data);
    data = [data(:,1)*7 data(:,2)];
    data = data(1:end-1,:);

    % other options for individual mice
    %data = readmatrix('dat6.csv');
    %data = readmatrix('dat3.csv');
    %data = readmatrix('dat11.csv');
    
    figure(1); clf
    plot(data(:,1), data(:,2), '-*','color', 'k', 'LineWidth', 2)
    hold on
    set(gca, 'FontSize', 15)
    xlabel('days from birth');
    ylabel('glucose (mg/dl)');

    
elseif dataset == 2 %  simulated data
    
    % extract initial parameter values
    p = [1.00000000000000e-08, 1.00000000000000e-08, 0.0161290322580645,...
         0.100000000000000,14.142135623731000,1.000000000000000e-06,...
         1.000000000000000e-06,0.010000000000000,21,...
         0.0100, 4.769999999999999e+05, 300, 100, 10]; 
    % Simulate data
    data = makeSimData(modeltype, parVer, wave, p, disease, dataset);
    data = [data(:,1)*7, data(:,2)]; % convert from weeks to days
    
    % plot
    figure(1); clf
    plot(data(:,1), data(:,2), '-o')
    hold on
    set(gca, 'FontSize', 15)
    xlabel('days');
    ylabel('glucose (mg/dl)');
    title('Mean glucose measurements for simulated data')
    
else
    msg = 'No such data set';
    error(msg);
    
end
end