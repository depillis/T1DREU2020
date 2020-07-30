%% Create simulated mouse glucose data
%  Author:       M. Watanabe
%  Date:         June 2020
%  Desc:         Create simulated mouse glucose data by evaluating the T1D
%                ODE and adding 5% Guassian noise.


function simData = makeSimData(modeltype, parVer, waveOn, params, disease, dataset)
% Load data
if dataset == 0
    % MATHEWS ET AL DATA
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

    elseif disease == 1
    % acute
    mu=onsetdistY.mu;
    sigma=onsetdistY.sigma;
    onset=round(exp(mu+sigma^2/2)); % convert lognormal to normal onset time
    data=[[onset-24:1:onset]', data_acute]; % shifted data
    distCI=paramci(onsetdistY); % range of mu and sigma
    end

elseif dataset == 1 || dataset == 2
    % Load Li et al data
    addpath('T1D_data_files/Lietal_data');
    data = readmatrix('avg_lietal_2.csv');
    %data = readmatrix('dat6.csv');
    data = sortrows(data);
%     data_all = lietal_meandat(disease);
%     [~, onsetdist] = onset_dist_Li(disease);
%     
%     mu=onsetdist.mu;
%     sigma=onsetdist.sigma;
%     onset=round(exp(mu+sigma^2/2)); % convert lognormal to normal onset time
%     data=[7*[onset-22:1:onset]', data_all']; % shifted data
%     distCI=paramci(onsetdist); % range of mu and sigma
    
else
    msg = 'No such data set';
    error(msg);
    
end
% ODE solver
[tMod, yMod] = T1Dfun_r(modeltype, parVer, waveOn, params, data); % only produces glucose data

% Determine noise (10%  Gaussian)
noiseSigma = .05*yMod;
noise = noiseSigma' .* randn(1, length(tMod));
% Add noise to model predictions
wNoise = yMod + noise';
% Create, save simulated data
simData = [tMod wNoise];
simData = sortrows(simData);
save('makesimData.mat', 'simData');
end