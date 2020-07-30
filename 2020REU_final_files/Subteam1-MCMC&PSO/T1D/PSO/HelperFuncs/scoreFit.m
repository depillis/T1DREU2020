%% Scoring final PSO fit (Mathews & Li Data)
%  Author:       C. Catlett
%  Date:         June 25, 2020
%  Desc:         Produce a time and shape score for the final best fit by
%                segmenting data and model at respective points of onset,
%                assinging score based on onset time difference, and shape
%                score based on aligned behavior pre/post onset time
%                (sum-of-squares)

function [timeScore, res] = scoreFit(Tpre, Tpost, Ypre, Ypost, Tdata, Ydata)
% Combine pre/post model predictions together
Tall = [Tpre; Tpost];
Yall = [Ypre; Ypost];

%%%%%%%%%%%%%%%%%%%%%%
%% Compare Onset Times
%%%%%%%%%%%%%%%%%%%%%%
% Onset: 250mgdl for 2 consecutive readings

% Determine onset time of diabetes (MODEL)
for i = 1:(length(Yall)-1)
  if Yall(i) >=250 && Yall(i+1) >=250
      onsettime_mod = Tall(i+1);
      onsettime_mod_ind = i+1;
      break
  elseif i == length(Yall)-1
      onsettime_mod = Tall(i+1);
      onsettime_mod_ind = i+1;
  end
end

% Determine onset time of diabetes (DATA)
for i = 1:(length(Yall)-1)
  if Ydata(i) >=250 && Ydata(i+1) >=250
      onsettime_data = Tdata(i+1);
      onsettime_data_ind = i+1;
      break
  end
end

% Compare onset times: absoulute difference
timeScore = abs(onsettime_mod - onsettime_data);

%%%%%%%%%%%%%%%%%%%%%%
%% Compare Model Shape
%%%%%%%%%%%%%%%%%%%%%%

% Isolate pre/post onset data/model predictions
modRange_pre =  Yall(1:onsettime_mod_ind);
modRange_post = Yall(onsettime_mod_ind+1:end);

datRange_pre =  Ydata(1:onsettime_data_ind);
datRange_post = Ydata(onsettime_data_ind+1:end);

% Align pre-onset data/model (largest common chunk)
if length(modRange_pre) < length(datRange_pre)
    datRange_pre = datRange_pre(onsettime_data_ind-length(modRange_pre)+1:end);
elseif length(modRange_pre) > length(datRange_pre)
    modRange_pre = modRange_pre(onsettime_mod_ind-length(datRange_pre)+1:end);
end

% Align post-onset data/model (largest common chunk)
if length(modRange_post) < length(datRange_post)
    datRange_post = datRange_post(1:length(modRange_post));
elseif length(modRange_post) > length(datRange_post)
    modRange_post = modRange_post(1:length(datRange_post));
end

% Find aligned sum-of-squares difference
res_pre = sum((datRange_pre-modRange_pre).*(datRange_pre-modRange_pre));
res_post = sum((datRange_post-modRange_post).*(datRange_post-modRange_post));

res = res_pre + res_post;
end