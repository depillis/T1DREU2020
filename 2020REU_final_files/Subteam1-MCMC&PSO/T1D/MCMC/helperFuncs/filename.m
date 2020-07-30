%% Specify file name for saving workspaces
%  Author:       M. Watanabe
%  Date:         July 2020
%  Desc:         Function to create a unique filename for the
%                parameterization routine to save the workspace and
%                variables.
function filename = filename(head, modeltype, wave, disease, dataset)

% set options in file name
if modeltype == 0
    strMod = 'NOD_';
else
    strMod = 'WT_';
end

if wave == 0
    strWave = 'waveOff_';
else
    strWave = 'waveOn_';
end

if disease == 0
    strDis = 'prog_';
elseif disease == 1
    strDis = 'acute_';
else
    strDis = '';
    
end
if dataset == 0
    strDat = 'mathewsetal';
  
elseif dataset == 1
    strDat = 'lietal';
else
    strDat = 'simDat';
end

filename = strcat(head, strDis, strMod, strWave, strDat);


end