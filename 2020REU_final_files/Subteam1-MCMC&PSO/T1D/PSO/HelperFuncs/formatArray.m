%% Formatting param struct to array
%  Author:       C. Catlett
%  Date:         June 10, 2020
%  Desc:         Given a structure of params, convert to length n row
%                vector

function paramArray = formatArray(paramStruct)
fn = fieldnames(paramStruct); % Field names
len = numel(fn);              % Number of fields
% For each field, place param in array
for i = 1:len
    paramArray(i) = paramStruct.(fn{i});
end
end

