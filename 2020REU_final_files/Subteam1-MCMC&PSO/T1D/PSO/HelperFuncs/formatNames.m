%% Assign an array to values to a copy of a struct object
%  Author:       C. Catlett
%  Date:         June 10, 2020
%  Desc:         Given a structure of params, assign given vals to
%                corresponsding params in new struct object

function newstruct = formatNames(paramStruct, vals)
fn = fieldnames(paramStruct); % Field names
len = numel(fn);              % Number of fields
% For each field, place param in array
for i = 1:len
    newstruct.(fn{i}) = vals(i);
end
end