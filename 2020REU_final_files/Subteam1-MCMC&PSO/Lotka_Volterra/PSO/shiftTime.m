%% Helper for plotAll.m
%  Author:      Christina Catlett
%  Date:        July 10, 2020
%  Desc:        Shifts time of data sets that start at relative year 0, not
%               absolute year 1845; takes time column of results vector

function t_shift = shiftTime(t)
% Shift from 0 -> 1845
t_shift = t.VarName1 + 1845;
end