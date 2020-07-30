
%CODE TO CALCULATE MSE VALUES FOR ALL MICE (CAN BE ADAPTED FOR JOINT
%ALGORITHM)

clear all;

file = 'all_values_IMPORTANT.csv';
vals = readtable(file);
vals = vals{:,:};
vals = vals(2:42,:);

for i = 1:11
    i_str = int2str(i);
    mouse_dat_file = strcat('dat', i_str,'.csv');
    mouse_dat = readtable(mouse_dat_file);
    mouse_dat = mouse_dat{:,:};
    mouse_dat = flip(mouse_dat);
    mouse_dat = mouse_dat';
    num_Readings = size(mouse_dat, 2);
    final_params = vals(:,i);
    error = determineError(mouse_dat, final_params, 0.75);
    i
    error / num_Readings
end