%% Part1 left_right plot
score_item = 'Fluid_intelligence_score_2_0';
% The csv files in the file are arranged in ascending name order
cd('E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age_fdr\score_table_lower\left_right');
% Traverse all subdirectories
subDirs = dir('*.*'); 
columnname = {'TOTAL','VIS','SMN','DAN','SAN','limbic','FPN','DMN'};
rowname = {'R','L','R+L','R-L'};
subDirs(1) = []; % Remove '.' from the current directory
subDirs(1) = []; % Remove '..'from the upper-level directory
subDirs = subDirs([7,8,6,2,5,4,3,1],:);

red_t = zeros(length(columnname),length(rowname));
red_pfdr = zeros(length(columnname),length(rowname));
syn_t = zeros(length(columnname),length(rowname));
syn_pfdr = zeros(length(columnname),length(rowname));

for i = 1:length(subDirs)
    if subDirs(i).isdir
        disp(i);
        % If it is a directory, iterate recursively
        cd(subDirs(i).name); % Go to subdirectory
        csvFiles = dir('*.csv'); % Only the CSV file in the current directory is found
        % Gets the struct containing '_red_'
        redFiles = csvFiles(contains({csvFiles.name}, '_red_'));
        redFiles = redFiles([4,2,3,1],:);
        % Gets the struct containing '_syn_'
        synFiles = csvFiles(contains({csvFiles.name}, '_syn_'));
        synFiles = synFiles([4,2,3,1],:);
        for j = 1:length(synFiles)
            redTable = readtable(redFiles(j).name); 
            redRows = redTable(strcmp(redTable.X, score_item), :);
            red_t(i,j) = table2array(redRows(:, {'TVal'}));
            red_pfdr(i,j) = table2array(redRows(:, {'Pfdr'}));

            synTable = readtable(synFiles(j).name);
            synRows = synTable(strcmp(synTable.X, score_item), :);
            syn_t(i,j) = table2array(synRows(:, {'TVal'}));
            syn_pfdr(i,j) = table2array(synRows(:, {'Pfdr'}));
        end
        cd('..'); % Returns the upper-level directory
    end
end

clear csvFiles i j redRows synRows redTable synTable subDirs redFiles synFiles

red_pfdr(red_pfdr > 0.05) = nan;
id_nan = isnan(red_pfdr);
red_t(id_nan) = nan;

syn_pfdr(syn_pfdr > 0.05) = nan;
id_nan = isnan(syn_pfdr);
syn_t(id_nan) = nan;

figure(1)
subplot(1,2,1)
h=heatmap(red_t);
title(['Redundancy-', strrep(score_item, '_', ' '),' < 64 years t value'])
% Set the color bar to 'turbo'
colormap(turbo);
max_val = max(red_t(:));
min_val = min(red_t(:));
% Calculate the range of the color bar so that 0 is in the center
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = rowname;
h.YDisplayLabels = columnname;

subplot(1,2,2)
h=heatmap(syn_t);
title(['Synergy-', strrep(score_item, '_', ' '),' < 64 years t value'])
% Set the color bar to 'turbo'
colormap(turbo);
max_val = max(syn_t(:));
min_val = min(syn_t(:));
% Calculate the range of the color bar so that 0 is in the center
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = rowname;
h.YDisplayLabels = columnname;


%% Part2 7net extra net plot
score_item = 'Fluid_intelligence_score_2_0';
cd('E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age_fdr\score_table_lower\extra_intra\ukb_7net_extra_net_mean_red_syn_glasser');
csvFiles = dir('*.csv'); % Only the CSV file in the current directory is found
columnname = {'VIS','SMN','DAN','SAN','limbic','FPN','DMN'};
rowname = {'Redundancy','Synergy'};
redFiles = csvFiles(contains({csvFiles.name}, '_red_'));
redFiles = redFiles([7,6,2,5,4,3,1],:);
synFiles = csvFiles(contains({csvFiles.name}, '_syn_'));
synFiles = synFiles([7,6,2,5,4,3,1],:);

mat_t = zeros(length(columnname),length(rowname));
mat_pfdr = zeros(length(columnname),length(rowname));
for i = 1:length(synFiles)
    redTable = readtable(redFiles(i).name); 
    redRows = redTable(strcmp(redTable.X, score_item), :);
    mat_t(i,1) = table2array(redRows(:, {'TVal'}));
    mat_pfdr(i,1) = table2array(redRows(:, {'Pfdr'}));

    synTable = readtable(synFiles(i).name);
    synRows = synTable(strcmp(synTable.X, score_item), :);
    mat_t(i,2) = table2array(synRows(:, {'TVal'}));
    mat_pfdr(i,2) = table2array(synRows(:, {'Pfdr'}));
end
clear csvFiles i redFiles redRows redTable synFiles synRows synTable

mat_pfdr(mat_pfdr > 0.05) = nan;
id_nan = isnan(mat_pfdr);
mat_t(id_nan) = nan;

figure(1)
h=heatmap(mat_t');
title([strrep(score_item, '_', ' '),' Extra yeo7nets < 64 years t value'])
% Set the color bar to 'turbo'
colormap(turbo);
max_val = max(mat_t(:));
min_val = min(mat_t(:));
% Calculate the range of the color bar so that 0 is in the center
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = columnname;
h.YDisplayLabels = rowname;


%% Part3 net matrix plot
score_item = 'Fluid_intelligence_score_2_0';

% Redundancy
cd('E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age_fdr\score_table_lower\extra_intra\ukb_mean_net_redundancy_matrix_table');
csvFiles = dir('*.csv');% Only the CSV file in the current directory is found
% ** IMPORTANT: The csv files in the file are arranged in deascending name order
[~,index] = sortrows({csvFiles.name}.'); csvFiles = csvFiles(index(end:-1:1)); clear index 
% columnname = {'VIS','SMN','DAN','SAN','limbic','FPN','DMN'};

red_array_t = zeros(length(csvFiles),1);
red_array_pfdr = zeros(length(csvFiles),1);
for i = 1:length(csvFiles)
    disp(i);
    redTable = readtable(csvFiles(i).name);
    redRows = redTable(strcmp(redTable.X, score_item), :);
    red_array_t(i) = table2array(redRows(:, {'TVal'}));
    red_array_pfdr(i) = table2array(redRows(:, {'Pfdr'}));
end
clear i redRows redTable

red_array_t(red_array_pfdr > 0.05) = nan;
red_net = array2network(red_array_t,1);
OrderId = [1,2,6,3,4,5,7];
red_net_t = red_net(OrderId,OrderId);
clear csvFiles red_array_pfdr red_net red_array_t

% Synergy
cd('E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age_fdr\score_table_lower\extra_intra\ukb_mean_net_synergy_matrix_table');
csvFiles = dir('*.csv');% Only the CSV file in the current directory is found
% ** IMPORTANT: The csv files in the file are arranged in deascending name order
[~,index] = sortrows({csvFiles.name}.'); csvFiles = csvFiles(index(end:-1:1)); clear index 
columnname = {'VIS','SMN','DAN','SAN','limbic','FPN','DMN'};

syn_array_t = zeros(length(csvFiles),1);
syn_array_pfdr = zeros(length(csvFiles),1);
for i = 1:length(csvFiles)
    disp(i);
    synTable = readtable(csvFiles(i).name);
    synRows = synTable(strcmp(synTable.X, score_item), :);
    syn_array_t(i) = table2array(synRows(:, {'TVal'}));
    syn_array_pfdr(i) = table2array(synRows(:, {'Pfdr'}));
end
clear i synRows synTable

syn_array_t(syn_array_pfdr > 0.05) = nan;
syn_net = array2network(syn_array_t,1);
OrderId = [1,2,6,3,4,5,7];
syn_net_t = syn_net(OrderId,OrderId);
clear csvFiles syn_array_pfdr syn_net syn_array_t

% plot net
figure(1)
subplot(1,2,1)
h=heatmap(red_net_t);
title([strrep(score_item, '_', ' '),' Redundancy < 64 years t value'])
colormap(turbo);
max_val = max(red_net_t(:));
min_val = min(red_net_t(:));
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = columnname;
h.YDisplayLabels = columnname;

subplot(1,2,2)
h=heatmap(syn_net_t);
title([strrep(score_item, '_', ' '),' Synergy < 64 years t value'])
colormap(turbo);
max_val = max(syn_net_t(:));
min_val = min(syn_net_t(:));
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = columnname;
h.YDisplayLabels = columnname;


