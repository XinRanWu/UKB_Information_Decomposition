%% Prepare the fMRI BOLD data
% load the BOLD data
cd('/public/home/zhangjie/DataAnalysis/wxr_UKB')
func = load('ukb_glasser_ts_2_0.mat');

% load the ukb_2_0 info table
cd('/public/home/yangsy/Brain/Brain_project/PhiID/UKB_PhiID_All/result/glasser/table');
ukb_2_0 = readtable('ukb_2_0_glasser_has_struc_func_info.csv');
ukb_2_0_ID = 10 * ukb_2_0.eid + 2;

% get the ID of the fMRI_2_0 
fMRI_2_0_ID = func.fMRI_Glasser_ID_2_0;
all_parts = cellfun(@(x) strsplit(x, '_'), fMRI_2_0_ID, 'UniformOutput', false);
first_parts = cellfun(@(x) x{1}, all_parts, 'UniformOutput', false);
third_parts = cellfun(@(x) x{3}, all_parts, 'UniformOutput', false);
first_parts = str2double(first_parts);
third_parts = str2double(third_parts);
fMRI_2_0_ID = 10 * first_parts + third_parts;
clear all_parts first_parts third_parts

% Find the indices of elements in 'func.fMRI_aparc_eID' that exist in 'ukb_2_0.eid'.
[~, indice] = ismember(fMRI_2_0_ID,ukb_2_0_ID);
clear fMRI_2_0_ID ukb_2_0_ID 
idx = find(indice ~= 0);
clear indice
ukb_2_0_BOLD_eID = func.fMRI_Glasser_eID_2_0(idx);
ukb_2_0_BOLD = func.fMRI_Glasser_TS_2_0(idx);
ukb_2_0_BOLD_ID = func.fMRI_Glasser_ID_2_0(idx);

% Exclude participants with missing values.
empty_indices = cellfun(@isempty, ukb_2_0_BOLD);
ukb_2_0_BOLD_eID = ukb_2_0_BOLD_eID(~empty_indices);
ukb_2_0_BOLD_ID = ukb_2_0_BOLD_ID(~empty_indices);
ukb_2_0_BOLD = ukb_2_0_BOLD(~empty_indices);
ukb_2_0 = ukb_2_0(~empty_indices,:);
clear empty_indices idx 

% Define a checking function to inspect if there exists a 
% row in the matrix that is entirely filled with zeros
checkFunction = @(matrix) any(all(matrix == 0,2));
% Apply the checking function to each cell and record the index array
empty_indices = cellfun(checkFunction, ukb_2_0_BOLD);

ukb_2_0_BOLD_eID = ukb_2_0_BOLD_eID(~empty_indices);
ukb_2_0_BOLD_ID = ukb_2_0_BOLD_ID(~empty_indices);
ukb_2_0_BOLD = ukb_2_0_BOLD(~empty_indices);
ukb_2_0 = ukb_2_0(~empty_indices,:);
clear checkFunction empty_indices func 

% Smoothing BOLD signal to remove noise.
n_subject = size(ukb_2_0,1);
ukb_2_0_BOLD_smooths_1 = cell(20000,1);
% size(ukb_2_0_BOLD,1)
for i = 1: 20000
    disp(i);
    originalData = ukb_2_0_BOLD{i, 1}; 
    smoothedData = smoothdata(originalData, 2, "movmean", 15);
    ukb_2_0_BOLD_smooths_1{i, 1} = smoothedData;
    clear originalData smoothedData; 
end
cd('/public/home/yangsy/Brain/Brain_project/PhiID/UKB_PhiID_All/data/glasser/middle_data');
save('ukb_2_0_BOLD_smooths_1.mat','ukb_2_0_BOLD_smooths_1','-v7.3');
clear ukb_2_0_BOLD_smooths_1

ukb_2_0_BOLD_smooths_2 = cell(n_subject-20000,1);
for i = 20001: size(ukb_2_0_BOLD,1)
    disp(i);
    originalData = ukb_2_0_BOLD{i, 1}; 
    smoothedData = smoothdata(originalData, 2, "movmean", 15);
    ukb_2_0_BOLD_smooths_2{i-20000, 1} = smoothedData;
    clear originalData smoothedData; 
end
save('ukb_2_0_BOLD_smooths_2.mat','ukb_2_0_BOLD_smooths_2','-v7.3');
clear ukb_2_0_BOLD_smooths_2

% Filtering operation on the BOLD signal.
% Remove NaN and Inf values from ukb_2_0_BOLD_smooths.
load('ukb_2_0_BOLD_smooths_1.mat');
load('ukb_2_0_BOLD_smooths_2.mat');
ukb_2_0_BOLD_smooths = [ukb_2_0_BOLD_smooths_1;ukb_2_0_BOLD_smooths_2];
clear ukb_2_0_BOLD_smooths_1 ukb_2_0_BOLD_smooths_2
save('ukb_2_0_BOLD_smooths.mat','ukb_2_0_BOLD_smooths','-v7.3');
delete('ukb_2_0_BOLD_smooths_1.mat');
delete('ukb_2_0_BOLD_smooths_2.mat');

checkFunction = @(matrix) any(isnan(matrix), 'all');
nan_indices = cellfun(checkFunction, ukb_2_0_BOLD_smooths);
ukb_2_0_BOLD_eID = ukb_2_0_BOLD_eID(~nan_indices);
ukb_2_0_BOLD_ID = ukb_2_0_BOLD_ID(~nan_indices);
ukb_2_0_BOLD = ukb_2_0_BOLD(~nan_indices);
ukb_2_0_BOLD_smooths = ukb_2_0_BOLD_smooths(~nan_indices);
ukb_2_0 = ukb_2_0(~nan_indices,:);

checkFunction = @(matrix) any(isinf(matrix), 'all');
inf_indices = cellfun(checkFunction, ukb_2_0_BOLD_smooths);
ukb_2_0_BOLD_eID = ukb_2_0_BOLD_eID(~inf_indices);
ukb_2_0_BOLD_ID = ukb_2_0_BOLD_ID(~inf_indices);
ukb_2_0_BOLD = ukb_2_0_BOLD(~inf_indices);
ukb_2_0_BOLD_smooths = ukb_2_0_BOLD_smooths(~inf_indices);
ukb_2_0 = ukb_2_0(~inf_indices,:);


D = 360;
fs = 1; % Sampling frequency, measured in Hz
fpass = [0.008, 0.09]; % Passband frequency range of the bandpass filter, measured in Hz
order = 5; % Order of the filter
% Set up the bandpass filter
[b, a] = butter(order, fpass / (fs / 2), 'bandpass');
dt = datetime('now', 'TimeZone', 'Asia/Shanghai');
% parpool('local', 18);
% size(ukb_2_0_BOLD_smooths)
n_subject = size(ukb_2_0,1);
ukb_2_0_BOLD_smooths_filter_1 = cell(10000,1);
for i = 1:10000
    disp(i)
    BOLD = ukb_2_0_BOLD_smooths{i,1};
    for row = 1:D
        BOLD(row,:) = filtfilt(b, a, BOLD(row,:));
    end
    ukb_2_0_BOLD_smooths_filter_1{i,1} = BOLD;
end
save('ukb_2_0_BOLD_smooths_filter_1.mat','ukb_2_0_BOLD_smooths_filter_1','-v7.3');
clear ukb_2_0_BOLD_smooths_filter_1

ukb_2_0_BOLD_smooths_filter_2 = cell(10000,1);
for i = 10001:20000
    disp(i)
    BOLD = ukb_2_0_BOLD_smooths{i,1};
    for row = 1:D
        BOLD(row,:) = filtfilt(b, a, BOLD(row,:));
    end
    ukb_2_0_BOLD_smooths_filter_2{i-10000,1} = BOLD;
end
save('ukb_2_0_BOLD_smooths_filter_2.mat','ukb_2_0_BOLD_smooths_filter_2','-v7.3');
clear ukb_2_0_BOLD_smooths_filter_2

ukb_2_0_BOLD_smooths_filter_3 = cell(10000,1);
for i = 20001:30000
    disp(i)
    BOLD = ukb_2_0_BOLD_smooths{i,1};
    for row = 1:D
        BOLD(row,:) = filtfilt(b, a, BOLD(row,:));
    end
    ukb_2_0_BOLD_smooths_filter_3{i-20000,1} = BOLD;
end
save('ukb_2_0_BOLD_smooths_filter_3.mat','ukb_2_0_BOLD_smooths_filter_3','-v7.3');
clear ukb_2_0_BOLD_smooths_filter_3

ukb_2_0_BOLD_smooths_filter_4 = cell(n_subject-30000,1);
for i = 30001:n_subject
    disp(i)
    BOLD = ukb_2_0_BOLD_smooths{i,1};
    for row = 1:D
        BOLD(row,:) = filtfilt(b, a, BOLD(row,:));
    end
    ukb_2_0_BOLD_smooths_filter_4{i-30000,1} = BOLD;
end

clear ukb_2_0_BOLD_smooths
load('ukb_2_0_BOLD_smooths_filter_1.mat');
load('ukb_2_0_BOLD_smooths_filter_2.mat');
ukb_2_0_BOLD_12 = [ukb_2_0_BOLD_smooths_filter_1;ukb_2_0_BOLD_smooths_filter_2];
clear ukb_2_0_BOLD_smooths_filter_1 ukb_2_0_BOLD_smooths_filter_2
load('ukb_2_0_BOLD_smooths_filter_3.mat');
ukb_2_0_BOLD_34 = [ukb_2_0_BOLD_smooths_filter_3;ukb_2_0_BOLD_smooths_filter_4];
clear ukb_2_0_BOLD_smooths_filter_3 ukb_2_0_BOLD_smooths_filter_4

ukb_2_0_BOLD_smooths_filter = [ukb_2_0_BOLD_12;ukb_2_0_BOLD_34];
clear ukb_2_0_BOLD_12 ukb_2_0_BOLD_34

cd('/public/home/yangsy/Brain/Brain_project/PhiID/UKB_PhiID_All/data/glasser');
writetable(ukb_2_0,'ukb_2_0_glasser_basic.csv');
save('ukb_2_0_glasser_BOLD.mat', 'ukb_2_0_BOLD_smooths_filter',...,
    'ukb_2_0', 'ukb_2_0_BOLD_eID','ukb_2_0_BOLD_ID','-v7.3');


