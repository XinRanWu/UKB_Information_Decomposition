%% Prepare the data
cd('/public/home/zhangjie')
currentDirectory = pwd;
ukb_2 = standardizeMissing(readtable('UKB_Basic.csv'), 'NA'); % Convert the string 'NA' to NaN missing values
% load the fMRI and Structure data to get the data_ID and data_eID
cd('/public/home/zhangjie/DataAnalysis/wxr_UKB');
func = load('ukb_glasser_ts_2_0.mat');
struc = load('ukb_glasser+s1_tract_2_0_streamline_count.mat');

% find the subjects who participated in the 2_0 fMRI data
fMRI_ID = func.fMRI_Glasser_ID_2_0;
clear func
fmri_eid_2 = [];
for i = 1:size(fMRI_ID,1)
    disp(i);
    parts = split(fMRI_ID{i},'_');
    if str2double(parts{3}) == 2
        fmri_eid_2 = [fmri_eid_2;str2double(parts{1})];
    end
end


struc_ID = struc.dMRI_Glasser_S1_ID_2_0;
clear struc
struc_eid_2 = [];
for i = 1:size(struc_ID,1)
    disp(i);
    parts = split(struc_ID{i},'_');
    if str2double(parts{3}) == 2
        struc_eid_2 = [struc_eid_2;str2double(parts{1})];
    end
end

func_struc_eid = struc_eid_2(ismember(struc_eid_2,fmri_eid_2));

% Filter out the list of subjects who participated in collection sessions
ukb_2 = rmmissing(ukb_2,'DataVariables',{'AgeAttend_2_0'});
column_name = {'eid','Sex_0_0','BirthYr_0_0','Centre_0_0',...,
    'AgeAttend_2_0','BMI_2_0','Handedness_0_0'};
ukb_2 = ukb_2(:,column_name);
ukb_2.AgeAttend_2_0 = str2double(ukb_2.AgeAttend_2_0);
ukb_2.BMI_2_0 = str2double(ukb_2.BMI_2_0);

[~, indice] = ismember(ukb_2.eid,func_struc_eid);
struc_func_idx = find(indice ~= 0);
ukb_2_0 = ukb_2(struc_func_idx,:);

cd('/public/home/yangsy/Brain/Brain_project/PhiID/UKB_PhiID_All/result/glasser/table');
writetable(ukb_2_0,'ukb_2_0_glasser_has_struc_func_info.csv');

