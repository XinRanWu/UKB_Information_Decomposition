%% Part 1 Prepare Data
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/data');
load('red_syn/ukb_2_0_glasser_redundancy.mat');
load('red_syn/ukb_2_0_glasser_synergy.mat');

n_subject = size(ukb_2_0_eID,1);
glasser_ukb = readtable('atlas/Glasser_info_ukb.csv');
glasser_ukb = glasser_ukb(1:360,:);


%% Part2 Compute 7nets mean redundancy & synergy of left and right globe
yeo_7 = readtable('atlas/yeo7CommunityNames.txt');
yeo_7 = yeo_7.name;
id_yeo_7 = arrayfun(@(val) find(strcmpi(glasser_ukb.Yeo_7net_name, val)), yeo_7, 'UniformOutput', false);

cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/left_right');

for i = 1: length(yeo_7)
    disp(i);
    result_table = mean_res_syn_intra_extra_globe(ukb_2_0,id_yeo_7{i},redundancy,synergy)
    writetable(result_table, ['ukb_' yeo_7{i} '_extra_intra_red_syn_table.csv']);
end




