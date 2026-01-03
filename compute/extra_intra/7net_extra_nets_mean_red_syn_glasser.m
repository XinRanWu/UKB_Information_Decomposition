%% Part 1 Prepare Data
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/data');
load('red_syn/ukb_2_0_glasser_redundancy.mat');
load('red_syn/ukb_2_0_glasser_synergy.mat');

n_subject = size(ukb_2_0_eID,1);
glasser_ukb = readtable('atlas/Glasser_info_ukb.csv');
glasser_ukb = glasser_ukb(1:360,:);


%% Part2 Compute 7nets mean extra redundancy & synergy 
yeo_7 = readtable('atlas/yeo7CommunityNames.txt');
yeo_7 = yeo_7.name;
id_yeo_7 = arrayfun(@(val) find(strcmpi(glasser_ukb.Yeo_7net_name, val)), yeo_7, 'UniformOutput', false);

n_subject = size(redundancy,3);
mean_extra_red_mat = zeros(n_subject,length(yeo_7));
mean_extra_syn_mat = zeros(n_subject,length(yeo_7));
for i = 1:length(yeo_7)
    disp(i);
    mean_extra_red_mat(:,i) = mean_extra_net(redundancy,id_yeo_7{i});
    mean_extra_syn_mat(:,i) = mean_extra_net(synergy,id_yeo_7{i});
end

column_name = yeo_7';
red_column_name = cellfun(@(val) ['extra_' val '_red'],column_name,'UniformOutput',false);
syn_column_name = cellfun(@(val) ['extra_' val '_syn'],column_name,'UniformOutput',false);
mean_extra_red_table = array2table(mean_extra_red_mat,'VariableNames',red_column_name);
mean_extra_syn_table = array2table(mean_extra_syn_mat,'VariableNames',syn_column_name);

result_table = [ukb_2_0,mean_extra_red_table,mean_extra_syn_table];
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/extra_intra');
writetable(result_table,'ukb_7net_extra_net_mean_red_syn_glasser.csv');



