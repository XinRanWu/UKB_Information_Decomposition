%% Part 1 Prepare Data
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/data');
load('red_syn/ukb_2_0_glasser_redundancy.mat');
load('red_syn/ukb_2_0_glasser_synergy.mat');

n_subject = size(ukb_2_0_eID,1);
glasser_ukb = readtable('atlas/Glasser_info_ukb.csv');
glasser_ukb = glasser_ukb(1:360,:);

%% Part2 Compute mean redundancy & synergy 7nets matrix
yeo_7 = readtable('atlas/yeo7CommunityNames.txt');
yeo_7 = yeo_7.name;
id_yeo_7 = arrayfun(@(val) find(strcmpi(glasser_ukb.Yeo_7net_name, val)), yeo_7, 'UniformOutput', false);

cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/left_right');

matrix_name = {};
for i = 1:length(yeo_7)
    for j = i:length(yeo_7)
        matrix_name{end+1} = [yeo_7{i} '_' yeo_7{j}];
    end
end

mean_red_mat = zeros(n_subject,length(matrix_name));
mean_syn_mat = zeros(n_subject,length(matrix_name));
for i = 1:n_subject
    disp(i);
    num = 0;
    for j = 1:length(yeo_7)
        for k = j:length(yeo_7)
            num = num +1;
            red_mat = redundancy(id_yeo_7{j},id_yeo_7{k},i);
            syn_mat = synergy(id_yeo_7{j},id_yeo_7{k},i);
            mean_red_mat(i,num) = mean(red_mat,'all');
            mean_syn_mat(i,num) = mean(syn_mat,'all');
        end
    end
    clear red_mat syn_mat
end

mean_red_table = array2table(mean_red_mat, 'VariableNames', matrix_name);
mean_syn_table = array2table(mean_syn_mat, 'VariableNames', matrix_name);

mean_red_matrix_table = [ukb_2_0,mean_red_table];
mean_syn_matrix_table = [ukb_2_0,mean_syn_table];

cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/extra_intra');
writetable(mean_red_matrix_table,'ukb_mean_net_redundancy_matrix_table.csv');
writetable(mean_syn_matrix_table,'ukb_mean_net_synergy_matrix_table.csv');


