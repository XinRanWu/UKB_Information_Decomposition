% This code is going to compute the mean redundancy and synergy of total
% extra and intra Yeo 7 net;
%% Part 1 Prepare Data
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/data');
load('red_syn/ukb_2_0_glasser_redundancy.mat');
load('red_syn/ukb_2_0_glasser_synergy.mat');

n_subject = size(ukb_2_0_eID,1);
glasser_ukb = readtable('atlas/Glasser_info_ukb.csv');
glasser_ukb = glasser_ukb(1:360,:);

%% Part2 Compute total mean redundancy & synergy of extra--intra 7net 
yeo_7 = readtable('atlas/yeo7CommunityNames.txt');
yeo_7 = yeo_7.name;
id_yeo_7 = arrayfun(@(val) find(strcmpi(glasser_ukb.Yeo_7net_name, val)), yeo_7, 'UniformOutput', false);

n_region = 360;
intra_net_idx = zeros(n_region); 
for i = 1:length(id_yeo_7)
    intra_net_idx(id_yeo_7{i},id_yeo_7{i}) = 1;
end

extra_net_idx = ~intra_net_idx;
intra_net_idx = ~extra_net_idx;

extra_net_red = zeros(n_subject,1);
extra_net_syn = zeros(n_subject,1);
intra_net_red = zeros(n_subject,1);
intra_net_syn = zeros(n_subject,1);
for i = 1:n_subject
    disp(i);
    red_mat = redundancy(:,:,i);
    syn_mat = synergy(:,:,i);
    extra_red_mat = red_mat.*extra_net_idx;
    extra_syn_mat = syn_mat.*extra_net_idx;
    intra_red_mat = red_mat.*intra_net_idx;
    intra_syn_mat = syn_mat.*intra_net_idx;
    extra_net_red(i) = mean(extra_red_mat(extra_red_mat~=0),'all');
    extra_net_syn(i) = mean(extra_syn_mat(extra_syn_mat~=0),'all');
    intra_net_red(i) = mean(intra_red_mat(intra_red_mat~=0),'all');
    intra_net_syn(i) = mean(intra_syn_mat(intra_syn_mat~=0),'all');
end

ukb_total_extra_intra_result_table = addvars(ukb_2_0, intra_net_red,intra_net_syn,...,
    extra_net_red,extra_net_syn,'NewVariableNames', {'total_intra_7net_red','total_intra_7net_syn','total_extra_7net_red','total_extra_7net_syn'});
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/extra_intra');
writetable(ukb_total_extra_intra_result_table, 'ukb_total_extra_intra_red_syn_table.csv');

fitlme(ukb_total_extra_intra_result_table,'total_intra_7net_red~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 2.7684e-05
fitlme(ukb_total_extra_intra_result_table,'total_intra_7net_syn~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 8.0857e-09
fitlme(ukb_total_extra_intra_result_table,'total_extra_7net_red~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 0.16667
fitlme(ukb_total_extra_intra_result_table,'total_extra_7net_syn~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 0.38364






