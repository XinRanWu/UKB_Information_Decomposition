% This code is going to compute the mean redundancy and synergy of total
% left and right globe;
%% Part 1 Prepare Data
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/data');
load('red_syn/ukb_2_0_glasser_redundancy.mat');
load('red_syn/ukb_2_0_glasser_synergy.mat');

n_subject = size(ukb_2_0_eID,1);
glasser_ukb = readtable('atlas/Glasser_info_ukb.csv');
glasser_ukb = glasser_ukb(1:360,:);

%% Part2 Compute total mean redundancy & synergy of left and right globe
right_left = {'R';'L'};
id_right_left = arrayfun(@(val) find(strcmpi(glasser_ukb.Hemisphere, val)), right_left, 'UniformOutput', false);

n_region = 360;
intra_right_idx = zeros(n_region); 
intra_right_idx(id_right_left{1},id_right_left{1}) = 1;

intra_left_idx = zeros(n_region); 
intra_left_idx(id_right_left{2},id_right_left{2}) = 1;

intra_left_right_idx = zeros(n_region);
intra_left_right_idx(id_right_left{1},id_right_left{1}) = 1;
intra_left_right_idx(id_right_left{2},id_right_left{2}) = 1;

extra_left_right_idx = zeros(n_region);
extra_left_right_idx(id_right_left{1},id_right_left{2}) = 1;

intra_right_red = zeros(n_subject,1);
intra_right_syn = zeros(n_subject,1);
intra_left_red = zeros(n_subject,1);
intra_left_syn = zeros(n_subject,1);
intra_left_right_red = zeros(n_subject,1);
intra_left_right_syn = zeros(n_subject,1);
extra_left_right_red = zeros(n_subject,1);
extra_left_right_syn = zeros(n_subject,1);
for i = 1:n_subject
    disp(i);
    red_mat = redundancy(:,:,i);
    syn_mat = synergy(:,:,i);
    
    intra_right_red_mat = red_mat.*intra_right_idx;
    intra_right_syn_mat = syn_mat.*intra_right_idx;
    
    intra_left_red_mat = red_mat.*intra_left_idx;
    intra_left_syn_mat = syn_mat.*intra_left_idx;
    
    intra_left_right_red_mat = red_mat.*intra_left_right_idx;
    intra_left_right_syn_mat = syn_mat.*intra_left_right_idx;
    
    extra_left_right_red_mat = red_mat.*extra_left_right_idx;
    extra_left_right_syn_mat = syn_mat.*extra_left_right_idx;

    intra_right_red(i) = mean(intra_right_red_mat(intra_right_red_mat~=0),'all');
    intra_right_syn(i) = mean(intra_right_syn_mat(intra_right_syn_mat~=0),'all');
    intra_left_red(i) = mean(intra_left_red_mat(intra_left_red_mat~=0),'all');
    intra_left_syn(i) = mean(intra_left_syn_mat(intra_left_syn_mat~=0),'all');
    intra_left_right_red(i) = mean(intra_left_right_red_mat(intra_left_right_red_mat~=0),'all');
    intra_left_right_syn(i) = mean(intra_left_right_syn_mat(intra_left_right_syn_mat~=0),'all');
    extra_left_right_red(i) = mean(extra_left_right_red_mat(extra_left_right_red_mat~=0),'all');
    extra_left_right_syn(i) = mean(extra_left_right_syn_mat(extra_left_right_syn_mat~=0),'all');
end

ukb_total_left_right_result_table = addvars(ukb_2_0, intra_right_red,intra_right_syn,...,
    intra_left_red,intra_left_syn,intra_left_right_red,intra_left_right_syn,extra_left_right_red,...,
    extra_left_right_syn,'NewVariableNames', {'total_intra_right_red','total_intra_right_syn',...,
    'total_intra_left_red','total_intra_left_syn','total_intra_left_right_red','total_intra_left_right_syn',...,
    'total_extra_left_right_red','total_extra_left_right_syn'});

cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/left_right');
writetable(ukb_total_left_right_result_table, 'ukb_total_left_right_result_table.csv');

fitlme(ukb_total_left_right_result_table,'total_intra_right_red~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 0.76567
fitlme(ukb_total_left_right_result_table,'total_intra_right_syn~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 1.4266
fitlme(ukb_total_left_right_result_table,'total_intra_left_red~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = -0.15878
fitlme(ukb_total_left_right_result_table,'total_intra_left_syn~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 1.8065 
fitlme(ukb_total_left_right_result_table,'total_intra_left_right_red~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 0.30762
fitlme(ukb_total_left_right_result_table,'total_intra_left_right_syn~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 1.6344
fitlme(ukb_total_left_right_result_table,'total_extra_left_right_red~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 0.26399
fitlme(ukb_total_left_right_result_table,'total_extra_left_right_syn~AgeAttend_2_0+Sex_0_0+BMI_2_0+(1|Centre_0_0)') % p = 1.5771



