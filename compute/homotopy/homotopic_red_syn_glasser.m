%% Part 1 Prepare Data
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/data');
load('red_syn/ukb_2_0_glasser_redundancy.mat');
load('red_syn/ukb_2_0_glasser_synergy.mat');

n_subject = size(ukb_2_0_eID,1);
glasser_ukb = readtable('atlas/Glasser_info_ukb.csv');
glasser_ukb = glasser_ukb(1:360,:);

%% Part2 compute
roi_name = glasser_ukb(1:180,:).ROI_Name';
roi_name = cellfun(@(x) strrep(x, '-', '_'), roi_name, 'UniformOutput', false);
column_name = cellfun(@(x) ['homo_' x], roi_name, 'UniformOutput', false);
homo_red = zeros(n_subject,n_region/2);
homo_syn = zeros(n_subject,n_region/2);
for i = 1:n_subject
    disp(i);
    red_mat = redundancy(:,:,i);
    syn_mat = synergy(:,:,i);
    homo_red(i,:) = diag(red_mat,n_region/2)';
    homo_syn(i,:) = diag(syn_mat,n_region/2)';
end

homo_red_table = array2table(homo_red,'VariableNames',column_name);
homo_syn_table = array2table(homo_syn,'VariableNames',column_name);
homo_red_table = [ukb_2_0,homo_red_table];
homo_syn_table = [ukb_2_0,homo_syn_table];

cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/homotopy');
writetable(homo_red_table,'ukb_homotopic_redundancy_glasser_table.csv');
writetable(homo_syn_table,'ukb_homotopic_synergy_glasser_table.csv');



