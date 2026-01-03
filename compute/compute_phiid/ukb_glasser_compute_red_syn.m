%% Compute the full PhiID decomposition for participants with 02 sections
% load data
cd('/public/home/yangsy/Brain/Brain_project/PhiID/UKB_PhiID_All/data/glasser');
load('ukb_2_0_glasser_BOLD.mat');
% Compute

BOLD = ukb_2_0_BOLD_smooths_filter;
clear ukb_2_0_BOLD_smooths_filter
n_subject = size(BOLD,1);
D = 360;

dt = datetime('now', 'TimeZone', 'Asia/Shanghai');
parpool('local', 20);
synergy_35001_38106 = zeros(D,D,3106);
redundancy_35001_38106 = zeros(D,D,3106);
for i = 35001:38106
    disp(i);
    % Computes full PhiID decomposition of input data
    synergy_mat = zeros(D,D);
    redundancy_mat = zeros(D,D);
    BOLD_i = BOLD{i,1};
    parfor row = 1:D
        for col = 1:D
            if row ~= col
                atoms = PhiIDFull([BOLD_i(row,:); BOLD_i(col,:)]);
                synergy_mat(row, col) = atoms.sts;
                redundancy_mat(row, col) = atoms.rtr;
            end
        end
    end
    synergy_35001_38106(:,:,i-35000) = synergy_mat;
    redundancy_35001_38106(:,:,i-35000) = redundancy_mat;
    clear synergy_mat redundancy_mat i BOLD_i
end
cd('/public/home/yangsy/Brain/Brain_project/PhiID/UKB_PhiID_All/data/glasser/middle_data');
save('ukb_2_0_glasser_red_syn_35001_38106.mat','synergy_35001_38106','redundancy_35001_38106','-v7.3');

