%% Loading Data

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/Schaefer/redundancy_synergy_con_mmi.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/Schaefer/redundancy_synergy_con_ccs.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/Schaefer/redundancy_synergy_dis_mmi.mat')

emptyCells = cellfun(@isempty, redundancy_con_mmi);
idx = find(emptyCells);
for i = 1:length(idx)
    redundancy_con_mmi{idx(i)} = nan(100,100);
end
emptyCells = cellfun(@isempty, synergy_con_mmi);
idx = find(emptyCells);
for i = 1:length(idx)
    synergy_con_mmi{idx(i)} = nan(100,100);
end
redundancy_con_mmi = cat(3,redundancy_con_mmi{:});
synergy_con_mmi = cat(3,synergy_con_mmi{:});


emptyCells = cellfun(@isempty, redundancy_dis_mmi);
idx = find(emptyCells);
for i = 1:length(idx)
    redundancy_dis_mmi{idx(i)} = nan(100,100);
end
emptyCells = cellfun(@isempty, synergy_dis_mmi);
idx = find(emptyCells);
for i = 1:length(idx)
    synergy_dis_mmi{idx(i)} = nan(100,100);
end
redundancy_dis_mmi = cat(3,redundancy_dis_mmi{:});
synergy_dis_mmi = cat(3,synergy_dis_mmi{:});


emptyCells = cellfun(@isempty, redundancy_con_ccs);
idx = find(emptyCells);
for i = 1:length(idx)
    redundancy_con_ccs{idx(i)} = nan(100,100);
end
emptyCells = cellfun(@isempty, synergy_con_ccs);
idx = find(emptyCells);
for i = 1:length(idx)
    synergy_con_ccs{idx(i)} = nan(100,100);
end
redundancy_con_ccs = cat(3,redundancy_con_ccs{:});
synergy_con_ccs = cat(3,synergy_con_ccs{:});


load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/rotate_parcellation-master/perm_centroid_info_HCPMMP1.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP')
conte69_32k_fc_gradient1 = load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/BrainSpace-master/BrainSpace-master/shared/template_gradients/conte69_32k_fc_gradient1.csv');
conte69_32k_fc_gradient2 = load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/BrainSpace-master/BrainSpace-master/shared/template_gradients/conte69_32k_fc_gradient2.csv');
glasser_360_conte69 = load('glasser_360_conte69.csv');

Xgradients = netMean1D(conte69_32k_fc_gradient1',glasser_360_conte69');
r = corr(X,Xgradients','Rows','pairwise','Type','spearman')
p_perm = perm_sphere_p(X,Xgradients',perm_id,'spearman');
scatter(X,Xgradients)

%% Individual Repeatablity

load('/public/home/zhangjie/DataAnalysis/wxr_toolbox/toolbox_matlab/cbrewer2-master/cbrewer2/colorbrewer.mat', 'colorbrewer');
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Data.mat', 'ukb_RedSynYeo7n_2_0')

load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/Parcellaitions_files/Scheafer7n_100p_Info.mat');

Yeo7Names = {'Visual';'Somatomotor';'DorsalAttention';'VentralAttention';'Limbic';'Frontoparietal';'Default'};

[RedYeo7n_Scheafer7n100p_con_mmi] = netMeans2D(redundancy_con_mmi,Scheafer7n_100p_Info.NetworkID);
[SynYeo7n_Scheafer7n100p_con_mmi] = netMeans2D(synergy_con_mmi,Scheafer7n_100p_Info.NetworkID);
RedYeo7n_In = wxr_mat2dia3d(RedYeo7n_Scheafer7n100p_con_mmi);
RedYeo7n_Btw = wxr_mat2vec3d(RedYeo7n_Scheafer7n100p_con_mmi);
SynYeo7n_In = wxr_mat2dia3d(SynYeo7n_Scheafer7n100p_con_mmi);
SynYeo7n_Btw = wxr_mat2vec3d(SynYeo7n_Scheafer7n100p_con_mmi);
UKB_RedSyn_Scheafer7n100p_con_mmi_2_0 = [array2table(fMRI_Scheafer_eID),array2table(RedYeo7n_In),array2table(RedYeo7n_Btw),...
    array2table(SynYeo7n_In),array2table(SynYeo7n_Btw)];
UKB_RedSyn_Scheafer7n100p_con_mmi_2_0.Properties.VariableNames{1} = 'eid';

[RedYeo7n_Scheafer7n100p_dis_mmi] = netMeans2D(redundancy_dis_mmi,Scheafer7n_100p_Info.NetworkID);
[SynYeo7n_Scheafer7n100p_dis_mmi] = netMeans2D(synergy_dis_mmi,Scheafer7n_100p_Info.NetworkID);
RedYeo7n_In = wxr_mat2dia3d(RedYeo7n_Scheafer7n100p_dis_mmi);
RedYeo7n_Btw = wxr_mat2vec3d(RedYeo7n_Scheafer7n100p_dis_mmi);
SynYeo7n_In = wxr_mat2dia3d(SynYeo7n_Scheafer7n100p_dis_mmi);
SynYeo7n_Btw = wxr_mat2vec3d(SynYeo7n_Scheafer7n100p_dis_mmi);
UKB_RedSyn_Scheafer7n100p_dis_mmi_2_0 = [array2table(fMRI_Scheafer_eID),array2table(RedYeo7n_In),array2table(RedYeo7n_Btw),...
    array2table(SynYeo7n_In),array2table(SynYeo7n_Btw)];
UKB_RedSyn_Scheafer7n100p_dis_mmi_2_0.Properties.VariableNames{1} = 'eid';

[RedYeo7n_Scheafer7n100p_con_ccs] = netMeans2D(redundancy_con_ccs,Scheafer7n_100p_Info.NetworkID);
[SynYeo7n_Scheafer7n100p_con_ccs] = netMeans2D(synergy_con_ccs,Scheafer7n_100p_Info.NetworkID);
RedYeo7n_In = wxr_mat2dia3d(RedYeo7n_Scheafer7n100p_con_ccs);
RedYeo7n_Btw = wxr_mat2vec3d(RedYeo7n_Scheafer7n100p_con_ccs);
SynYeo7n_In = wxr_mat2dia3d(SynYeo7n_Scheafer7n100p_con_ccs);
SynYeo7n_Btw = wxr_mat2vec3d(SynYeo7n_Scheafer7n100p_con_ccs);
UKB_RedSyn_Scheafer7n100p_con_ccs_2_0 = [array2table(fMRI_Scheafer_eID),array2table(RedYeo7n_In),array2table(RedYeo7n_Btw),...
    array2table(SynYeo7n_In),array2table(SynYeo7n_Btw)];
UKB_RedSyn_Scheafer7n100p_con_ccs_2_0.Properties.VariableNames{1} = 'eid';


ukb_RedSynYeo7n_2_0 = removevars(ukb_RedSynYeo7n_2_0, {'SynYeo7n_Total','RedYeo7n_Total'});

[b2f,f2b] = idFinderNum(brain_data_360p.eid,UKB_RedSyn_Scheafer7n100p_con_mmi_2_0.eid);
[r_ind p] = corr(UKB_RedSyn_Scheafer7n100p_con_mmi_2_0{f2b,[30:57,2:29]}, ukb_RedSynYeo7n_2_0{b2f,[2:29,30:57]});
[r_ind p] = corr(UKB_RedSyn_Scheafer7n100p_dis_mmi_2_0{f2b,[30:57,2:29]}, ukb_RedSynYeo7n_2_0{b2f,[2:29,30:57]});
[r_ind p] = corr(UKB_RedSyn_Scheafer7n100p_con_ccs_2_0{f2b,[30:57,2:29]}, ukb_RedSynYeo7n_2_0{b2f,[2:29,30:57]});



rVec = diag(r_ind);resVec = rVec([1:28]+28);
pVec = diag(p);

figure;
resMat = icatb_vec2mat(resVec(8:28));
resMat(find(eye(7))) = resVec(1:7);
imagesc(resMat,[-1,1]);set(gca,'DataAspectRatio',[1 1 1]);colorbar;colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)))
set(gca, 'XTick', 1:7, 'XTickLabel', Yeo7Names);set(gca, 'YTick', 1:7, 'YTickLabel', Yeo7Names);xtickangle(90)
h = heatmap(resMat);
h.ColorbarVisible = 'on'; % ÏÔÊ¾ÑÕÉ«Ìõ
h.Title = 'Trans-Parcellation Correlation';
h.YDisplayLabels = Yeo7Names;
h.XDisplayLabels = Yeo7Names;
colorbar;colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)))
h.ColorLimits = [-1, 1];

load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/rotate_parcellation-master/perm_centroid_info_fs5.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'Red360p_NC', 'Syn360p_NC')

brain_data_360p = [ukb_RedSynYeo7n_2_0(:,1), array2table(Red360p_NC),array2table(Syn360p_NC)];

X_Scheafer7n100p_con_mmi = nanmean(nanmean(redundancy_con_mmi(:,:,f2b)),3);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,1:360),1);

scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_con_mmi,'schaefer_100_fsa5');
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
scheafer7n_100_fsa5(scheafer7n_100_fsa5==0)=NaN;
glasser_360_fsa5(glasser_360_fsa5==0)=NaN;
% redundancy_aparc_68_fsa5 = parcel_to_surface(1:68,'aparc_fsa5');
% glasser_360_fsa5 = parcel_to_surface(1:360,'glasser_360_fsa5');
r_perm_red = corr(scheafer7n_100_fsa5',glasser_360_fsa5','Rows','pairwise','Type','spearman');
p_perm_red = perm_sphere_p(scheafer7n_100_fsa5',glasser_360_fsa5',perm_id(:,1:1000),'spearman');

X_Scheafer7n100p_con_mmi = nanmean(nanmean(synergy_con_mmi(:,:,f2b)),3);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,361:720),1);

scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_con_mmi,'schaefer_100_fsa5');
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
scheafer7n_100_fsa5(scheafer7n_100_fsa5==0)=NaN;
glasser_360_fsa5(glasser_360_fsa5==0)=NaN;
% redundancy_aparc_68_fsa5 = parcel_to_surface(1:68,'aparc_fsa5');
% glasser_360_fsa5 = parcel_to_surface(1:360,'glasser_360_fsa5');
r_perm_syn = corr(scheafer7n_100_fsa5',glasser_360_fsa5','Rows','pairwise','Type','spearman');
p_perm_syn = perm_sphere_p(scheafer7n_100_fsa5',glasser_360_fsa5',perm_id(:,1:1000),'spearman');





X_Scheafer7n100p_dis_mmi = nanmean(nanmean(redundancy_dis_mmi(:,:,f2b)),3);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,1:360),1);

scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_dis_mmi,'schaefer_100_fsa5');
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
scheafer7n_100_fsa5(scheafer7n_100_fsa5==0)=NaN;
glasser_360_fsa5(glasser_360_fsa5==0)=NaN;
% redundancy_aparc_68_fsa5 = parcel_to_surface(1:68,'aparc_fsa5');
% glasser_360_fsa5 = parcel_to_surface(1:360,'glasser_360_fsa5');
r_perm_red = corr(scheafer7n_100_fsa5',glasser_360_fsa5','Rows','pairwise','Type','spearman');
p_perm_red = perm_sphere_p(scheafer7n_100_fsa5',glasser_360_fsa5',perm_id,'spearman');

X_Scheafer7n100p_dis_mmi = nanmean(nanmean(synergy_dis_mmi(:,:,f2b)),3);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,361:720),1);

scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_dis_mmi,'schaefer_100_fsa5');
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
scheafer7n_100_fsa5(scheafer7n_100_fsa5==0)=NaN;
glasser_360_fsa5(glasser_360_fsa5==0)=NaN;
% redundancy_aparc_68_fsa5 = parcel_to_surface(1:68,'aparc_fsa5');
% glasser_360_fsa5 = parcel_to_surface(1:360,'glasser_360_fsa5');
r_perm_syn = corr(scheafer7n_100_fsa5',glasser_360_fsa5','Rows','pairwise','Type','spearman');
p_perm_syn = perm_sphere_p(scheafer7n_100_fsa5',glasser_360_fsa5',perm_id,'spearman');




X_Scheafer7n100p_con_ccs = nanmean(nanmean(redundancy_con_ccs(:,:,f2b)),3);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,1:360),1);

scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_con_ccs,'schaefer_100_fsa5');
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
scheafer7n_100_fsa5(scheafer7n_100_fsa5==0)=NaN;
glasser_360_fsa5(glasser_360_fsa5==0)=NaN;
% redundancy_aparc_68_fsa5 = parcel_to_surface(1:68,'aparc_fsa5');
% glasser_360_fsa5 = parcel_to_surface(1:360,'glasser_360_fsa5');
r_perm_red = corr(scheafer7n_100_fsa5',glasser_360_fsa5','Rows','pairwise','Type','spearman');
p_perm_red = perm_sphere_p(scheafer7n_100_fsa5',glasser_360_fsa5',perm_id,'spearman');

X_Scheafer7n100p_con_ccs = nanmean(nanmean(synergy_con_ccs(:,:,f2b)),3);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,361:720),1);

scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_con_ccs,'schaefer_100_fsa5');
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
scheafer7n_100_fsa5(scheafer7n_100_fsa5==0)=NaN;
glasser_360_fsa5(glasser_360_fsa5==0)=NaN;
% redundancy_aparc_68_fsa5 = parcel_to_surface(1:68,'aparc_fsa5');
% glasser_360_fsa5 = parcel_to_surface(1:360,'glasser_360_fsa5');
r_perm_syn = corr(scheafer7n_100_fsa5',glasser_360_fsa5','Rows','pairwise','Type','spearman');
p_perm_syn = perm_sphere_p(scheafer7n_100_fsa5',glasser_360_fsa5',perm_id,'spearman');


%% Spatial Correlation (Scatter Maps)

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP');

scatter(scheafer7n_100_fsa5',glasser_360_fsa5')

Yeo7_fsa5 = parcel_to_surface(Yeo7MMP,'glasser_360_fsa5');
Yeo7Color = [
    120 18 134;  % 1. Visual (×ÏÉ«)
    70 130 180;  % 2. Somatomotor (À¶É«)
    0 118 14;    % 3. Dorsal Attention (ÂÌÉ«)
    196 58 250;  % 4. Ventral Attention (Ç³×Ï/·Û)
    220 248 164; % 5. Limbic (Ç³»Æ)
    230 148 34;  % 6. Frontoparietal (³ÈÉ«)
    205 62 78    % 7. Default (ºìÉ«)
    ] ./ 255;
figure;

subplot(2,3,1);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,1:360),1);
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
X_Scheafer7n100p_con_mmi = nanmean(nanmean(redundancy_con_mmi(:,:,f2b)),3);
scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_con_mmi,'schaefer_100_fsa5');
z = find((scheafer7n_100_fsa5.* glasser_360_fsa5)~=0);
scatter(scheafer7n_100_fsa5(z), glasser_360_fsa5(z), 36, Yeo7_fsa5(z), 'filled');
colormap(Yeo7Color);caxis([0.5 7.5]);
xlabel('Schaefer7n100p con-MMI Redundancy (fsa5)');
ylabel('Glasser360p con-MMI Redundancy (fsa5)');
grid on;hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
x_data = scheafer7n_100_fsa5(z);y_data = glasser_360_fsa5(z);
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;

subplot(2,3,2);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,1:360),1);
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
X_Scheafer7n100p_dis_mmi = nanmean(nanmean(redundancy_dis_mmi(:,:,f2b)),3);
scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_dis_mmi,'schaefer_100_fsa5');
z = find((scheafer7n_100_fsa5.* glasser_360_fsa5)~=0);
scatter(scheafer7n_100_fsa5(z), glasser_360_fsa5(z), 36, Yeo7_fsa5(z), 'filled');
colormap(Yeo7Color);caxis([0.5 7.5]);
xlabel('Schaefer7n100p con-MMI Redundancy (fsa5)');
ylabel('Glasser360p dis-MMI Redundancy (fsa5)');
grid on;hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
x_data = scheafer7n_100_fsa5(z);y_data = glasser_360_fsa5(z);
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;

subplot(2,3,3);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,1:360),1);
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
X_Scheafer7n100p_con_ccs = nanmean(nanmean(redundancy_con_ccs(:,:,f2b)),3);
scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_con_ccs,'schaefer_100_fsa5');
z = find((scheafer7n_100_fsa5.* glasser_360_fsa5)~=0);
scatter(scheafer7n_100_fsa5(z), glasser_360_fsa5(z), 36, Yeo7_fsa5(z), 'filled');
colormap(Yeo7Color);caxis([0.5 7.5]);
xlabel('Schaefer7n100p con-MMI Redundancy (fsa5)');
ylabel('Glasser360p con-CCS Redundancy (fsa5)');
grid on;hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
x_data = scheafer7n_100_fsa5(z);y_data = glasser_360_fsa5(z);
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;

subplot(2,3,4);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,[1:360]+360),1);
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
X_Scheafer7n100p_con_mmi = nanmean(nanmean(synergy_con_mmi(:,:,f2b)),3);
scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_con_mmi,'schaefer_100_fsa5');
z = find((scheafer7n_100_fsa5.* glasser_360_fsa5)~=0);
scatter(scheafer7n_100_fsa5(z), glasser_360_fsa5(z), 36, Yeo7_fsa5(z), 'filled');
colormap(Yeo7Color);caxis([0.5 7.5]);
xlabel('Schaefer7n100p Synergy con-MMI (fsa5)');
ylabel('Glasser360p con-MMI Synergy (fsa5)');
grid on;hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
x_data = scheafer7n_100_fsa5(z);y_data = glasser_360_fsa5(z);
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;


subplot(2,3,5);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,[1:360]+360),1);
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
X_Scheafer7n100p_dis_mmi = nanmean(nanmean(synergy_dis_mmi(:,:,f2b)),3);
scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_dis_mmi,'schaefer_100_fsa5');
z = find((scheafer7n_100_fsa5.* glasser_360_fsa5)~=0);
scatter(scheafer7n_100_fsa5(z), glasser_360_fsa5(z), 36, Yeo7_fsa5(z), 'filled');
colormap(Yeo7Color);caxis([0.5 7.5]);
xlabel('Schaefer7n100p con-MMI Synergy (fsa5)');
ylabel('Glasser360p dis-MMI Synergy (fsa5)');
grid on;hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
x_data = scheafer7n_100_fsa5(z);y_data = glasser_360_fsa5(z);
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;


subplot(2,3,6);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,[1:360]+360),1);
glasser_360_fsa5 = parcel_to_surface(X_Glasser,'glasser_360_fsa5');
X_Scheafer7n100p_con_ccs = nanmean(nanmean(synergy_con_ccs(:,:,f2b)),3);
scheafer7n_100_fsa5 = parcel_to_surface(X_Scheafer7n100p_con_ccs,'schaefer_100_fsa5');
z = find((scheafer7n_100_fsa5.* glasser_360_fsa5)~=0);
scatter(scheafer7n_100_fsa5(z), glasser_360_fsa5(z), 36, Yeo7_fsa5(z), 'filled');
colormap(Yeo7Color);caxis([0.5 7.5]);
xlabel('Schaefer7n100p con-MMI Synergy (fsa5)');
ylabel('Glasser360p con-CCS Synergy (fsa5)');
grid on;hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
x_data = scheafer7n_100_fsa5(z);y_data = glasser_360_fsa5(z);
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;

%% Spatial Correlation (Brain Maps)

X_Glasser = zscore(nanmean(brain_data_360p{b2f,2:end}(:,[1:360]),1));
[a, cb] = plot_cortical(parcel_to_surface(X_Glasser([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-3.2,3.2]);

X_Glasser = zscore(nanmean(brain_data_360p{b2f,2:end}(:,[1:360]+360),1));
[a, cb] = plot_cortical(parcel_to_surface(X_Glasser([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-3.2,3.2]);


X_Scheafer7n100p = zscore(nanmean(nanmean(redundancy_con_mmi(:,:,f2b)),3));
[a, cb] = plot_cortical(parcel_to_surface(X_Scheafer7n100p,'schaefer_100_fsa5'),'color_range',[-3.2,3.2]);

X_Scheafer7n100p = zscore(nanmean(nanmean(redundancy_dis_mmi(:,:,f2b)),3));
[a, cb] = plot_cortical(parcel_to_surface(X_Scheafer7n100p,'schaefer_100_fsa5'),'color_range',[-3.2,3.2]);

X_Scheafer7n100p = zscore(nanmean(nanmean(redundancy_con_ccs(:,:,f2b)),3));
[a, cb] = plot_cortical(parcel_to_surface(X_Scheafer7n100p,'schaefer_100_fsa5'),'color_range',[-3.2,3.2]);


X_Scheafer7n100p = zscore(nanmean(nanmean(synergy_con_mmi(:,:,f2b)),3));
[a, cb] = plot_cortical(parcel_to_surface(X_Scheafer7n100p,'schaefer_100_fsa5'),'color_range',[-3.2,3.2]);

X_Scheafer7n100p = zscore(nanmean(nanmean(synergy_dis_mmi(:,:,f2b)),3));
[a, cb] = plot_cortical(parcel_to_surface(X_Scheafer7n100p,'schaefer_100_fsa5'),'color_range',[-3.2,3.2]);

X_Scheafer7n100p = zscore(nanmean(nanmean(synergy_con_ccs(:,:,f2b)),3));
[a, cb] = plot_cortical(parcel_to_surface(X_Scheafer7n100p,'schaefer_100_fsa5'),'color_range',[-3.2,3.2]);



%% Gradients

conte69_32k_fc_gradient1 = load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/BrainSpace-master/BrainSpace-master/shared/template_gradients/conte69_32k_fc_gradient1.csv');
conte69_32k_fc_gradient2 = load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/BrainSpace-master/BrainSpace-master/shared/template_gradients/conte69_32k_fc_gradient2.csv');
glasser_360_conte69 = load('glasser_360_conte69.csv');

load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/rotate_parcellation-master/perm_centroid_info_HCPMMP1.mat')

Xgradients = netMean1D(conte69_32k_fc_gradient1',glasser_360_conte69');

figure;
[a, cb] = plot_cortical(parcel_to_surface(Xgradients([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-8.5,8.5]);


X_Glasser = nanmean(brain_data_360p{:,2:end}(:,[1:360]),1);
r_perm_red = corr(X_Glasser',Xgradients','Rows','pairwise','Type','spearman');
p_perm_red = perm_sphere_p(X_Glasser',Xgradients',perm_id,'spearman');

X_Glasser = nanmean(brain_data_360p{:,2:end}(:,[1:360]+360),1);
r_perm_syn = corr(X_Glasser',Xgradients','Rows','pairwise','Type','spearman');
p_perm_syn = perm_sphere_p(X_Glasser',Xgradients',perm_id,'spearman');

v = unique(results_demo_360p.xname);
for i = 1:length(v)
    X = results_demo_360p.t(contains(results_demo_360p.xname,v{i}));
    X_Glasser = X([1:360])';
    r_perm_red(i,1) = corr(X_Glasser',Xgradients','Rows','pairwise','Type','spearman');
    p_perm_red(i,1) = perm_sphere_p(X_Glasser',Xgradients',perm_id,'spearman');
    
    X_Glasser = X([1:360]+360)';
    r_perm_syn(i,1) = corr(X_Glasser',Xgradients','Rows','pairwise','Type','spearman');
    p_perm_syn(i,1) = perm_sphere_p(X_Glasser',Xgradients',perm_id,'spearman');
    
    disp(['Permuation Test of  ',v{i}])
end



formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];

% 360 regions
UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_Outcomes_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = {'x20016_2_0','x4526_2_0','x2178_2_0'}';
results_pheno_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

v = {'x20016_2_0','x4526_2_0','x2178_2_0'}';
for i = 1:length(v)
    X = results_pheno_360p.t(contains(results_pheno_360p.xname,v{i}));
    X_Glasser = X([1:360])';
    r_perm_red(i,1) = corr(X_Glasser',Xgradients','Rows','pairwise','Type','spearman');
    p_perm_red(i,1) = perm_sphere_p(X_Glasser',Xgradients',perm_id,'spearman');
    
    X_Glasser = X([1:360]+360)';
    r_perm_syn(i,1) = corr(X_Glasser',Xgradients','Rows','pairwise','Type','spearman');
    p_perm_syn(i,1) = perm_sphere_p(X_Glasser',Xgradients',perm_id,'spearman');
    
    disp(['Permuation Test of  ',v{i}])
end





figure;

X = brain_data_360p{b2f,2:end};

subplot(1,2,1);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,[1:360]),1);
scatter(Xgradients ,X_Glasser, 36, Yeo7Color(Yeo7MMP,:), 'filled');
colormap(Yeo7Color);caxis([0.5 7.5]);
xlabel('Sensorimotor-Association Axis');
ylabel('Redundancy');
grid on;hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
z = ~isnan(sum([Xgradients', X_Glasser'],2));
x_data = Xgradients(z)';y_data = X_Glasser(z)';
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;

subplot(1,2,2);
X_Glasser = nanmean(brain_data_360p{b2f,2:end}(:,[1:360]+360),1);
scatter(Xgradients ,X_Glasser, 36, Yeo7Color(Yeo7MMP,:), 'filled');
colormap(Yeo7Color);caxis([0.5 7.5]);
xlabel('Sensorimotor-Association Axis');
ylabel('Synergy');
grid on;hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
z = ~isnan(sum([Xgradients', X_Glasser'],2));
x_data = Xgradients(z)';y_data = X_Glasser(z)';
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;


%% Association analysis (removing low SNR Subj & including more covaraibles)

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Data.mat', 'ukb_RedSynYeo7n_2_0')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'Red360p_NC', 'Syn360p_NC')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'ukb_RedYeo7n_2_0', 'ukb_SynYeo7n_2_0')
load('/public/home/zhangjie/DataAnalysis/wxr_toolbox/toolbox_matlab/cbrewer2-master/cbrewer2/colorbrewer.mat', 'colorbrewer');
UKB_Outcomes_2_0 = readtable('/public/home/zhangjie/UKB_Outcomes_2_0.csv');
UKB_Lifestyle_0_0 = readtable('/public/home/zhangjie/UKB_Lifestyle_MRI_0_0.csv');
UKB_Lifestyle_0_0.x1160_L = double(UKB_Lifestyle_0_0.x1160 >= 9);
UKB_Lifestyle_0_0.x1160_S = double(UKB_Lifestyle_0_0.x1160 < 6);
UKB_Lifestyle_0_0.x1160_S(isnan(UKB_Lifestyle_0_0.x1160)) = NaN;
UKB_Lifestyle_0_0.x1160_L(isnan(UKB_Lifestyle_0_0.x1160)) = NaN;
UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, 'x1160');
x709 = UKB_Lifestyle_0_0.x709;
UKB_Lifestyle_0_0.x709 = double(UKB_Lifestyle_0_0.x709 == 1);
UKB_Lifestyle_0_0.x709(isnan(x709)) = NaN;

UKB_BrainMRICov = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BrainMRICov.txt');
UKB_BrainMRICov = UKB_BrainMRICov(...
    UKB_BrainMRICov.SNR_2_0 <(nanmean(UKB_BrainMRICov.SNR_2_0)+...
    3*nanstd(UKB_BrainMRICov.SNR_2_0)) & ...
    UKB_BrainMRICov.HeadMotion_2_0 < 0.2,:);



%% Segmented Regression

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'Red360p_NC', 'Syn360p_NC')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'ukb_2_0_eID')
UKB_rfMRIcov = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_rfMRIcov.csv');
Redundancy = netMean1D(Red360p_NC,Yeo7MMP);
Synergy = netMean1D(Syn360p_NC,Yeo7MMP);
UKB_PID_Yeo7n = [array2table(ukb_2_0_eID),array2table(Redundancy), array2table(Synergy)];
UKB_PID_Yeo7n.Properties.VariableNames{1} = 'eid';
UKB_PID_Yeo7n = innerjoin(UKB_rfMRIcov(UKB_rfMRIcov.HeadMotion_2_0 < 0.2,1),UKB_PID_Yeo7n);
writetable(UKB_PID_Yeo7n,'/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_PID_Yeo7n_2_0.txt','Delimiter' ,',')



%% Long& Short Range Connection
[short_idx,middle_idx,long_idx] = partition_distance_bins(Glasser360DistMatrix, 181:360, 1:180);
fprintf('\nLoading redundancy data (37 GB)...\n');
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/red_syn/ukb_2_0_glasser_redundancy.mat', 'ukb_2_0_eID', 'redundancy');
num_subs = size(redundancy, 3);
fprintf('Extracting means for %d subjects...\n', num_subs);
red_mean3_short_HY = zeros(num_subs, 1);
red_mean3_mid_HY   = zeros(num_subs, 1);
red_mean3_long_HY  = zeros(num_subs, 1);
for i = 1:num_subs
    sub_mat = redundancy(:, :, i);
    red_mean3_short_HY(i) = mean(sub_mat(short_idx));
    red_mean3_mid_HY(i)   = mean(sub_mat(middle_idx));
    red_mean3_long_HY(i)  = mean(sub_mat(long_idx));
end
clear redundancy;
fprintf('Redundancy processing done. Memory cleared!\n');
redundancy_3bin_HY = [red_mean3_long_HY, red_mean3_mid_HY, red_mean3_short_HY];
fprintf('\nLoading synergy data (37 GB)...\n');
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/red_syn/ukb_2_0_glasser_synergy.mat', 'synergy');
syn_mean3_short_HY = zeros(num_subs, 1);
syn_mean3_mid_HY   = zeros(num_subs, 1);
syn_mean3_long_HY  = zeros(num_subs, 1);
for i = 1:num_subs
    sub_mat = synergy(:, :, i);
    syn_mean3_short_HY(i) = mean(sub_mat(short_idx));
    syn_mean3_mid_HY(i)   = mean(sub_mat(middle_idx));
    syn_mean3_long_HY(i)  = mean(sub_mat(long_idx));
end
clear synergy;
fprintf('Synergy processing done. Memory cleared!\n');
synergy_3bin_HY = [syn_mean3_long_HY, syn_mean3_mid_HY, syn_mean3_short_HY];
brain_data_3bin_HY = [array2table(ukb_2_0_eID),array2table(redundancy_3bin_HY),array2table(synergy_3bin_HY)];
brain_data_3bin_HY.Properties.VariableNames{1} = 'eid';
writetable(brain_data_3bin_HY,'UKB_PID_HY3bin_2_0.txt')

%% Functional Connection
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/fMRI/FC/ukb_glasser_FC_2_0.mat')
num_subs = size(glasser360_FC, 3);
fprintf('Extracting means for %d subjects...\n', num_subs);
FC360p_NC = squeeze(sum(glasser360_FC,2));
brain_data_360p_FC = [array2table(fMRI_Glasser_eID_2_0), array2table(FC360p_NC')];
[FCYeo7n_glasser] = NetMeanMatrix3D(glasser360_FC,Yeo7MMP);
FCYeo7n_In = wxr_mat2dia3d(FCYeo7n_glasser);
FCYeo7n_Btw = wxr_mat2vec3d(FCYeo7n_glasser);
brain_data_7n_matrix_FC = [array2table(fMRI_Glasser_eID_2_0), array2table(FCYeo7n_In),array2table(FCYeo7n_Btw)];
FC7n_NC = netMean1D(brain_data_360p_FC{:,2:end},Yeo7MMP);
brain_data_7n_FC = [array2table(fMRI_Glasser_eID_2_0), array2table(FC7n_NC)];
brain_data_7n_FC.Properties.VariableNames{1} = 'eid';
brain_data_360p_FC.Properties.VariableNames{1} = 'eid';
brain_data_7n_matrix_FC.Properties.VariableNames{1} = 'eid';
clear glasser360_FC;
fprintf('functional connectivity processing done. Memory cleared!\n');

save /public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_FC_Data.mat ...
    brain_data_360p_FC brain_data_7n_FC brain_data_7n_matrix_FC -v7.3

Yeo7Names = {'Visual';'Somatomotor';'DorsalAttention';'VentralAttention';'Limbic';'Frontoparietal';'Default'};
load('/public/home/zhangjie/DataAnalysis/wxr_toolbox/toolbox_matlab/cbrewer2-master/cbrewer2/colorbrewer.mat', 'colorbrewer');
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Results.mat', 'brain_data_360p', 'brain_data_7n', 'brain_data_7n_matrix')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_FC_Data.mat');

[b2f,f2b] = idFinderNum(brain_data_7n.eid,brain_data_7n_FC.eid);
[r_ind p] = corr(brain_data_7n_FC{f2b,2:8}, brain_data_7n{b2f,2:end});


[r_ind p] = corr(brain_data_7n_matrix_FC{f2b,2:29}, brain_data_7n_matrix{b2f,[2:29]});
rVec = diag(r_ind);

[r_ind p] = corr(brain_data_7n_matrix_FC{f2b,2:29}, brain_data_7n_matrix{b2f,[31:58]});
rVec = diag(r_ind);

figure;
resMat = icatb_vec2mat(rVec(8:28));
resMat(find(eye(7))) = rVec(1:7);
imagesc(resMat,[-1,1]);set(gca,'DataAspectRatio',[1 1 1]);colorbar;colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)))
set(gca, 'XTick', 1:7, 'XTickLabel', Yeo7Names);set(gca, 'YTick', 1:7, 'YTickLabel', Yeo7Names);xtickangle(90)
h = heatmap(resMat);
h.ColorbarVisible = 'on'; 
h.Title = 'Trans-Parcellation Correlation';
h.YDisplayLabels = Yeo7Names;
h.XDisplayLabels = Yeo7Names;
colorbar;colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)))
h.ColorLimits = [-1, 1];



load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Results.mat', 'UKB_Basic', 'UKB_BrainMRICov', 'UKB_Lifestyle_0_0', 'UKB_Outcomes_2_0')
UKB_BrainMRICov = removevars(UKB_BrainMRICov, {'Income_0_0','Income_1_0','Income_2_0','Income_3_0','EduAge_0_0','EduAge_1_0','EduAge_2_0','Sex_0_0','BirthYr_0_0','BirthMon_0_0','Centre_0_0','Centre_1_0','Centre_2_0','Centre_3_0','AgeRecruit_0_0','AgeAttend_0_0','AgeAttend_1_0','AgeAttend_2_0','AgeAttend_3_0','BMI_0_0','BMI_1_0','BMI_2_0','BMI_3_0','BirthWeight_0_0','BirthWeight_1_0','BirthWeight_2_0','Handedness_0_0','Handedness_1_0','Handedness_2_0','BloodPress_0_0','BloodPress_0_1','BloodPress_1_0','BloodPress_1_1','BloodPress_2_0','BloodPress_2_1','BloodPress_3_0','BloodPress_3_1','SmokingStatus_0_0','SmokingStatus_1_0','SmokingStatus_2_0','SmokingStatus_3_0','CigarettesNum_0_0','CigarettesNum_1_0','CigarettesNum_2_0','CigarettesNum_3_0','AlcoholStatus_0_0','AlcoholStatus_1_0','AlcoholStatus_2_0','AlcoholStatus_3_0','AlcoholIntakeFreq_0_0','AlcoholIntakeFreq_1_0','AlcoholIntakeFreq_2_0','AlcoholIntakeFreq_3_0','SleepHour_0_0','SleepHour_1_0','SleepHour_2_0','SleepHour_3_0','Insomnia_0_0','Insomnia_1_0','Insomnia_2_0','Insomnia_3_0','EduYr_0_0','EduYr_1_0','EduYr_2_0','Race_1','Race_2','Race_3','Race_4'});

% 360 regions * age
UKB_BrainXY = innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),brain_data_360p_FC);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p_FC.Properties.VariableNames(2:end)'; 

formula = ['Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_360p_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

formula = ['AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_360p_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_360p_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_360p_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

% 7*7 network-modules
UKB_BrainXY = innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),brain_data_7n_matrix_FC);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix_FC.Properties.VariableNames(2:end)'; 

formula = ['Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_7n_matrix_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_7n_matrix_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_7n_matrix_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_7n_matrix_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

% 7 networks
UKB_BrainXY = innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),brain_data_7n_FC);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_FC.Properties.VariableNames(2:end)'; 
formula = ['Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_7n_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_7n_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_7n_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_7n_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


results_demo_360p_FC = [results_age_360p_FC;results_sex_360p_FC;results_bmi_360p_FC;results_tiv_360p_FC];
results_demo_7n_matrix_FC = [results_age_7n_matrix_FC;results_sex_7n_matrix_FC;results_bmi_7n_matrix_FC;results_tiv_7n_matrix_FC];
results_demo_7n_FC = [results_age_7n_FC;results_sex_7n_FC;results_bmi_7n_FC;results_tiv_7n_FC];


formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];

% 7*7 network-modules
UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),UKB_Outcomes_2_0),brain_data_7n_matrix_FC);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix_FC.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
results_pheno_7n_matrix_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

% 7 networks

UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),UKB_Outcomes_2_0),brain_data_7n_FC);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_FC.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
results_pheno_7n_FC = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

results_pheno_7n(isnan(results_pheno_7n.t),:) = [];
results_pheno_7n_matrix_FC_FC(isnan(results_pheno_7n_matrix.t),:) = [];

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/red_syn/ukb_3_0_glasser_red_syn.mat', 'redundancy', 'synergy')

FC = nanmean(glasser360_FC_3_0,3);
Xred = nanmean(redundancy,3);
Xsyn = nanmean(synergy,3);
subplot(1,2,1)
density_scatter(icatb_mat2vec(FC),icatb_mat2vec(Xred))
subplot(1,2,2)
density_scatter(icatb_mat2vec(FC),icatb_mat2vec(Xsyn))


[a, cb] = plot_cortical(parcel_to_surface(results_age_360p_FC.t([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-24,24],'layout','grid');
[a, cb] = plot_cortical(parcel_to_surface(results_demo_360p.t([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-24,24],'layout','grid');
[a, cb] = plot_cortical(parcel_to_surface(results_demo_360p.t(360+[181:360,1:180]),'glasser_360_fsa5'),'color_range',[-24,24],'layout','grid');
X = icatb_vec2mat(results_age_7n_matrix_FC.t(8:28));X(find(eye(7))) = results_age_7n_matrix_FC.t(1:7);
plot_matrix_Yeo7n(X,34)


X = icatb_vec2mat(results_demo_7n_matrix.t(8:28));X(find(eye(7))) = results_demo_7n_matrix.t(1:7);
plot_matrix_Yeo7n(X,34)


X = icatb_vec2mat(results_demo_7n_matrix.t(37:57));X(find(eye(7))) = results_demo_7n_matrix.t(30:36);
plot_matrix_Yeo7n(X,34)


%% Control Covariables

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/sMRI/UKB_Freesurfer_Glasser360_2_0.mat','UKB_Glasser360_thickness_2_0','UKB_Glasser360_area_2_0', 'roiNames')
UKB_Glasser360_thickness_2_0 = removevars(UKB_Glasser360_thickness_2_0, {'L_NA','R_NA'});
UKB_Glasser360_area_2_0 = removevars(UKB_Glasser360_area_2_0, {'L_NA','R_NA'});

UKB_Freesurfer_ASEG = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/sMRI/UKB_Freesurfer_ASEG.csv');
UKB_Freesurfer_ASEG.HBF_R_HBF_3_0 = str2double(UKB_Freesurfer_ASEG.HBF_R_HBF_3_0);

UKB_Freesurfer_APARC = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/sMRI/UKB_Freesurfer_APARC.csv','TreatAsEmpty','NA');
UKB_Freesurfer_APARC.SA_WB_Total_2_0 = UKB_Freesurfer_APARC.SA_L_Total_2_0 + UKB_Freesurfer_APARC.SA_R_Total_2_0;
UKB_Freesurfer_APARC.CT_WB_Total_2_0 = (UKB_Freesurfer_APARC.CT_L_Total_2_0 + UKB_Freesurfer_APARC.CT_R_Total_2_0)./2;
UKB_Freesurfer_APARC.Euler_Num_L_2_0 = 2 - 2 .* UKB_Freesurfer_APARC.HBF_L_HBF_2_0;
UKB_Freesurfer_APARC.Euler_Num_R_2_0 = 2 - 2 .* UKB_Freesurfer_APARC.HBF_R_HBF_2_0;
UKB_Freesurfer_APARC.Euler_Num_WB_2_0 = UKB_Freesurfer_APARC.Euler_Num_R_2_0 + UKB_Freesurfer_APARC.Euler_Num_L_2_0;

UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_BrainMRICov,UKB_Freesurfer_APARC(:,[1,740,741,744])),UKB_Glasser360_thickness_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
UKB_BrainXY(UKB_BrainXY.Euler_Num_WB_2_0<-217,:) = [];



YList = brain_data_7n.Properties.VariableNames(2:end)'; 

formula = ['Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'CT_WB_Total_2_0+SA_WB_Total_2_0+Euler_Num_WB_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'}; 
results_age_7n_rmStr = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_BrainMRICov,UKB_Freesurfer_APARC(:,[1,740,741,744])),UKB_Glasser360_thickness_2_0),brain_data_7n_matrix);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
UKB_BrainXY(UKB_BrainXY.Euler_Num_WB_2_0<-217,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; 

results_age_7n_matrix_rmStr = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_BrainMRICov,UKB_Freesurfer_APARC(:,[1,740,741,744])),UKB_Glasser360_thickness_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
UKB_BrainXY(UKB_BrainXY.Euler_Num_WB_2_0<-217,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; 

results_age_360p_rmStr = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


UKB_BP_2_0 = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BP_2_0.csv');
UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];

YList = brain_data_7n.Properties.VariableNames(2:end)'; 

formula = ['Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'x2443_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'}; 
results_age_7n_rmDB = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_Freesurfer_ASEG),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
UKB_BrainXY(UKB_BrainXY.Euler_Num_WB_2_0<-217,:) = [];









