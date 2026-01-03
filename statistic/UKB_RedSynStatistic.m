ukb_RedSyn = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/result/yangsy/ukb_aparc_redsyn_2_0.csv');

UKB_BrainXY = innerjoin(UKB_UsedData_Brain,ukb_RedSyn(:,[1,8:end]));
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = [UKB_AllVars.Properties.VariableNames(3:end)';...
    UKB_Biomarker.Properties.VariableNames(2:end)'];

XList = XList(contains(XList,'_2_'));
XList = XList(~(contains(XList,'_2_0_1')|contains(XList,'Intensity_')|...
    contains(XList,'VolRatio_')|contains(XList,'HBF_')|contains(XList,'Vol_WB_TIV_2_0')));
XList = XList(~(contains(XList,'_0_0')|contains(XList,'_1_0')|contains(XList,'_3_0')));

XListAll = [{'TDI_0_0';'Income_2_0';'EduAge_0_0';'Handedness_0_0'};...
    XList;BrainTotalList;BrainAsegList;PRSList;MentalList(sum(UKB_UsedData_Brain{:,2195:2210},1)>200)];
YList = ukb_RedSyn.Properties.VariableNames(8:end)'; % Y.Properties.RedSyniableNames(2:end)';

[RedSynT,RedSynP,RedSynDF,RedSynBeta,RedSynCI1,RedSynCI2] = LMM_withFormula2(UKB_BrainXY,formula,XListAll,YList);

for i = 1:12
    RedSynLMMRes = [table(XListAll),array2table([RedSynBeta(:,i),RedSynT(:,i),RedSynDF(:,i),RedSynP(:,i)])];
    RedSynLMMRes.Properties.VariableNames(1:5) = {'X','Beta','TVal','DF','PVal'};
    writetable(RedSynLMMRes,[YList{i},'_forManhattan.csv']);
end

formula = ['Sex_0_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
YList = ukb_RedSyn.Properties.VariableNames(8:end)';
[RedSynAgeT,RedSynAgeP,RedSynAgeDF,RedSynAgeBeta] = LMM_withFormula2(UKB_BrainXY,formula,{'AgeAttend_2_0'},YList);





XRedSyn = [array2table(ukb_RedSyn.eid),ukb_RedSyn(:,[1,8:16])];
XRedSyn.Properties.VariableNames{2} = 'IID';XRedSyn.Properties.VariableNames{1} = 'FID';

writetable(XRedSyn,'/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/redsyn/UKB_redsyn_pheno.txt','Delimiter' ,' ')

cd('/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/redsyn')
phenofile = '/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/redsyn/UKB_redsyn_pheno.txt';

covarfile = '/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/UKB_GeneCov.txt';
covarname = 'sex,AgeAttend_2_0,gene_pc1,gene_pc2,gene_pc3,gene_pc4,gene_pc5,gene_pc6,gene_pc7,gene_pc8,gene_pc9,gene_pc10';

parfor i = 1:8
    for chr = 1:22
        phenoname = PList{i};
        bfile_qc = ['/home1/ISTBI_data/K_J_J/BackupData/UKBimputeV3_wbs/UKB_gene_v3_imp_qc_chr',...
            num2str(chr)];
        outfile = ['gwas_',phenoname,'_chr',num2str(chr)];
        plink_command = ...
            ['/public/home/zhangjie/DataAnalysis/wxr_toolbox/toobox_other/plink2_linux_x86_64_20191030/plink2 ',...
            '--bfile ',bfile_qc,' ',...
            '--linear  hide-covar --debug ',...
            '--pheno ',phenofile,' ',...
            '--pheno-name ',phenoname,' '...
            '--covar ',covarfile,' ',...
            '--covar-name ',covarname,' ',...
            '--out ',outfile];
        unix(plink_command);
        disp(['GWAS of ',phenoname,' in chr',num2str(chr)])
    end
end

for i = 1:8
    
    unix(['head -n 1 gwas_',PList{i},'_chr1.',PList{i},'.glm.linear > gwas.',PList{i},'.glm.linear']);
    unix(['tail -n +2 gwas_',PList{i},'_chr1.',PList{i},'.glm.linear >> gwas.',PList{i},'.glm.linear']);
    
    for chr = 2:22
        unix(['tail -n +2 gwas_',PList{i},'_chr',num2str(chr),'.',PList{i},'.glm.linear >> gwas.',PList{i},'.glm.linear']);
    end
    
    unix(['sed -i ','''','1s/.*/CHR POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P/',...
        '''',' gwas.',PList{i},'.glm.linear']);
    
    unix(['gzip gwas.',PList{i},'.glm.linear'])
end



load('/public/home/zhangjie/ZJLab/UKBiobank_Project/code/Redundancy_Synergy_Information/UKB_PhiID_All/data/ukb_2_0_dk_syn_red.mat')
DK_label_ukb = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/code/Redundancy_Synergy_Information/UKB_PhiID_All/data/DK_label_ukb.csv');
DK_Yeo7 = nan(68,1);
DK_Yeo7(strcmp(DK_label_ukb.yeo_7, 'VIS')) = 1;
DK_Yeo7(strcmp(DK_label_ukb.yeo_7, 'SOM')) = 2;
DK_Yeo7(strcmp(DK_label_ukb.yeo_7, 'DAN')) = 3;
DK_Yeo7(strcmp(DK_label_ukb.yeo_7, 'VAN')) = 4;
DK_Yeo7(strcmp(DK_label_ukb.yeo_7, 'LIM')) = 5;
DK_Yeo7(strcmp(DK_label_ukb.yeo_7, 'FPN')) = 6;
DK_Yeo7(strcmp(DK_label_ukb.yeo_7, 'DMN')) = 7;




load('/public/home/zhangjie/ZJLab/UKBiobank_Project/code/Redundancy_Synergy_Information/UKB_PhiID_All/data/ukb_2_0_dk_syn_red.mat')

[RedYeo7n] = netMeans2D(redundancy,DK_Yeo7);
[SynYeo7n] = netMeans2D(synergy,DK_Yeo7);
RedYeo7n_In = wxr_mat2dia3d(RedYeo7n);
RedYeo7n_Btw = nanmean(wxr_mat2vec3d(RedYeo7n),2);
SynYeo7n_In = wxr_mat2dia3d(SynYeo7n);
SynYeo7n_Btw = nanmean(wxr_mat2vec3d(SynYeo7n),2);

ukb_RedSyn = removevars(ukb_RedSyn, {'mean_redundancy_without_smn','mean_redundancy_within_smn',...
    'mean_synergy_without_smn','mean_synergy_within_smn','intra_net_red_without_smn','intra_net_red_within_smn',...
    'intra_net_syn_without_smn','intra_net_syn_within_smn','extra_net_red_without_smn',...
    'extra_net_red_within_smn','extra_net_syn_without_smn','extra_net_syn_within_smn'});

ukb_RedSyn = [ukb_RedSyn,array2table(RedYeo7n_In),array2table(RedYeo7n_Btw),...
    array2table(SynYeo7n_In),array2table(SynYeo7n_Btw)];


x1 = netMeans1D(Red360p_NC,Yeo7MMP)./360;
x2 = netMeans1D(Syn360p_NC,Yeo7MMP)./360;
Y360p_redsyn_NC = [array2table(ukb_2_0_eID),array2table(x1),array2table(x2)];
Y360p_redsyn_NC.Properties.VariableNames{1} = 'eid';
Y360p_redsyn_NC.Properties.VariableNames{2} = 'Redundancy_VN';
Y360p_redsyn_NC.Properties.VariableNames{3} = 'Redundancy_SMN';
Y360p_redsyn_NC.Properties.VariableNames{4} = 'Redundancy_DAN';
Y360p_redsyn_NC.Properties.VariableNames{5} = 'Redundancy_VAN';
Y360p_redsyn_NC.Properties.VariableNames{6} = 'Redundancy_LN';
Y360p_redsyn_NC.Properties.VariableNames{7} = 'Redundancy_FPN';
Y360p_redsyn_NC.Properties.VariableNames{8} = 'Redundancy_DMN';
Y360p_redsyn_NC.Properties.VariableNames{9} = 'Synergy_VN';
Y360p_redsyn_NC.Properties.VariableNames{10} = 'Synergy_SMN';
Y360p_redsyn_NC.Properties.VariableNames{11} = 'Synergy_DAN';
Y360p_redsyn_NC.Properties.VariableNames{12} = 'Synergy_VAN';
Y360p_redsyn_NC.Properties.VariableNames{13} = 'Synergy_LN';
Y360p_redsyn_NC.Properties.VariableNames{14} = 'Synergy_FPN';
Y360p_redsyn_NC.Properties.VariableNames{15} = 'Synergy_DMN';
UKB_BrainXY = innerjoin(UKB_BrainMRICov,Y360p_redsyn_NC);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];




UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,innerjoin(ukb_RedSynYeo7n_2_0,UKB_AllVars)),UKB_BrainImagingVars_2_0);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+HeadMotion_2_0+SmokingStatus_2_0+AlcoholIntakeFreq_2_0+SleepHour_2_0+SNR_clean_2_0+',...
    'HBF_L_HBF_2_0+HBF_R_HBF_2_0+Handedness_LtoR+IntensityScaling_2_0+BloodPress_2_0+BloodPress_2_1+TDI_0_0+',...
    'BMI_2_0+Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];

XList = UKB_AllVars.Properties.VariableNames';
XList = XList(contains(XList,'_2_'));

XListAll = [{'Income_2_0';'EduAge_0_0'};XList];
YList = ukb_RedSyn.Properties.VariableNames(8:end)'; % Y.Properties.RedSyniableNames(2:end)';

results_Pheno = LMM_withFormula2(UKB_BrainXY,formula,XListAll,YList);




UKB_BrainXY = innerjoin(UKB_BrainMRICov,innerjoin(ukb_RedSynYeo7n_2_0,UKB_BrainImagingVars_2_0));
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = ukb_RedSynYeo7n_2_0.Properties.VariableNames(2:end)';

formula = ['Sex_0_0+HeadMotion_2_0+SmokingStatus_2_0+AlcoholIntakeFreq_2_0+SleepHour_2_0+SNR_clean_2_0+',...
    'HBF_L_HBF_2_0+HBF_R_HBF_2_0+Handedness_LtoR+IntensityScaling_2_0+BloodPress_2_0+BloodPress_2_1+TDI_0_0+',...
    'BMI_2_0+Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];

results_Age = LMM_withFormula2(UKB_BrainXY,formula,{'AgeAttend_2_0'},YList);


imagesc(RedSynT(1:193,1:8)'.*(RedSynT(1:193,1:8)'<0.05/193),[-15,15])
xticks(1:193);xticklabels(strrep(XListAll(1:193),'_',' '));xtickangle(90);
yticks(1:8);yticklabels([Yeo7Names;'BtwNetworks'])

imagesc(RedSynT(1:193,9:16)'.*(RedSynT(1:193,9:16)'<0.05/193),[-15,15])
xticks(1:193);xticklabels(strrep(XListAll(1:193),'_',' '));xtickangle(90);
yticks(1:8);yticklabels([Yeo7Names;'BtwNetworks'])

IQ_id = [124,151,152,158,160];
subplot(1,2,1)
imagesc(RedSynT(IQ_id,1:8)'.*(RedSynT(IQ_id,1:8)'<0.05/193),[-8,8])
yticks(1:8);yticklabels([Yeo7Names;'BtwNetworks'])
xticks(1:length(IQ_id));xticklabels(strrep(XListAll(IQ_id),'_',' '));xtickangle(90);
title('Red');colorbar;colormap(redbluecmap)
subplot(1,2,2)
imagesc(RedSynT(IQ_id,9:16)'.*(RedSynT(IQ_id,9:16)'<0.05/193),[-8,8])
yticks(1:8);yticklabels([Yeo7Names;'BtwNetworks'])
xticks(1:length(IQ_id));xticklabels(strrep(XListAll(IQ_id),'_',' '));xtickangle(90);
title('Syn');colorbar;colormap(redbluecmap)

plot_connmatrix(nanmean(synergy,3),DK_Yeo7,Yeo7Names,[2.43,3.8],'jet');colormap(jet(200))
plot_connmatrix(nanmean(redundancy,3),DK_Yeo7,Yeo7Names,[.04,1.18],'jet');colormap(jet(200))

[a, cb] = plot_cortical(parcel_to_surface(DK_Yeo7,'aparc_fsa5'),'color_range',[1,7]);
Yeo7Colormap = [0.471,0.0710,0.522;0.275,0.510,0.706;0,0.463,0.0550;...
    0.769,0.224,0.976;0.863,0.973,0.639;0.902,0.576,0.129;0.804,0.239,0.306];

colormap(Yeo7Colormap)

for i = 1:size(redundancy,3)
    [Rval(i,1),Pval(i,1)] = corr(icatb_mat2vec(redundancy(:,:,i)),icatb_mat2vec(synergy(:,:,i)));
end


for i = 1:7
    subplot(2,4,i)
    scatter(Age)
end


glasser_360_conte69 = load('glasser_360_fsa5.csv');% glasser_360_conte69 = glasser_360_conte69';
aparc_conte69 = load('aparc_fsa5.csv');
for j = 1:68
    for i = 1:360
        dice_matrix(i,j) = dice(aparc_conte69==j, glasser_360_conte69==i);
    end
end
for i = 1:360
    if max(dice_matrix(i,:)) ~= 0
        ROI2ROI(i,1) = find(dice_matrix(i,:) == max(dice_matrix(i,:)));
    end
end
ROI2ROI = nan(360,1);
for i = 1:360
    if max(dice_matrix(i,:)) ~= 0
        ROI2ROI(i,1) = find(dice_matrix(i,:) == max(dice_matrix(i,:)));
    end
end


% age, sex, body mass index (BMI), head motion, signal-to-noise ratio (SNR) 
% race, Townsend deprivation index (TDI)
% intracranial volume (ICV), data acquisition site as covariables

UKB_BrainImagingVars = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BrainImagingVars.csv');
UKB_GWC = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_GWC.csv');
UKB_BrainImagingVars = innerjoin(UKB_BrainImagingVars, UKB_GWC);
UKB_BrainImagingVars_2_0 = UKB_BrainImagingVars(:,[1,find(contains(UKB_BrainImagingVars.Properties.VariableNames,'_2_0'))]);

UKB_BrainXY = innerjoin(UKB_BrainMRICov,innerjoin(Y360p_redsyn_NC,UKB_BrainImagingVars_2_0));
y0y = UKB_BrainXY.Properties.VariableNames(contains(UKB_BrainXY.Properties.VariableNames,'Thickness'))';
y0y = y0y([2:34,36:68]);


Centre = [AsDummy(UKB_BrainXY.Centre_0_0)];

COV0 = [UKB_BrainXY{:,[2,10,20,24,29,32,67:70,77,81,83]},Centre];
COV0(:,13) = zscore( COV0(:,13));

x1 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'GWC')};
% x2 = netMeans1D(UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Syn')},ROI2ROI);
% x3 = netMeans1D(UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Red')},ROI2ROI);

x2 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Syn')}*dice_matrix;
x3 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Red')}*dice_matrix;

% x2(:,67:end)= [];x3(:,67:end)= [];x1(:,[4,39]) = [];

[r_syn p_syn] = partialcorr(x1,x2,COV0,'Rows','pairwise');
[r_red p_red] = partialcorr(x1,x3,COV0,'Rows','pairwise');
% indx = [1:3,5:34,35:38,40:68];
% rx = nan(68,68);rx(indx,indx) = r;
% px = nan(68,68);px(indx,indx) = p;

COV_LIST = UKB_BrainXY.Properties.VariableNames([2,10,20,24,29,32,67:70,77,81,83])';

max(abs(diag(r_red)))
[a, cb] = plot_cortical(parcel_to_surface(diag(r_red),'aparc_fsa5'),'color_range',[-0.05,0.05]);


DK_info_ukb = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/atlas/DK_info_ukb.csv');
DK_Yeo7 = nan(84,1);
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'visual')) = 1;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'somatomotor')) = 2;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'dorsal attention')) = 3;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'ventral attention')) = 4;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'limbic')) = 5;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'frontoparietal')) = 6;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'default mode')) = 7;
DK_Yeo7(69:84) = 8;


Centre = [AsDummy(UKB_BrainXY.Centre_0_0)];
Centre(:,[sum(Centre)<3500]) = [];

COV0 = [UKB_BrainXY{:,[2,10,20,24,29,32,67:70,77,81,83]},Centre];
COV0(:,13) = zscore( COV0(:,13));

x1 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Area')};
x3 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'BloodPress_0_0')};
[r p] = corr([nanmean(x1,2),nanmean(x2,2),x3],'Rows','pairwise');


x1 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Area')};
x2 = netMeans1D(UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Syn')},ROI2ROI);
x2(:,67:end)= [] ;
[r p] = partialcorr(x1,x2,COV0,'Rows','pairwise');
r = r([2:34,36:68],:);p = p([2:34,36:68],:);
rx = nan(66,68);rx(:,unique(ROI2ROI(~isnan(ROI2ROI)))) = r;
px = nan(66,68);px(:,unique(ROI2ROI(~isnan(ROI2ROI)))) = p;
rxx = nan(68,68);rxx([1:31,33:65,67:68],:) = rx;
pxx = nan(68,68);pxx([1:31,33:65,67:68],:) = px;
[a, cb] = plot_cortical(parcel_to_surface(diag(rxx),'aparc_fsa5'),'color_range',[-0.07,0.07]);

[a, cb] = plot_cortical(parcel_to_surface(diag(rxx).*(diag(pxx)<0.05/64),'aparc_fsa5'),'color_range',[-0.07,0.07]);


parpool('local',20);

UKB_BrainXY = innerjoin(UKB_BrainMRICov,innerjoin(ukb_RedSynYeo7n_2_0,UKB_BrainImagingVars_2_0));

YList = {'SynYeo7n_Total','RedYeo7n_Total'};
MList = {'WMHyperInt_WB_Total_2_0';'FA_WB_Total_2_0';'MD_WB_Total_2_0';'MO_WB_Total_2_0';...
    'ICVF_WB_Total_2_0';'ISOVF_WB_Total_2_0'};
COVList = {'Sex_0_0','HeadMotion_2_0','SNR_2_0','TDI_0_0',...
    'BMI_2_0','Race_1','Race_2','Race_3','Race_4','Vol_WB_TIV_2_0'};
COV = UKB_BrainXY{:,COVList};
for m = 1:length(MList)
    for i = 1:length({'AgeAttend_2_0'})
        for j = 1:length(YList)
            % clear MediationRes
            Xvar = UKB_BrainXY.AgeAttend_2_0;
            Yvar_Syn = UKB_BrainXY{:,YList{j}};
            Mvar = UKB_BrainXY{:,MList{m}};
            [X2,Y2,M2,COV2] = nanrm(Xvar,Yvar_Syn,Mvar,COV0);
            [paths, stats] = mediation(X2, Y2, M2,'covs',COV2);%
            MediateResult = [stats.paths;stats.mean;stats.ste;stats.t;stats.df;stats.p];
            MediateResults(:,:,i,j,m) = MediateResult;
            disp(['Mediate Analysis: X ','AgeAttend_2_0',' ','Y ',YList{j},'M ',MList{m}]);
        end
    end
end

MediateResults = squeeze(MediateResults);
size(MediateResults)
MediateResult = squeeze(MediateResults(:,:,1,1));

tM = squeeze(MediateResults(4,5,:,:,:));
pM = squeeze(MediateResults(5,5,:,:,:));

[r p] = partialcorr(UKB_BrainXY{:,YList([2,1])},UKB_BrainXY{:,MList},COV0,'Rows','pairwise');

figure;imagesc(r.*(p<0.001),[-0.08,0.08])
set(gca,'DataAspectRatio',[1 1 1]);colorbar;
set(gca,'DataAspectRatio',[1 1 1]);colorbar;colormap(redbluecmap)
Mx = {'WMHyper','FA','MD','MO','ICVF','ISOVF'}
set(gca, 'XTick', 1:12, 'XTickLabel', Mx);xtickangle(90)
set(gca, 'YTick', 1:2, 'YTickLabel', {'Redundancy','Synergy'});

figure;

r0 = r(2,:);

h = barh(r0);h.FaceColor = 'flat';h.CData = r0;
caxis([-0.05, 0.05]);xlim([-0.05, 0.05]);ylim([0.5, 6.5]);
original_cmap = redbluecmap;
continuous_cmap = interp1(1:size(original_cmap,1), original_cmap, linspace(1, size(original_cmap,1), 256));
colormap(continuous_cmap);

colorbar;grid on;
set(gca, 'YDir', 'reverse');

Mx = {'WMHyper','FA','MD','MO','ICVF','ISOVF'};
set(gca, 'YTick', 1:6, 'YTickLabel', Mx);


r0 = r(1,:);
figure;

theta = (0:5) * (2*pi/6);
width = 2*pi/6 * 0.6;
radius = abs(r0);
original_cmap = redbluecmap;
continuous_cmap = interp1(1:size(original_cmap,1), original_cmap, linspace(1, size(original_cmap,1), 256));

hold on;

for i = 1:length(r0)
  theta_start = theta(i) - width/2;
  theta_end = theta(i) + width/2;
  
  theta_sector = linspace(theta_start, theta_end, 50);
  
  for j = 1:length(theta_sector)-1
      theta_fill = [theta_sector(j), theta_sector(j+1), theta_sector(j+1), theta_sector(j), theta_sector(j)];
      r_fill = [0, 0, radius(i), radius(i), 0];
      
      [x, y] = pol2cart(theta_fill, r_fill);
      color_index = round((r0(i) + 0.05) / 0.1 * 255) + 1;
      color_index = max(1, min(256, color_index));
      
      fill(x, y, continuous_cmap(color_index, :), 'EdgeColor', 'none');
  end
end

axis equal;
xlim([-0.06, 0.06]);ylim([-0.06, 0.06]);

colormap(continuous_cmap);caxis([-0.05, 0.05]);% colorbar;

theta_labels = theta;
label_radius = 0.055;
Mx = {'WMHyper','FA','MD','MO','ICVF','ISOVF'};

for i = 1:length(theta_labels)
  [x_label, y_label] = pol2cart(theta_labels(i), label_radius);
  text(x_label, y_label, Mx{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 15);
end

for rs = 0.01:0.01:0.05
  theta_circle = linspace(0, 2*pi, 100);
  [x_circle, y_circle] = pol2cart(theta_circle, rs);
  plot(x_circle, y_circle, 'k--', 'LineWidth', 0.5, 'Color', [0.7 0.7 0.7]);
  text(0.002, rs, sprintf('%.2f', rs), 'FontSize', 12, 'Color', [0.1 0.1 0.1]);
end
for i = 1:length(theta)
  [x_line, y_line] = pol2cart([theta(i), theta(i)], [0, 0.05]);
  plot(x_line, y_line, 'k--', 'LineWidth', 0.5, 'Color', [0.7 0.7 0.7]);
end

grid off;axis off;

for i = 1:length(r0)
  color_idx = round((r0(i) + 0.05) / 0.1 * 255) + 1;
  color_idx = max(1, min(256, color_idx));
  fprintf('r0(%d) = %.3f, color', i, r0(i), color_idx);
end








UKB_AllVars = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_AllVars_NumCol.csv');
UKB_AllVarDic = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_AllVarUsed.xlsx');
UKB_AllVarDic.Properties.VariableNames{5} = 'Description_CN';
UKB_AllVarDic.Description_to_Colnames = replace_column_names(UKB_AllVarDic.Description);
UKB_AllVars = removevars(UKB_AllVars, {'MRI_completed_2_0','MRI_completed_3_0'});
XList= UKB_AllVars.Properties.VariableNames([29,33,37,43:45,198:4:279,452,547,551,558]-2)';
UKB_BrainXY = innerjoin(UKB_BrainMRICov,innerjoin(ukb_RedSynYeo7n_2_0,UKB_AllVars));
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];




UKB_GWC = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_GWC.csv');
UKB_BrainImagingVars = innerjoin(UKB_BrainImagingVars, UKB_GWC);
UKB_BrainImagingVars_2_0 = UKB_BrainImagingVars(:,[1,find(contains(UKB_BrainImagingVars.Properties.VariableNames,'_2_0'))]);

UKB_BrainXY = innerjoin(UKB_BrainMRICov,innerjoin(Y360p_redsyn_NC,UKB_BrainImagingVars_2_0));


x1 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'GWC')};
x2 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Syn')}*dice_matrix;
x3 = UKB_BrainXY{:,contains(UKB_BrainXY.Properties.VariableNames,'Red')}*dice_matrix;
Centre = [AsDummy(UKB_BrainXY.Centre_0_0)];
Centre(:,[sum(Centre)<3500]) = [];
COV0 = [UKB_BrainXY{:,[2,10,20,24,29,32,67:70,77,81,83]},Centre];
COV0(:,13) = zscore( COV0(:,13));

for i = 1
    for m = 1:68
        % clear MediationRes
        Xvar = COV0(:,3);
        Yvar_Syn = x2(:,m);Yvar_Red = x3(:,m);
        Mvar = x1(:,m);
        [X2,Y2,M2,COV2] = nanrm(Xvar,Yvar_Red,Mvar,COV0(:,[1,2,4:16]));
        [paths, stats] = mediation(X2, Y2, M2,'covs',COV2);%
        MediateResult = [stats.paths;stats.mean;stats.ste;stats.t;stats.df;stats.p];
        MediateResults_Red(:,:,m) = MediateResult;
        [X2,Y2,M2,COV2] = nanrm(Xvar,Yvar_Syn,Mvar,COV0(:,[1,2,4:16]));
        [paths, stats] = mediation(X2, Y2, M2,'covs',COV2);%
        MediateResult = [stats.paths;stats.mean;stats.ste;stats.t;stats.df;stats.p];
        MediateResults_Syn(:,:,m) = MediateResult;
        disp(['Mediate Analysis: X(AgeAttend_2_0) - M - Y ',num2str(m)]);
    end
end




