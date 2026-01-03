%% 2.0
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/atlas/HCPMMP_atlas_info.mat')

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/red_syn/ukb_2_0_glasser_synergy.mat')

for i = 1:size(synergy,3)
    Syn360p_PC(i,:) = participation_coef(synergy(:,:,i),Yeo7MMP);
    Syn360p_NC(i,:) = sum(synergy(:,:,i));
end

SynYeo7n_Total = mean(icatb_mat2vec3d(synergy),2);
[SynYeo7n] = NetMeanMatrix3D(synergy,Yeo7MMP);
SynYeo7n_Btw = icatb_mat2vec3d(SynYeo7n);
SynYeo7n_In = icatb_mat2dia3d(SynYeo7n);


clear synergy

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/red_syn/ukb_2_0_glasser_redundancy.mat')

for i = 1:size(redundancy,3)
    Red360p_PC(i,:) = participation_coef(redundancy(:,:,i),Yeo7MMP);
    Red360p_NC(i,:) = sum(redundancy(:,:,i));
end



RedYeo7n_Total = mean(icatb_mat2vec3d(redundancy),2);
[RedYeo7n] = NetMeanMatrix3D(redundancy,Yeo7MMP);
RedYeo7n_Btw = icatb_mat2vec3d(RedYeo7n);
RedYeo7n_In = icatb_mat2dia3d(RedYeo7n);

clear redundancy

ukb_SynYeo7n_2_0 = [array2table(ukb_2_0_eID),array2table(SynYeo7n_In),array2table(SynYeo7n_Btw),array2table(SynYeo7n_Total)];ukb_SynYeo7n_2_0.Properties.VariableNames{1} = 'eid';
ukb_RedYeo7n_2_0 = [array2table(ukb_2_0_eID),array2table(RedYeo7n_In),array2table(RedYeo7n_Btw),array2table(RedYeo7n_Total)];ukb_RedYeo7n_2_0.Properties.VariableNames{1} = 'eid';

ukb_RedSynYeo7n_2_0 = innerjoin(ukb_SynYeo7n_2_0,ukb_RedYeo7n_2_0);

SynYeo7n_PC = netMean1D(Syn360p_PC,Yeo7MMP);
RedYeo7n_PC = netMean1D(Red360p_PC,Yeo7MMP);



ukb_RedSynPCYeo7n_2_0.Properties.VariableNames{1} = 'eid';

UKB_BrainXY = innerjoin(UKB_UsedData_Brain,ukb_RedSynPCYeo7n_2_0);
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
%     ;
YList = ukb_RedSynPCYeo7n_2_0.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

[RedSynT,RedSynP,RedSynDF,RedSynBeta,RedSynCI1,RedSynCI2] = LMM_withFormula2(UKB_BrainXY,formula,XListAll,YList);

for i = 1:12
    RedSynLMMRes = [table(XListAll),array2table([RedSynBeta(:,i),RedSynT(:,i),RedSynDF(:,i),RedSynP(:,i)])];
    RedSynLMMRes.Properties.VariableNames(1:5) = {'X','Beta','TVal','DF','PVal'};
    writetable(RedSynLMMRes,[YList{i},'_forManhattan.csv']);
end




UKB_BrainXY = innerjoin(UKB_BrainMRICov,ukb_RedSynYeo7n_2_0);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
formula = ['Sex_0_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
YList = ukb_RedSynYeo7n_2_0.Properties.VariableNames(8:end)';
[RedSynAgeT,RedSynAgeP,RedSynAgeDF,RedSynAgeBeta] = LMM_withFormula2(UKB_BrainXY,formula,{'AgeAttend_2_0'},YList);





XRedSyn = [array2table(ukb_RedSynYeo7n_2_0.eid),ukb_RedSynYeo7n_2_0(:,[1,8:16])];
XRedSyn.Properties.VariableNames{2} = 'IID';XRedSyn.Properties.VariableNames{1} = 'FID';

writetable(XRedSyn,'/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/redsyn/UKB_redsyn_pheno.txt','Delimiter' ,' ')


XRedSyn = [brain_data(:,1),array2table(Red7n_NC),array2table(Syn7n_NC)];
XRedSyn.Properties.VariableNames{2} = 'IID';XRedSyn.Properties.VariableNames{1} = 'FID';
writetable(XRedSyn,'/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/redsyn/UKB_redsyn_pheno.txt','Delimiter' ,' ')





cd('/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/redsyn')
phenofile = '/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/redsyn/UKB_redsyn_pheno.txt';

covarfile = '/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/UKB_GeneCov.txt';
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

ukb_RedSynYeo7n_2_0 = removevars(ukb_RedSynYeo7n_2_0, {'mean_redundancy_without_smn','mean_redundancy_within_smn',...
    'mean_synergy_without_smn','mean_synergy_within_smn','intra_net_red_without_smn','intra_net_red_within_smn',...
    'intra_net_syn_without_smn','intra_net_syn_within_smn','extra_net_red_without_smn',...
    'extra_net_red_within_smn','extra_net_syn_without_smn','extra_net_syn_within_smn'});

ukb_RedSynYeo7n_2_0 = [ukb_RedSynYeo7n_2_0,array2table(RedYeo7n_In),array2table(RedYeo7n_Btw),...
    array2table(SynYeo7n_In),array2table(SynYeo7n_Btw)];





UKB_BrainXY = innerjoin(UKB_UsedData_Brain,ukb_RedSynYeo7n_2_0(:,[1,8:end]));
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
YList = ukb_RedSynYeo7n_2_0.Properties.VariableNames(8:end)'; % Y.Properties.RedSyniableNames(2:end)';

[RedSynT,RedSynP,RedSynDF,RedSynBeta,RedSynCI1,RedSynCI2] = LMM_withFormula2(UKB_BrainXY,formula,XListAll,YList);

formula = ['Sex_0_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
YList = ukb_RedSynYeo7n_2_0.Properties.VariableNames(8:end)';
[RedSynAgeT,RedSynAgeP,RedSynAgeDF,RedSynAgeBeta] = LMM_withFormula2(UKB_BrainXY,formula,{'AgeAttend_2_0'},YList);


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



%% 3.0
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/red_syn/ukb_3_0_glasser_red_syn.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/red_syn/ukb_2_0_glasser_synergy.mat', 'ukb_2_0_eID')

[x2y,y2x] = idFinderNum(ukb_2_0_eID, ukb_3_0_BOLD_eID);

clear Syn360p_PC_3_0 Syn360p_NC_3_0 Red360p_PC_3_0 Red360p_NC_3_0
for i = 1:size(synergy,3)
    Syn360p_PC_3_0(i,:) = participation_coef(synergy(:,:,i),Yeo7MMP);
    Syn360p_NC_3_0(i,:) = sum(synergy(:,:,i));
end
SynYeo7n_Total = mean(icatb_mat2vec3d(synergy),2);
[SynYeo7n] = NetMeanMatrix3D(synergy,Yeo7MMP);
SynYeo7n_Btw = icatb_mat2vec3d(SynYeo7n);
SynYeo7n_In = icatb_mat2dia3d(SynYeo7n);
clear synergy
for i = 1:size(redundancy,3)
    Red360p_PC_3_0(i,:) = participation_coef(redundancy(:,:,i),Yeo7MMP);
    Red360p_NC_3_0(i,:) = sum(redundancy(:,:,i));
end
RedYeo7n_Total = mean(icatb_mat2vec3d(redundancy),2);
[RedYeo7n] = NetMeanMatrix3D(redundancy,Yeo7MMP);
RedYeo7n_Btw = icatb_mat2vec3d(RedYeo7n);
RedYeo7n_In = icatb_mat2dia3d(RedYeo7n);
clear redundancy
ukb_SynYeo7n_3_0 = [array2table(ukb_3_0_BOLD_eID),array2table(SynYeo7n_In),array2table(SynYeo7n_Btw),array2table(SynYeo7n_Total)];ukb_SynYeo7n_3_0.Properties.VariableNames{1} = 'eid';
ukb_RedYeo7n_3_0 = [array2table(ukb_3_0_BOLD_eID),array2table(RedYeo7n_In),array2table(RedYeo7n_Btw),array2table(RedYeo7n_Total)];ukb_RedYeo7n_3_0.Properties.VariableNames{1} = 'eid';
ukb_RedSynYeo7n_3_0 = innerjoin(ukb_SynYeo7n_3_0,ukb_RedYeo7n_3_0);
SynYeo7n_PC_3_0 = netMean1D(Syn360p_PC_3_0,Yeo7MMP);
RedYeo7n_PC_3_0 = netMean1D(Red360p_PC_3_0,Yeo7MMP);
save UKB_glasser_RedSyn_Yeo7n
writetable(ukb_RedSynYeo7n_3_0,'ukb_RedSynYeo7n_3_0.txt')
ukb_SynPCYeo7n_3_0 = [array2table(ukb_3_0_BOLD_eID),array2table(SynYeo7n_PC_3_0)];ukb_SynYeo7n_3_0.Properties.VariableNames{1} = 'eid';
ukb_RedPCYeo7n_3_0 = [array2table(ukb_3_0_BOLD_eID),array2table(RedYeo7n_PC_3_0)];ukb_RedYeo7n_3_0.Properties.VariableNames{1} = 'eid';
ukb_RedPCSynYeo7n_3_0 = innerjoin(ukb_SynPCYeo7n_3_0,ukb_RedPCYeo7n_3_0);
writetable(ukb_RedPCSynYeo7n_3_0,'ukb_RedSynPCYeo7n_3_0.txt')





