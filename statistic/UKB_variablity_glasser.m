% Figure 1
% Calculate the variability 
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/fMRI/ukb_glasser_ts_2_0.mat')
for i = 1:length(fMRI_Glasser_TS_2_0)
    NaID(i,1) = (sum(sum(isnan(fMRI_Glasser_TS_2_0{i})))>0)|(isempty(fMRI_Glasser_TS_2_0{i}));
end
fMRI_Glasser_TS_2_0(NaID==1) = [];
fMRI_Glasser_ID_2_0(NaID==1) = [];
fMRI_Glasser_eID_2_0(NaID==1) = [];

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/fMRI/ukb_glasser_ts_3_0.mat')
clear NaID
for i = 1:length(fMRI_Glasser_TS_3_0)
    NaID(i,1) = (sum(sum(isnan(fMRI_Glasser_TS_3_0{i})))>0)|(isempty(fMRI_Glasser_TS_3_0{i}));
end
fMRI_Glasser_TS_3_0(NaID==1) = [];
fMRI_Glasser_ID_3_0(NaID==1) = [];
fMRI_Glasser_eID_3_0(NaID==1) = [];

parpool(20);
parfor i = 1:length(fMRI_Glasser_eID_2_0)
    [var_glasser_2_0(i,:)] = func_arch_var(fMRI_Glasser_TS_2_0{i}',20:10:50);
    disp(['Temporal Variability of UK Biobank: sub-',num2str(fMRI_Glasser_eID_2_0(i)),'.'])
end
save /public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_glasser_variability_2_0.mat var_glasser_2_0 fMRI_Glasser_eID_2_0

parpool(20);
parfor i = 1:length(fMRI_Glasser_eID_3_0)
    [var_glasser_3_0(i,:)] = func_arch_var(fMRI_Glasser_TS_3_0{i}',20:10:50);
    disp(['Temporal Variability of UK Biobank: sub-',num2str(fMRI_Glasser_eID_3_0(i)),'.'])
end
save /public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_glasser_variability_3_0.mat ...
    var_glasser_3_0 fMRI_Glasser_eID_3_0


%%

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/imaging/fMRI/ukb_Tian_Subcortex_ts.mat', 'fMRI_Tian_Subcortex_ID', 'fMRI_Tian_Subcortex_S1_TS', 'fMRI_Tian_Subcortex_eID')

fMRI_Tian_Subcortex_eID_2_0 = fMRI_Tian_Subcortex_eID(contains(fMRI_Tian_Subcortex_ID,"_2_0"));
fMRI_Tian_Subcortex_S1_TS_2_0 = fMRI_Tian_Subcortex_S1_TS(contains(fMRI_Tian_Subcortex_ID,"_2_0"));
fMRI_Tian_Subcortex_ID_2_0 = fMRI_Tian_Subcortex_ID(contains(fMRI_Tian_Subcortex_ID,"_2_0"));

fMRI_Tian_Subcortex_eID_3_0 = fMRI_Tian_Subcortex_eID(contains(fMRI_Tian_Subcortex_ID,"_3_0"));
fMRI_Tian_Subcortex_S1_TS_3_0 = fMRI_Tian_Subcortex_S1_TS(contains(fMRI_Tian_Subcortex_ID,"_3_0"));
fMRI_Tian_Subcortex_ID_3_0 = fMRI_Tian_Subcortex_ID(contains(fMRI_Tian_Subcortex_ID,"_3_0"));

clear fMRI_Tian_Subcortex_S1_TS fMRI_Tian_Subcortex_ID fMRI_Tian_Subcortex_eID

clear NaID
for i = 1:length(fMRI_Tian_Subcortex_S1_TS_2_0)
    NaID(i,1) = (sum(sum(isnan(fMRI_Tian_Subcortex_S1_TS_2_0{i})))>0)|(isempty(fMRI_Tian_Subcortex_S1_TS_2_0{i}));
end
fMRI_Tian_Subcortex_S1_TS_2_0(NaID==1) = [];
fMRI_Tian_Subcortex_ID_2_0(NaID==1) = [];
fMRI_Tian_Subcortex_eID_2_0(NaID==1) = [];

clear NaID
for i = 1:length(fMRI_Tian_Subcortex_S1_TS_3_0)
    NaID(i,1) = (sum(sum(isnan(fMRI_Tian_Subcortex_S1_TS_3_0{i})))>0)|(isempty(fMRI_Tian_Subcortex_S1_TS_3_0{i}));
end
fMRI_Tian_Subcortex_S1_TS_3_0(NaID==1) = [];
fMRI_Tian_Subcortex_ID_3_0(NaID==1) = [];
fMRI_Tian_Subcortex_eID_3_0(NaID==1) = [];

% 2-0
fMRI_common_eID_2_0 = intersect(fMRI_Glasser_eID_2_0,fMRI_Tian_Subcortex_eID_2_0);
fMRI_Tian_Subcortex_S1_TS_2_0_sCx = fMRI_Tian_Subcortex_S1_TS_2_0(ismember(fMRI_Tian_Subcortex_eID_2_0,fMRI_common_eID_2_0));
fMRI_Glasser_TS_2_0_sCx = fMRI_Glasser_TS_2_0(ismember(fMRI_Glasser_eID_2_0,fMRI_common_eID_2_0));

parpool(20);
parfor i = 1:length(fMRI_common_eID_2_0)
    [var_glasser_s1_2_0(i,:)] = func_arch_var_cort2sub(...
        fMRI_Glasser_TS_2_0_sCx{i}',fMRI_Tian_Subcortex_S1_TS_2_0_sCx{i}',20:10:50);
    disp(['Temporal Variability (Cort-SubCort) of UK Biobank: sub-',num2str(fMRI_common_eID_2_0(i)),'.'])
end
save /public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_glasser+s1_variability_2_0.mat ...
    var_glasser_s1_2_0 fMRI_common_eID_2_0


% 3-0
fMRI_common_eID_3_0 = intersect(fMRI_Glasser_eID_3_0,fMRI_Tian_Subcortex_eID_3_0);
fMRI_Tian_Subcortex_S1_TS_3_0_sCx = fMRI_Tian_Subcortex_S1_TS_3_0(ismember(fMRI_Tian_Subcortex_eID_3_0,fMRI_common_eID_3_0));
fMRI_Glasser_TS_3_0_sCx = fMRI_Glasser_TS_3_0(ismember(fMRI_Glasser_eID_3_0,fMRI_common_eID_3_0));

parpool(20);
parfor i = 1:length(fMRI_common_eID_3_0)
    [var_glasser_s1_3_0(i,:)] = func_arch_var_cort2sub(...
        fMRI_Glasser_TS_3_0_sCx{i}',fMRI_Tian_Subcortex_S1_TS_3_0_sCx{i}',20:10:50);
    disp(['Temporal Variability (Cort-SubCort) of UK Biobank: sub-',num2str(fMRI_common_eID_3_0(i)),'.'])
end
save /public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_glasser+s1_variability_3_0.mat ...
    var_glasser_s1_3_0 fMRI_common_eID_3_0

%%
% Figure 2

UKB_BrainMRICov = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BasicMRI.csv');
UKB_BrainImagingVars = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BrainImagingVars.csv');
UKB_AllVars = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_AllVars_NumCol.csv');
UKB_rfMRICov = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_rfMRIcov.csv');
UKB_Biomarker = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Biomaker.csv');
UKB_PRS = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_PRS.csv');
UKB_MentalDiagnosis = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_MentalDiagnosis.csv');


UKB_Biomarker.Vitamin_D_1_0 = str2double(UKB_Biomarker.Vitamin_D_1_0);
UKB_BrainImagingVars.Vol_WB_CWM_3_0 = str2double(UKB_BrainImagingVars.Vol_WB_CWM_3_0);
UKB_PRS.venous_thromboembolic_disease_VTE_1 = str2double(UKB_PRS.venous_thromboembolic_disease_VTE_1);
UKB_rfMRICov.SNR_3_0 = str2double(UKB_rfMRICov.SNR_3_0);

sum(UKB_BrainMRICov.Handedness_0_0==1)
sum(UKB_BrainMRICov.Handedness_0_0==2)
sum(UKB_BrainMRICov.Handedness_0_0==3)

diffIdx = (UKB_BrainMRICov.Handedness_0_0 - UKB_BrainMRICov.Handedness_2_0);
find((~isnan(diffIdx))&(diffIdx~=0))
UKB_BrainMRICov.Handedness = UKB_BrainMRICov.Handedness_0_0;
UKB_BrainMRICov.Handedness(find((~isnan(diffIdx))&(diffIdx~=0))) = NaN;

UKB_BrainMRICov.Handedness_R = double(UKB_BrainMRICov.Handedness==1);
UKB_BrainMRICov.Handedness_R(find(isnan(UKB_BrainMRICov.Handedness)))=NaN;

UKB_BrainMRICov.Handedness_L = double(UKB_BrainMRICov.Handedness==2);
UKB_BrainMRICov.Handedness_L(find(isnan(UKB_BrainMRICov.Handedness)))=NaN;

UKB_BrainMRICov.Handedness_M = UKB_BrainMRICov.Handedness;
UKB_BrainMRICov.Handedness_M(UKB_BrainMRICov.Handedness==3) = 1;
UKB_BrainMRICov.Handedness_M(UKB_BrainMRICov.Handedness~=3) = 0;
UKB_BrainMRICov.Handedness_M(find(isnan(UKB_BrainMRICov.Handedness)))=NaN;

UKB_BrainMRICov.Handedness_LtoR = UKB_BrainMRICov.Handedness;
UKB_BrainMRICov.Handedness_LtoR(UKB_BrainMRICov.Handedness==2) = 1;
UKB_BrainMRICov.Handedness_LtoR(UKB_BrainMRICov.Handedness==3) = 2;
UKB_BrainMRICov.Handedness_LtoR(UKB_BrainMRICov.Handedness==1) = 3;
UKB_BrainMRICov.Handedness_LtoR(find(isnan(UKB_BrainMRICov.Handedness)))=NaN;



load('/public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_Scheafer7n100p_variability_2_0.mat')
fMRI_Scheafer7n100p_eID_2_0 = fMRI_Scheafer_eID_2_0;
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_Scheafer7n200p_variability_2_0.mat')
fMRI_Scheafer7n200p_eID_2_0 = fMRI_Scheafer_eID_2_0;
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_Scheafer7n400p_variability_2_0.mat')
fMRI_Scheafer7n400p_eID_2_0 = fMRI_Scheafer_eID_2_0;
clear fMRI_Scheafer_eID_2_0

[a, cb] = plot_cortical(parcel_to_surface(nanmean(var_Scheafer7n100p_2_0),'schaefer_100_fsa5'),'cmap','jet','color_range',[.47,.8]);
[a, cb] = plot_cortical(parcel_to_surface(nanmean(var_Scheafer7n200p_2_0),'schaefer_200_fsa5'),'cmap','jet','color_range',[.49,.85]);
[a, cb] = plot_cortical(parcel_to_surface(nanmean(var_Scheafer7n400p_2_0),'schaefer_400_fsa5'),'cmap','jet','color_range',[.53,.89]);

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_glasser_variability_2_0.mat')


Y360p_roi = [array2table(fMRI_Glasser_eID_2_0),array2table(var_glasser_2_0)];Y360p_roi.Properties.VariableNames{1} = 'eid';
clear xl;for i = 1:360;xl{i,1} = ['roi_',num2str(i,'%03d')];end
Y360p_roi.Properties.VariableNames(2:end) = xl;
Y360p_roi = innerjoin(UKB_rfMRICov(:,[1,6,10]),Y360p_roi);Y360p_roi(Y360p_roi.HeadMotion_2_0>0.2,:) = [];

Y100p_roi = [array2table(fMRI_Scheafer7n100p_eID_2_0),array2table(var_Scheafer7n100p_2_0)];Y100p_roi.Properties.VariableNames{1} = 'eid';
clear xl;for i = 1:100;xl{i,1} = ['roi_',num2str(i,'%03d')];end
Y100p_roi.Properties.VariableNames(2:end) = xl;
Y100p_roi = innerjoin(UKB_rfMRICov(:,[1,6,10]),Y100p_roi);Y100p_roi(Y100p_roi.HeadMotion_2_0>0.2,:) = [];

Y200p_roi = [array2table(fMRI_Scheafer7n200p_eID_2_0),array2table(var_Scheafer7n200p_2_0)];Y200p_roi.Properties.VariableNames{1} = 'eid';
clear xl;for i = 1:200;xl{i,1} = ['roi_',num2str(i,'%03d')];end
Y200p_roi.Properties.VariableNames(2:end) = xl;
Y200p_roi = innerjoin(UKB_rfMRICov(:,[1,6,10]),Y200p_roi);Y200p_roi(Y200p_roi.HeadMotion_2_0>0.2,:) = [];

Y400p_roi = [array2table(fMRI_Scheafer7n400p_eID_2_0),array2table(var_Scheafer7n400p_2_0)];Y400p_roi.Properties.VariableNames{1} = 'eid';
clear xl;for i = 1:400;xl{i,1} = ['roi_',num2str(i,'%03d')];end
Y400p_roi.Properties.VariableNames(2:end) = xl;
Y400p_roi = innerjoin(UKB_rfMRICov(:,[1,6,10]),Y400p_roi);Y400p_roi(Y400p_roi.HeadMotion_2_0>0.2,:) = [];

clear xl;

[coef400p,score400p,~,~,exp400p] = pca(Y400p_roi{:,4:end});
[coef200p,score200p,~,~,exp200p] = pca(Y200p_roi{:,4:end});
[coef100p,score100p,~,~,exp100p] = pca(Y100p_roi{:,4:end});
[coef360p,score360p,~,~,exp360p] = pca(Y360p_roi{:,4:end});

common_eID = intersect(intersect(intersect(Y360p_roi.eid, Y400p_roi.eid),Y200p_roi.eid),Y100p_roi.eid);

colorSet = redbluecmap(4)
plot(exp100p(1:10),'-*','Color',colorSet(1,:));hold on;
plot(exp200p(1:10),'-*','Color',colorSet(2,:));hold on;
plot(exp360p(1:10),'-*','Color',colorSet(3,:));hold on;
plot(exp400p(1:10),'-*','Color',colorSet(4,:));
legend('yeo100','yeo200','glasser360','yeo400')

[rVal1,pVal1] = corr([score100p(ismember(Y100p_roi.eid,common_eID),1),score200p(ismember(Y200p_roi.eid,common_eID),1),...
    score360p(ismember(Y360p_roi.eid,common_eID),1),score400p(ismember(Y400p_roi.eid,common_eID),1)],'Rows','pairwise')

[rVal2,pVal2] = corr([nanmean(Y100p_roi{ismember(Y100p_roi.eid,common_eID),4:end},2),...
    nanmean(Y200p_roi{ismember(Y200p_roi.eid,common_eID),4:end},2),...
    nanmean(Y360p_roi{ismember(Y360p_roi.eid,common_eID),4:end},2),...
    nanmean(Y400p_roi{ismember(Y400p_roi.eid,common_eID),4:end},2)],'Rows','pairwise')

[rVal3,pVal3] = corr([...
    parcel_to_surface(nanmean(Y100p_roi{ismember(Y100p_roi.eid,common_eID),4:end},1),'schaefer_100_fsa5')',...
    parcel_to_surface(nanmean(Y200p_roi{ismember(Y200p_roi.eid,common_eID),4:end},1),'schaefer_200_fsa5')',...
    parcel_to_surface(nanmean(Y360p_roi{ismember(Y360p_roi.eid,common_eID),4:end},1),'glasser_360_fsa5')',...
    parcel_to_surface(nanmean(Y400p_roi{ismember(Y400p_roi.eid,common_eID),4:end},1),'schaefer_400_fsa5')'],'Rows','pairwise')

[rVal4,pVal4] = corr([parcel_to_surface(coef360p(:,1),'glasser_360_fsa5')',...
    parcel_to_surface(coef400p(:,1),'schaefer_400_fsa5')',...
    parcel_to_surface(coef360p(:,2),'glasser_360_fsa5')',...
    parcel_to_surface(coef400p(:,2),'schaefer_400_fsa5')',...
    parcel_to_surface(coef360p(:,3),'glasser_360_fsa5')',...
    parcel_to_surface(coef400p(:,3),'schaefer_400_fsa5')'],...
    'Rows','pairwise')






UKB_UsedData_Brain = innerjoin(innerjoin(innerjoin(innerjoin(innerjoin(innerjoin(UKB_BrainMRICov,...
    UKB_rfMRICov),UKB_AllVars),UKB_Biomarker),UKB_BrainImagingVars),UKB_MentalDiagnosis),UKB_PRS);


BrainTotalList = {'FA_WB_Total_2_0','MD_WB_Total_2_0','MO_WB_Total_2_0','L1_WB_Total_2_0','L2_WB_Total_2_0',...
    'L3_WB_Total_2_0','ICVF_WB_Total_2_0','ISOVF_WB_Total_2_0','Area_WB_Total_2_0','Thickness_WB_Total_2_0',...
    'WMHyperInt_WB_Total_2_0','Vol_WB_SCG_2_0','Vol_WB_CSF_2_0','Vol_WB_Vent3_2_0','Vol_WB_Vent4_2_0',...
    'Vol_WB_Vent5_2_0','Vol_WB_Cort_2_0','Vol_WB_CWM_2_0'}';

PRSList = UKB_PRS.Properties.VariableNames(2:end)';
MentalList = UKB_MentalDiagnosis.Properties.VariableNames(2:end)';

BrainList = UKB_BrainImagingVars.Properties.VariableNames(2:end)';
BrainAsegList = BrainList((contains(BrainList,'Vol'))&(contains(BrainList,'_2_0')));
BrainAsegList = BrainAsegList([31:38,47:54]);

parpool(20);

clear xl
Y360p_roi = [array2table(fMRI_Glasser_eID_2_0),array2table(var_glasser_2_0)];
Y360p_roi.Properties.VariableNames{1} = 'eid';
for i = 1:360;xl{i,1} = ['roi_',num2str(i,'%03d')];end
Y360p_roi.Properties.VariableNames(2:361) = xl;
Y360p_roi.MeanVar = mean(Y360p_roi{:,2:361},2);
Y360p_roi = innerjoin(UKB_rfMRICov(:,[1,6,10]),Y360p_roi);
Y360p_roi(Y360p_roi.HeadMotion_2_0>0.2,:) = [];

X360p_roi = Y360p_roi{:,4:363};
[coef360p,score360p,~,~,exp360p] = pca(X360p_roi);
Y360p_pca = [Y360p_roi(:,1:3),array2table(score360p(:,1:5))];
Y360p_pca.Properties.VariableNames(4:end) = {'PC1','PC2','PC3','PC4','PC5'};

Y360p_band = [Y360p_roi(:,1:3),...
    array2table([nanmean(X360p_roi(:,coef360p(:,1)>0.06),2),...% B1: Low-order
    nanmean(X360p_roi(:,(coef360p(:,1)>0.04)&(coef360p(:,1)<=0.06)),2),...
    nanmean(X360p_roi(:,coef360p(:,1)<=0.04),2)])];% B3: High-order
Y360p_band.Properties.VariableNames(4:end) = {'Band1','Band2','Band3'};

X360p = nan(360,1);
X360p(coef360p(:,1)>0.06) = 1;
X360p((coef360p(:,1)>0.04)&(coef360p(:,1)<=0.06)) = 2;
X360p(coef360p(:,1)<=0.04) = 3;
[a, cb] = plot_cortical(parcel_to_surface(X360p,'glasser_360_fsa5'),'color_range',[1,3]);colormap([1,0,0;0,1,0;0,0,1])

set(gcf, 'Position', [650, 0, 1250, 980]);% Big Figure Size


Y400p_roi = [array2table(fMRI_Scheafer7n400p_eID_2_0),array2table(var_Scheafer7n400p_2_0)];
Y400p_roi.Properties.VariableNames{1} = 'eid';
for i = 1:400;xl{i,1} = ['roi_',num2str(i,'%03d')];end
Y400p_roi.Properties.VariableNames(2:401) = xl;
Y400p_roi.MeanVar = mean(Y400p_roi{:,2:401},2);
Y400p_roi = innerjoin(UKB_rfMRICov(:,[1,6,10]),Y400p_roi);
Y400p_roi(Y400p_roi.HeadMotion_2_0>0.2,:) = [];

X400p_roi = Y400p_roi{:,4:403};
[coef400p,score400p,~,~,exp400p] = pca(X400p_roi);
Y400p_pca = [Y400p_roi(:,1:3),array2table(score400p(:,1:5))];
Y400p_pca.Properties.VariableNames(4:end) = {'PC1','PC2','PC3','PC4','PC5'};

Y400p_band = [Y400p_roi(:,1:3),...
    array2table([nanmean(X400p_roi(:,coef400p(:,1)>0.06),2),...% B1: Low-order
    nanmean(X400p_roi(:,(coef400p(:,1)>0.04)&(coef400p(:,1)<=0.06)),2),...
    nanmean(X400p_roi(:,coef400p(:,1)<=0.04),2)])];% B3: High-order
Y400p_band.Properties.VariableNames(4:end) = {'Band1','Band2','Band3'};


X400p = nan(400,1);
X400p(coef400p(:,1)>0.06) = 1;
X400p((coef400p(:,1)>0.04)&(coef400p(:,1)<=0.06)) = 2;
X400p(coef400p(:,1)<=0.04) = 3;
[a, cb] = plot_cortical(parcel_to_surface(X400p,'schaefer_400_fsa5'),'color_range',[1,3]);colormap([1,0,0;0,1,0;0,0,1])

Y400p = [Y400p_pca,Y400p_band(:,4:6)];
Y360p = [Y360p_pca,Y360p_band(:,4:6)];

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_glasser+s1_variability_2_0.mat')

clear xl
YS1_roi = [array2table(fMRI_common_eID_2_0),array2table(var_glasser_s1_2_0)];
YS1_roi.Properties.VariableNames{1} = 'eid';
for i = 1:16;xl{i,1} = ['roi_',strrep(S1_ROI{i},'-','_')];end
YS1_roi.Properties.VariableNames(2:17) = xl;
YS1_roi = innerjoin(UKB_rfMRICov(:,[1,6,10]),YS1_roi);
YS1_roi(YS1_roi.HeadMotion_2_0>0.2,:) = [];


writetable(Y360p,'ukb_varpheno_glasser360p.csv')
writetable(Y400p,'ukb_varpheno_scheafer400p.csv')


UKB_BrainXY = innerjoin(UKB_UsedData_Brain,Y360p(:,[1,4:11]));
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = [UKB_AllVars.Properties.VariableNames(3:end)';...
    UKB_Biomarker.Properties.VariableNames(2:end)'];
XList = XList(contains(XList,'_2_'));
XList = XList(~(contains(XList,'_2_0_1')|contains(XList,'Intensity_')|...
    contains(XList,'VolRatio_')|contains(XList,'HBF_')|contains(XList,'Vol_WB_TIV_2_0')));
XList = XList(~(contains(XList,'_0_0')|contains(XList,'_1_0')|contains(XList,'_3_0')));
XListAll = [{'TDI_0_0';'Income_2_0';'EduAge_0_0';'Handedness_LtoR'};...
    XList;BrainTotalList;BrainAsegList;PRSList;MentalList(sum(UKB_UsedData_Brain{:,2195:2210},1)>200)];
YList = {'PC1','PC2','PC3','PC4','PC5','Band1','Band2','Band3'}; % Y.Properties.VariableNames(2:end)';
[VpcT,VpcP,VpcDF,VpcBeta,VpcCI1,VpcCI2] = LMM_withFormula2(UKB_BrainXY,formula,XListAll,YList);

UKB_BrainXY = innerjoin(UKB_UsedData_Brain,Y360p(:,[1,4:end]));
formula = ['Sex_0_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
YList = Y360p.Properties.VariableNames(4:end)'; 
[VpcTAge,VpcPAge,VpcDFAge,VpcBetaAge,VpcCI1Age,VpcCI2Age] = LMM_withFormula2(UKB_BrainXY,formula,{'AgeAttend_2_0'},YList);

UKB_BrainXY = innerjoin(UKB_UsedData_Brain,Y360p(:,[1,4:end]));
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
YList= Y360p.Properties.VariableNames(4:end)'; 
[VpcTSex,VpcPSex,VpcDFSex,VpcBetaSex,VpcCI1Sex,VpcCI2Sex] = LMM_withFormula(UKB_BrainXY,formula,{'Sex_0_0'},YList);



UKB_BrainXY = innerjoin(UKB_UsedData_Brain,YS1_roi(:,[1,4:end]));
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = [UKB_AllVars.Properties.VariableNames(3:end)';...
    UKB_Biomarker.Properties.VariableNames(2:end)'];
XList = XList(contains(XList,'_2_'));
XList = XList(~(contains(XList,'_2_0_1')|contains(XList,'Intensity_')|...
    contains(XList,'VolRatio_')|contains(XList,'HBF_')|contains(XList,'Vol_WB_TIV_2_0')));
XList = XList(~(contains(XList,'_0_0')|contains(XList,'_1_0')|contains(XList,'_3_0')));
XListAll = [{'TDI_0_0';'Income_2_0';'EduAge_0_0';'Handedness_LtoR'};...
    XList;BrainTotalList;BrainAsegList;PRSList;MentalList(sum(UKB_UsedData_Brain{:,2195:2210},1)>200)];
YList = xl; % Y.Properties.VariableNames(2:end)';
[Vs1T,Vs1P,Vs1DF,Vs1Beta,Vs1CI1,Vs1CI2] = LMM_withFormula2(UKB_BrainXY,formula,XListAll,YList);

UKB_BrainXY = innerjoin(UKB_UsedData_Brain,YS1_roi(:,[1,4:end]));
formula = ['Sex_0_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
YList = YS1_roi.Properties.VariableNames(4:end)'; 
[Vs1TAge,Vs1PAge,Vs1DFAge,Vs1BetaAge,Vs1CI1Age,Vs1CI2Age] = LMM_withFormula2(UKB_BrainXY,formula,{'AgeAttend_2_0'},YList);

UKB_BrainXY = innerjoin(UKB_UsedData_Brain,YS1_roi(:,[1,4:end]));
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
YList= YS1_roi.Properties.VariableNames(4:end)'; 
[Vs1TSex,Vs1PSex,Vs1DFSex,Vs1BetaSex,Vs1CI1Sex,Vs1CI2Sex] = LMM_withFormula(UKB_BrainXY,formula,{'Sex_0_0'},YList);


Vs1TSig = Vs1T'.*(Vs1P'<10^-8);SigID = (sum(Vs1TSig==0,1)~=16)&(sum(isnan(Vs1TSig),1)==0);
S1Order = ([1:8;9:16]);S1Order = S1Order(:);
imagesc(Vs1TSig(S1Order,SigID)',[-20,20]);
xticks(1:16);xticklabels(S1_ROI(S1Order))
yticks(1:sum(SigID));yticklabels(strrep(XListAll(SigID),'_',' '));
xtickangle(90);
set(gca,'DataAspectRatio',[1 1 1]);
title('Variability of Tian S1');colorbar;colormap(redbluecmap)


for i = 1:5
    VarLMMRes = [table(XList),array2table([VpcBeta(:,i),VpcT(:,i),VpcDF(:,i),VpcP(:,i)])];
    VarLMMRes.Properties.VariableNames(1:5) = {'X','Beta','TVal','DF','PVal'};
    writetable(VarLMMRes,['PC',num2str(i),'forManhattan_glasser360p.csv']);
end

for i = 1:3
    n = i + 5;
    VarLMMRes = [table(XList),array2table([VpcBeta(:,n),VpcT(:,n),VpcDF(:,n),VpcP(:,n)])];
    VarLMMRes.Properties.VariableNames(1:5) = {'X','Beta','TVal','DF','PVal'};
    writetable(VarLMMRes,['Band',num2str(i),'forManhattan_glasser360p.csv']);
end




for i = 1:16
    VarLMMRes = [table(XListAll),array2table([Vs1Beta(:,i),Vs1T(:,i),Vs1DF(:,i),Vs1P(:,i)])];
    VarLMMRes.Properties.VariableNames(1:5) = {'X','Beta','TVal','DF','PVal'};
    % VarLMMRes(VarLMMRes.PVal>10^-5,:) = [];
    writetable(VarLMMRes,['ROI',num2str(i),'forManhattan_tian16p.csv']);
end




UKB_BrainXY = innerjoin(UKB_UsedData_Brain,Y400p(:,[1,4:11]));
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
YList = {'PC1','PC2','PC3','PC4','PC5','Band1','Band2','Band3'}; % Y.Properties.VariableNames(2:end)';

[VarT400p,VarP400p,VarDF400p,VarBeta400p,VarCI1400p,VarCI2400p] = LMM_withFormula2(UKB_BrainXY,formula,XListAll,YList);



UKB_Protein = readtable('/public/home/zhangjie/UKB_UsedData_Protein.csv');
UKB_UsedData_Protein = innerjoin(UKB_rfMRICov,UKB_Protein);


UKB_UsedData_Protein.Race_1 = strrep(UKB_UsedData_Protein.Race_1,'NaN','');
UKB_UsedData_Protein.Race_2 = strrep(UKB_UsedData_Protein.Race_2,'NaN','');
UKB_UsedData_Protein.Race_3 = strrep(UKB_UsedData_Protein.Race_3,'NaN','');
UKB_UsedData_Protein.Race_4 = strrep(UKB_UsedData_Protein.Race_4,'NaN','');
UKB_UsedData_Protein.Race_1 = str2double(UKB_UsedData_Protein.Race_1);
UKB_UsedData_Protein.Race_2 = str2double(UKB_UsedData_Protein.Race_2);
UKB_UsedData_Protein.Race_3 = str2double(UKB_UsedData_Protein.Race_3);
UKB_UsedData_Protein.Race_4 = str2double(UKB_UsedData_Protein.Race_4);


UKB_BrainXY = innerjoin(UKB_UsedData_Protein,Y);
UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+Race_1+Race_2+Race_3+Race_4+AgeGap+',...
    'fastingtime_0_0+BMI_0_0+BloodPress_0_0+(1|Centre_0_0)'];
XList = UKB_Protein.Properties.VariableNames(71:end)';
YList = {'PC1','PC2','PC3','PC4','PC5'}; % Y.Properties.VariableNames(2:end)';
[VarPT,VarPP,VarPDF,VarPBeta,VarPCI1,VarPCI2] = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);

covar = [strcmp(UKB_BrainXY.Sex_0_0,'1'),UKB_BrainXY.AgeAttend_2_0,UKB_BrainXY.AgeAttend_2_0.^2,...
    UKB_BrainXY.HeadMotion_2_0,UKB_BrainXY.SNR_2_0,UKB_BrainXY.Race_1,UKB_BrainXY.Race_2,...
    UKB_BrainXY.Race_3,UKB_BrainXY.Race_4+UKB_BrainXY.AgeGap,AsDummy(UKB_BrainXY.batch),...
    UKB_BrainXY.fastingtime_0_0,UKB_BrainXY.BMI_0_0,UKB_BrainXY.BloodPress_0_0,AsDummy(UKB_BrainXY.Centre_0_0)];

parfor i = 1:size(UKB_BrainXY{:,81:3000},2)
    [VarPR2(i,:),VarPP2(i,:)] = partialcorr(UKB_BrainXY{:,80+i},UKB_BrainXY{:,3001:3005},covar,'Rows','pairwise');
    disp(['partical corr of ',XList{i},' & brain have been calculated.']);
end


% Variability PC1* Age (R)

% Variability by network * Age  (R)

% Advanced/Low level network (5 bands according to PC1) variability development trajectory (R)

% Figure 3
% Variability PC1* Cognition &  Mental Health
% Variability PC1* Physiological Variables & Disease 
% Variability PC1* Lifestyle & SES
% (Manhattan Plot: )
% Cox Regression of Disease
% (Manhattan Plot)

% SM 
% Variability ROI* Behavioral Phenotype 
% Variability ROI* Physiological Disease 
% Variability ROI* Lifestyle  
% (Heatmap)


% Mediating analysis (verifying the lifestyle-variable-health link)




UKB_BrainXY = innerjoin(UKB_UsedData_Brain,Y360p(:,[1,4:11]));
COV = [UKB_BrainXY.Sex_0_0,UKB_BrainXY.AgeAttend_2_0,(UKB_BrainXY.AgeAttend_2_0).^2,...
    UKB_BrainXY.HeadMotion_2_0,UKB_BrainXY.Race_1,UKB_BrainXY.Race_2,UKB_BrainXY.Race_3,...
    UKB_BrainXY.Race_4,UKB_BrainXY.Vol_WB_TIV_2_0,AsDummy(UKB_BrainXY.Centre_0_0)];

XList = [UKB_AllVars.Properties.VariableNames(3:end)'];
% XList = XList(contains(XList,'_0_'));XList = XList(~(contains(XList,'_0_0_1')));
XList = XList(~(contains(XList,'_3_')));
XNum = [1;2;3;4;5;6;7;8;9;10;12;14;16;17;18;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;...
    61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;...
    101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;...
    131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;220;221;222;223;224;225;226;227;228;229;230;...
    231;232;233;234;235;236;237;238;239;240;241;242;243;244;245;246;247;248;249;250;251;252;253;254;255;256;257;258;260;...
    264;267;270;273;276;279;282;285;286;287;288;289;290;291;292;293;294;295;296;297;298;299;300;...
    301;302;303;304;305;306;307;308;309;310;311;312;314;316;318;320;322;324;325;326;327;338;387;388;389;390;...
    391;392;393;394;395;396;397;398;399;400;401;408;409;410;411;412;413;414;415;416;417;418;419;420;...
    421;422;423;424;425;426;427;434;435;436;440;441;442;443;444;445;573;574;575;577;578;579;580;581;582;583;584;585;586;587;594;595;596;597];
XList = XList(XNum);
XListAll_0_0 = [{'TDI_0_0';'Income_0_0';'EduAge_0_0'};XList(contains(XList,'_0_0'))];
XListAll_2_0 = [XList(contains(XList,'_2_0'))];

MList = {'PC1'};%,'PC2','PC3','PC4','PC5','Band1','Band2','Band3'}; % Y.Properties.VariableNames(2:end)';
YNum = [13;14;15;16;17;18;19;20;21;22;23;24;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;...
    116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;...
    141;142;143;200;201;202;203;204;205;206;207;208;209;210;211;212;213;214;215;216;217;218;219;220;...
    221;222;223;224;225;226;227;228;229;230;231;232;233;234;235;236;237;238;239;240;241;242;252;253;258;259;...
    261;262;263;264;266;267;274;275;293;294];
YList = [UKB_AllVars.Properties.VariableNames(3:end)'];
YList = YList(~(contains(YList,'_0_')|contains(YList,'_1_')));
YList = YList(YNum);

YList_2_0 = [YList(contains(YList,'_2_0'))];
YList_3_0 = [YList(contains(YList,'_3_0'))];


for m = 1:length(MList)
    for j = 1:length(YList_3_0)
        parfor i = 1:length(XListAll_0_0) 
            % clear MediationRes
            Xvar = UKB_BrainXY{:,XListAll_0_0{i}};
            Yvar = UKB_BrainXY{:,YList_3_0{j}};
            Mvar = UKB_BrainXY{:,MList{m}};
            [X2,Y2,M2,COV2] = nanrm(Xvar,Yvar,Mvar,COV);
            [paths, stats] = mediation(X2, Y2, M2,'covs',COV2);%
            MedRes = [stats.paths;stats.mean;stats.ste;stats.t;stats.p];
            VarMed_BLtFU3(:,:,i,j,m) = MedRes;
            disp(['Mediation Analysis of UK Biobank: X(',XListAll_0_0{i}, ')-M(',MList{m},')-Y(',YList_3_0{j},').'])
        end
    end
end

V = squeeze(VarMed_BLtFU3(5,5,:,:,1));
imagesc(-1*log10(V(:,1:29))')
xticks(1:length(XListAll_0_0));xticklabels(strrep(XListAll_0_0,'_',' '));xtickangle(90);
yticks(1:length(YList_3_0));yticklabels(strrep(YList_3_0,'_',' '));

for m = 1:length(MList)
    for j = 1:length(YList_3_0)
        parfor i = 1:length(PRSList) 
            % clear MediationRes
            Xvar = UKB_BrainXY{:,PRSList{i}};
            Yvar = UKB_BrainXY{:,YList_3_0{j}};
            Mvar = UKB_BrainXY{:,MList{m}};
            [X2,Y2,M2,COV2] = nanrm(Xvar,Yvar,Mvar,COV);
            [paths, stats] = mediation(X2, Y2, M2,'covs',COV2);%
            MedRes = [stats.paths;stats.mean;stats.ste;stats.t;stats.p];
            VarMed_PRStFU3(:,:,i,j,m) = MedRes;
            disp(['Mediation Analysis of UK Biobank: X(',PRSList{i}, ')-M(',MList{m},')-Y(',YList_3_0{j},').'])
        end
    end
end

% Figure 4
% Variable PC1* Protein Brain Region * Protein (Annex) 
% Variable PC1* gene 
% (Manhanttan Plot + Gene loci + Hertiability + Genetic Overlap)
% (FUMA gene mapping + Overlap with Disease + MR Plot + Venn)

% Variability mode *AHBA 
% (Bar Plot + Cell Type + Enrichment Bubble)
X360p = [array2table(Y.eid),Y400p_roi,Y(:,2:6),Y_band(:,2:4),Y_net(:,2:8)];
X360p.Properties.VariableNames{2} = 'IID';X360p.Properties.VariableNames{1} = 'FID';

writetable(X360p,'/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/variability/UKB_variability_pheno.txt','Delimiter' ,' ')

cd('/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/variability')
phenofile = '/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/variability/UKB_variability_pheno.txt';

covarfile = '/public/home/zhangjie/DataAnalysis/wxr_UKB/gene_data/UKB_GeneCov.txt';
covarname = 'sex,AgeAttend_2_0,gene_pc1,gene_pc2,gene_pc3,gene_pc4,gene_pc5,gene_pc6,gene_pc7,gene_pc8,gene_pc9,gene_pc10';

parfor i = 1:5
    for chr = 1:22
        phenoname = ['PC',num2str(i)];
        bfile_qc = ['/home1/ISTBI_data/K_J_J/BackupData/UKBimputeV3_wbs/UKB_gene_v3_imp_qc_chr',...
            num2str(chr)];
        outfile = ['gwas_PC',num2str(i),'_chr',num2str(chr)];
        plink_command = ...
            ['/public/home/zhangjie/DataAnalysis/wxr_toolbox/toobox_other/plink2_linux_x86_64_20191030/plink2 ',...
            '--bfile ',bfile_qc,' ',...
            '--linear  hide-covar --debug ',...
            '--pheno ',phenofile,' ',...
            '--pheno-name ','PC',num2str(i),' '...
            '--covar ',covarfile,' ',...
            '--covar-name ',covarname,' ',...
            '--out ',outfile];
        unix(plink_command);
        disp(['GWAS of PC',num2str(i),' in chr',num2str(chr)])
    end
end

for i = 1:5
    
    unix(['head -n 1 gwas_PC',num2str(i),'_chr1.PC',num2str(i),'.glm.linear > gwas.PC',num2str(i),'.glm.linear']);
    unix(['tail -n +2 gwas_PC',num2str(i),'_chr1.PC',num2str(i),'.glm.linear >> gwas.PC',num2str(i),'.glm.linear']);
    
    for chr = 2:22
        unix(['tail -n +2 gwas_PC',num2str(i),'_chr',num2str(chr),'.PC',num2str(i),'.glm.linear >> gwas.PC',num2str(i),'.glm.linear']);
    end
    
    unix(['sed -i ','''','1s/.*/CHR POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P/',...
        '''',' gwas.PC',num2str(i),'.glm.linear']);
    
    unix(['gzip gwas.PC',num2str(i),'.glm.linear'])
end


parfor i = 1:3
    for chr = 1:22
        phenoname = ['Band',num2str(i)];
        bfile_qc = ['/home1/ISTBI_data/K_J_J/BackupData/UKBimputeV3_wbs/UKB_gene_v3_imp_qc_chr',...
            num2str(chr)];
        outfile = ['gwas_Band',num2str(i),'_chr',num2str(chr)];
        plink_command = ...
            ['/public/home/zhangjie/DataAnalysis/wxr_toolbox/toobox_other/plink2_linux_x86_64_20191030/plink2 ',...
            '--bfile ',bfile_qc,' ',...
            '--linear  hide-covar --debug ',...
            '--pheno ',phenofile,' ',...
            '--pheno-name ','Band',num2str(i),' '...
            '--covar ',covarfile,' ',...
            '--covar-name ',covarname,' ',...
            '--out ',outfile];
        unix(plink_command);
        disp(['GWAS of Band',num2str(i),' in chr',num2str(chr)])
    end
end


for i = 1:3
    
    unix(['head -n 1 gwas_Band',num2str(i),'_chr1.Band',num2str(i),'.glm.linear > gwas.Band',num2str(i),'.glm.linear']);
    unix(['tail -n +2 gwas_Band',num2str(i),'_chr1.Band',num2str(i),'.glm.linear >> gwas.Band',num2str(i),'.glm.linear']);
    
    for chr = 2:22
        unix(['tail -n +2 gwas_Band',num2str(i),'_chr',num2str(chr),'.Band',num2str(i),'.glm.linear >> gwas.Band',num2str(i),'.glm.linear']);
    end

    unix(['sed -i ','''','1s/.*/CHR POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P/',...
        '''',' gwas.Band',num2str(i),'.glm.linear']);
    
    unix(['gzip gwas.Band',num2str(i),'.glm.linear'])
end

% Gene + protein Protein interaction network





% Variability PC1* Mendelian randomization (validation of lifestyle-variability - health link)
 
% Figure 5
% Variability PC *Gradient mode 
% (5 Scatter Plot)
% Variability mode *Neurotransmitter
% (bar Plot)
% Variability mode *Neuromaps
% (Neurosynth Heatmap + Word Cloud + Other Plot Scatter)

% Verification 
% Whole brain signal regression/non-regression

load('/public/home/zhangjie/DataAnalysis/wxr_UKB/ukb_GS_ts.mat', 'fMRI_GS_eID_2_0', 'fMRI_GS_TS_2_0');
[common_eID,rRm,gRm] = findMissingNumbersIndices(fMRI_Glasser_eID_2_0, fMRI_GS_eID_2_0);

fMRI_Glasser_TS_2_0(rRm) = [];fMRI_GS_TS_2_0(gRm) = [];

fMRI_Glasser_eID_2_0(rRm) = [];fMRI_GS_eID_2_0(gRm) = [];

parfor i = 1:length(fMRI_Glasser_eID_2_0)
    fMRI_Glasser_TS_2_0{i,1} = linearFitResiduals(fMRI_Glasser_TS_2_0{i,1}',fMRI_GS_TS_2_0{i,1}')';
    disp(['Removing Global Signal of UK Biobank: sub-',num2str(fMRI_Glasser_eID_2_0(i)),'.'])
end

parfor i = 1:length(fMRI_Glasser_eID_2_0)
    [var_glasser_noGlobal(i,:)] = func_arch_var(fMRI_Glasser_TS_2_0{i}',30:40:60);
    disp(['Temporal Variability of UK Biobank: sub-',num2str(fMRI_Glasser_eID_2_0(i)),'.'])
end




% Longitudual

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_glasser_variability_3_0.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/result/xinranwu/ukb_glasser+s1_variability_3_0.mat')

score360p_followup = var_glasser_3_0*coef360p;

Y360p_roi_followup = [array2table(fMRI_Glasser_eID_3_0),array2table(var_glasser_3_0)];
Y360p_roi_followup.Properties.VariableNames{1} = 'eid';
for i = 1:360;xl{i,1} = ['roi_',num2str(i,'%03d')];end
Y360p_roi_followup.Properties.VariableNames(2:361) = xl;
Y360p_roi_followup.MeanVar = mean(Y360p_roi_followup{:,2:361},2);
Y360p_roi_followup = innerjoin(UKB_rfMRICov(:,[1,7,11]),Y360p_roi_followup);

Y360p_pca_followup = [Y360p_roi_followup(:,1:3),array2table(score360p_followup(:,1:5))];
Y360p_pca_followup.Properties.VariableNames(4:end) = {'PC1','PC2','PC3','PC4','PC5'};

X360p_roi_followup = Y360p_roi_followup{:,4:end};
Y360p_band_followup = [Y360p_roi_followup(:,1:3),...
    array2table([nanmean(X360p_roi_followup(:,coef360p(:,1)>0.06),2),...% B1: Low-order
    nanmean(X360p_roi_followup(:,(coef360p(:,1)>0.04)&(coef360p(:,1)<=0.06)),2),...
    nanmean(X360p_roi_followup(:,coef360p(:,1)<=0.04),2)])];% B3: High-order
Y360p_band_followup.Properties.VariableNames(4:end) = {'Band1','Band2','Band3'};


Y360p_roi_followup(Y360p_roi_followup.HeadMotion_3_0>0.2,:) = [];
Y360p_pca_followup(Y360p_pca_followup.HeadMotion_3_0>0.2,:) = [];
Y360p_band_followup(Y360p_band_followup.HeadMotion_3_0>0.2,:) = [];

[i2t3,i3t2] = idFinderNum(Y360p_roi.eid, Y360p_roi_followup.eid);
X360p_pca_delta = Y360p_pca_followup{i3t2,4:end} - Y360p_pca{i2t3,4:end};

% Split template

% Test Set + Verification Set (PCA)
% Leave one site out of the forecast
