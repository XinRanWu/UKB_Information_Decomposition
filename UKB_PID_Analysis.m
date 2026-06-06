%% Data loading
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Data.mat', 'ukb_RedSynYeo7n_2_0')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'Red360p_NC', 'Syn360p_NC')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'ukb_RedYeo7n_2_0', 'ukb_SynYeo7n_2_0')
load('/public/home/zhangjie/DataAnalysis/wxr_toolbox/toolbox_matlab/cbrewer2-master/cbrewer2/colorbrewer.mat', 'colorbrewer');
% UKB_Outcomes_2_0 = readtable('/public/home/zhangjie/UKB_Outcomes_2_0.csv');
% UKB_Lifestyle_0_0 = readtable('/public/home/zhangjie/UKB_Lifestyle_MRI_0_0.csv');
% UKB_WebQ_0_0 = readtable("/public/home/zhangjie/UKB_WebQ_0_0.txt",'TreatAsEmpty','NA');
% UKB_WebQ_0_0 = removevars(UKB_WebQ_0_0, 'Var1');


UKB_Lifestyle_0_0 = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Lifestyle_0_0.csv');
UKB_Lifestyle_2_0 = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Lifestyle_2_0.csv');
UKB_Outcomes_2_0 = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Outcomes_2_0.csv');
UKB_Outcomes_2_0.x20012_2_5 = str2double(UKB_Outcomes_2_0.x20012_2_5);
UKB_Outcomes_2_0(:, mean(ismissing(UKB_Outcomes_2_0),1) > 0.9) = [];
UKB_Outcomes_2_0 = removevars(UKB_Outcomes_2_0, 'x20006_2_0');

is_nan_mask = isnan(UKB_Outcomes_2_0.x4803_2_0);
UKB_Outcomes_2_0.x4803_2_0 = double(UKB_Outcomes_2_0.x4803_2_0 > 0);
UKB_Outcomes_2_0.x4803_2_0(is_nan_mask) = NaN;

UKB_Lifestyle_2_0(:, mean(ismissing(UKB_Lifestyle_2_0),1) > 0.9) = [];

UKB_Lifestyle_2_0 = removevars(UKB_Lifestyle_2_0, {'x845_2_0','x1767_2_0','x1677_2_0','x1787_2_0','x1777_2_0'});
UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, {'x20117_0_0','x1558_0_0','x1568_0_0','x1578_0_0','x1588_0_0',...
    'x1598_0_0','x1608_0_0','x4407_0_0','x4418_0_0','x4429_0_0','x4440_0_0','x4451_0_0','x1309_0_0','x1319_0_0',...
    'x1289_0_0','x1299_0_0','x1329_0_0','x1339_0_0','x1349_0_0','x1369_0_0','x1379_0_0','x1389_0_0','x1438_0_0',...
    'x1448_0_0','x1458_0_0','x1468_0_0','x884_0_0','x894_0_0','x904_0_0','x914_0_0','x1160_0_0','x1200_0_0',...
    'x20116_0_0','x1070_0_0','x1080_0_0','x709_0_0','x1031_0_0'});


UKB_Lifestyle_0_0 = innerjoin(UKB_Lifestyle_0_0,UKB_Lifestyle_2_0);
UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, {'x1468_2_0','x1448_2_0'});

Yeo7Names = {'Visual';'Somatomotor';'DorsalAttention';'VentralAttention';'Limbic';'Frontoparietal';'Default'};

UKB_Lifestyle_0_0.x1160_2_0_L = double(UKB_Lifestyle_0_0.x1160_2_0 >= 9);
UKB_Lifestyle_0_0.x1160_2_0_S = double(UKB_Lifestyle_0_0.x1160_2_0 < 6);
UKB_Lifestyle_0_0.x1160_2_0_S(isnan(UKB_Lifestyle_0_0.x1160_2_0)) = NaN;
UKB_Lifestyle_0_0.x1160_2_0_L(isnan(UKB_Lifestyle_0_0.x1160_2_0)) = NaN;
UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, 'x1160_2_0');
x709_2_0 = UKB_Lifestyle_0_0.x709_2_0;
UKB_Lifestyle_0_0.x709_2_0 = double(UKB_Lifestyle_0_0.x709_2_0 == 1);
UKB_Lifestyle_0_0.x709_2_0(isnan(x709_2_0)) = NaN;

% UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, {'x1289','x1299','x1329','x1339','x1349','x1369','x1379','x1389','x1438','x1448','x1458'});
% UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, {'x1568','x1578','x1588','x1598','x1608'});

% UKB_Lifestyle_0_0 = innerjoin(UKB_Lifestyle_0_0,UKB_WebQ_0_0);

load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP')

brain_data_360p = [ukb_RedSynYeo7n_2_0(:,1), array2table(Red360p_NC),array2table(Syn360p_NC)];
Red7n_NC = netMean1D(brain_data_360p{:,2:361},Yeo7MMP);
Syn7n_NC = netMean1D(brain_data_360p{:,362:721},Yeo7MMP);
brain_data_7n = [brain_data_360p(:,1),array2table(Red7n_NC),array2table(Syn7n_NC)];
brain_data_7n_matrix = ukb_RedSynYeo7n_2_0;

clear Red360p_NC Red7n_NC results_7n results_lifestyle_7nx results_pheno_adult results_pheno_older ...
    Syn360p_NC Syn7n_NC T T_ab ukb_SynYeo7n_2_0 ukb_RedYeo7n_2_0 ukb_RedSynYeo7n_2_0

%% Demographic Variables Asociation
UKB_Basic = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Basic.csv','TreatAsEmpty','NA');
UKB_BrainMRICov = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BrainMRICov.csv');
UKB_BrainMRICov.t1_R_HBF_3_0 = str2double(UKB_BrainMRICov.t1_R_HBF_3_0);

UKB_BrainMRICov(sum(isnan(UKB_BrainMRICov{:,2:end}),2)==26,:) = [];

UKB_BrainMRICov.Properties.VariableNames{22} = 'Vol_WB_TIV_2_0';
UKB_BrainMRICov.Properties.VariableNames{23} = 'Vol_WB_TIV_3_0';
UKB_BrainMRICov.Properties.VariableNames{24} = 'HBF_L_HBF_2_0';
UKB_BrainMRICov.Properties.VariableNames{26} = 'HBF_R_HBF_2_0';
UKB_BrainMRICov.Properties.VariableNames{25} = 'HBF_L_HBF_3_0';
UKB_BrainMRICov.Properties.VariableNames{27} = 'HBF_R_HBF_3_0';
UKB_BrainMRICov.Properties.VariableNames{14} = 'HeadMotion_2_0';
UKB_BrainMRICov.Properties.VariableNames{15} = 'HeadMotion_3_0';
UKB_BrainMRICov.Properties.VariableNames{8} = 'SNR_2_0';
UKB_BrainMRICov.Properties.VariableNames{9} = 'SNR_3_0';
UKB_BrainMRICov.Euler_Num_L_2_0 = 2 - 2 .* UKB_BrainMRICov.HBF_L_HBF_2_0;
UKB_BrainMRICov.Euler_Num_R_2_0 = 2 - 2 .* UKB_BrainMRICov.HBF_R_HBF_2_0;
UKB_BrainMRICov.Euler_Num_WB_2_0 = UKB_BrainMRICov.Euler_Num_R_2_0 + UKB_BrainMRICov.Euler_Num_L_2_0;

% 360 regions * age
UKB_BrainXY = innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; 

formula = ['AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


UKB_BrainXY.AgeAttend_2_0_Sq = zscore(UKB_BrainXY.AgeAttend_2_0).^2;
formula = ['AgeAttend_2_0+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0_Sq'};
results_age2_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


% 
% UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_360p);
% formula = ['Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + fmri_mot_all_mean_2_0 + ',...
%     'fmri_time_Z_disp_med_2_0 + fmri_tSNR_clean_2_0 + TDI_0_0 + EduYr_0_0 + ',...
%     'Race_1 + Race_2 + Race_3 + Race_4 + Vol_WB_TIV_2_0 + (1|Centre_0_0)'];
% 
% snr_threshold = mean(UKB_BrainXY.SNR_2_0, 'omitnan') - 3 * std(UKB_BrainXY.SNR_2_0, 'omitnan');

% UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0 > 0.2, :) = []; 
% UKB_BrainXY(UKB_BrainXY.SNR_2_0 < snr_threshold, :) = []; 
% YList = brain_data_360p.Properties.VariableNames(2:end)'; 
% XList = {'AgeAttend_2_0'};
% results_age_360p_Control = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
% 

formula = ['AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

% 7*7 network-modules
UKB_BrainXY = innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),brain_data_7n_matrix);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; 

UKB_BrainXY.AgeAttend_2_0_Sq = zscore(UKB_BrainXY.AgeAttend_2_0).^2;

formula = ['Sex_0_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0_Sq'};
results_age2_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

% 7 networks
UKB_BrainXY = innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; 

UKB_BrainXY.AgeAttend_2_0_Sq = zscore(UKB_BrainXY.AgeAttend_2_0).^2;

% 
% UKB_BrainXY_Control = UKB_BrainXY;
% UKB_BrainXY_Control(zscore(UKB_BrainXY_Control.fmri_tSNR_clean_2_0)>3,:) = [];
% formula = ['Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + fmri_mot_all_mean_2_0 + ',...
%     'fmri_time_Z_disp_med_2_0 + fmri_tSNR_clean_2_0 + ',...
%     'Race_1 + Race_2 + Race_3 + Race_4 + Vol_WB_TIV_2_0 + (1|Centre_0_0)'];
% XList = {'AgeAttend_2_0'};
% results_age_7n_Control = LMM_with_EffectSizes_Full(UKB_BrainXY_Control,formula,XList,YList);


formula = ['Sex_0_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0_Sq'};
results_age2_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);



formula = ['AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


results_demo_360p = [results_age_360p;results_age2_360p;results_sex_360p;results_bmi_360p;results_tiv_360p];
results_demo_7n_matrix = [results_age_7n_matrix;results_age2_7n_matrix;results_sex_7n_matrix;results_bmi_7n_matrix;results_tiv_7n_matrix];
results_demo_7n = [results_age_7n;results_age2_7n;results_sex_7n;results_bmi_7n;results_tiv_7n];


lmm = fitlme(UKB_BrainXY,'Red7n_NC7~Sex_0_0+AgeAttend_2_0+BMI_2_0+AgeAttend_2_0^2+HeadMotion_2_0+Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)')
fe_table = lmm.Coefficients;
beta1 = fe_table.Estimate(strcmp(fe_table.Name, 'AgeAttend_2_0'));
beta2 = fe_table.Estimate(strcmp(fe_table.Name, 'AgeAttend_2_0^2'));
z_inflection = -beta1 / (2 * beta2);
age_mean = mean(UKB_BrainXY.AgeAttend_2_0, 'omitnan');
age_std  = std(UKB_BrainXY.AgeAttend_2_0, 'omitnan');
age_inflection_years = (z_inflection * age_std) + age_mean;

fprintf('--- Analysis of DMN Redundant Quadratic Term Inflection Points ---\n');
fprintf('First-order term Beta1: %.4f\n', beta1);
fprintf('Second-order term Beta2: %.4f\n', beta2);
fprintf('Z-score Inflection Point: %.4f\n', z_inflection); fprintf("Age of inflection in the real world: %.2f years\n", age_inflection_years);
min_age = min(UKB_BrainXY.AgeAttend_2_0);
max_age = max(UKB_BrainXY.AgeAttend_2_0);



%% Phenotype Association
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];

% 360 regions
UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_Outcomes_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = {'x20016_2_0','x4526_2_0','x2178_2_0'}';
% XList = {'x20016_2_0';'x6373_2_0';'x20018_2_0';'x20023_2_0';'x23324_2_0';'x21004_2_0';'x6348_2_0';'x20197_2_0'};
results_pheno_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_pheno_adult_360p = LMM_with_EffectSizes_Full(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0<65,:),formula,XList,YList);
results_pheno_older_360p = LMM_with_EffectSizes_Full(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0>=65,:),formula,XList,YList);

V = '21004';
Xt = results_pheno_360p.t(contains(results_pheno_360p.xname,V));
Xp = results_pheno_360p.p(contains(results_pheno_360p.xname,V));
X = Xt.*(Xp<(0.05/360));
[a, cb] = plot_cortical(parcel_to_surface(X([181:360,1:180]+360),'glasser_360_fsa5'),'color_range',[-7,7],'layout','grid');


% 7*7 network-modules
UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_Outcomes_2_0),brain_data_7n_matrix);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
results_pheno_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_pheno_adult_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0<65,:),formula,XList,YList);
results_pheno_older_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0>=65,:),formula,XList,YList);

% 7 networks

UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
results_pheno_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_pheno_adult_7n = LMM_with_EffectSizes_Full(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0<65,:),formula,XList,YList);
results_pheno_older_7n = LMM_with_EffectSizes_Full(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0>=65,:),formula,XList,YList);


results_pheno_7n(isnan(results_pheno_7n.t),:) = [];
results_pheno_adult_7n(isnan(results_pheno_adult_7n.t),:) = [];
results_pheno_older_7n(isnan(results_pheno_older_7n.t),:) = [];


results_pheno_7n_matrix(isnan(results_pheno_7n_matrix.t),:) = [];
results_pheno_adult_7n_matrix(isnan(results_pheno_adult_7n_matrix.t),:) = [];
results_pheno_older_7n_matrix(isnan(results_pheno_older_7n_matrix.t),:) = [];

res = results_pheno_7n(mafdr(results_pheno_7n.p,'BHFDR',true)<0.05,:);
resA = results_pheno_adult_7n(mafdr(results_pheno_adult_7n.p,'BHFDR',true)<0.05,:);
resO = results_pheno_older_7n(mafdr(results_pheno_older_7n.p,'BHFDR',true)<0.05,:);



% 7 networks

formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+AgeAttend_2_0^2+HeadMotion_2_0+x2443_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic,UKB_BrainMRICov),UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
XList(contains(XList,'x2443')) = [];
results_pheno_7n_Control = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

results_pheno_7n_Control(isnan(results_pheno_7n_Control.t),:) = [];


results_pheno_7n(contains(results_pheno_7n.xname,'x2443'),:) = [];




resRed = results_pheno_7n(contains(results_pheno_7n.yname,'Red'),:);
resSyn = results_pheno_7n(contains(results_pheno_7n.yname,'Syn'),:);
resRedC = results_pheno_7n_Control(contains(results_pheno_7n_Control.yname,'Red'),:);
resSynC = results_pheno_7n_Control(contains(results_pheno_7n_Control.yname,'Syn'),:);
t1 = resSyn.t;
t2 = resSynC.t;
xname_list = resRed.xname;
yname_list = resRed.yname;

id_keys = UKB_OutcomesList.FieldID;
if isnumeric(id_keys)
    id_keys_cell = cellfun(@num2str, num2cell(id_keys), 'UniformOutput', false);
else
    id_keys_cell = cellfun(@num2str, id_keys, 'UniformOutput', false);
end
abbr_values = UKB_OutcomesList.Abbreviation;
dict_map = containers.Map(id_keys_cell, abbr_values);


N_total = length(xname_list);
Label = cell(N_total, 1);
net_colors = zeros(N_total, 3); 

for i = 1:N_total

    raw_x = xname_list{i};
    tokens = regexp(raw_x, 'x(\d+)_', 'tokens');
    if ~isempty(tokens)
        field_id_str = tokens{1}{1};
    else
        field_id_str = '';
    end

    if isKey(dict_map, field_id_str)
        short_name = dict_map(field_id_str);
    else
        short_name = ['Field_' field_id_str]; 
    end
    raw_y = yname_list{i};
    net_tokens = regexp(raw_y, 'NC(\d+)', 'tokens');
    if ~isempty(net_tokens)
        net_num = str2double(net_tokens{1}{1});
    else
        net_num = 1; 
    end
 
    Label{i} = [short_name];
    net_colors(i, :) = Yeo7Color(net_num, :);
end

figure('Color', 'w');
scatter(t1, t2, 60, net_colors, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;

all_data = [t1(:); t2(:)];
global_min = min(all_data) - 0.5;
global_max = max(all_data) + 0.5;
plot([global_min, global_max], [global_min, global_max], 'k--', 'LineWidth', 1.5);
xlim([global_min, global_max]);
ylim([global_min, global_max]);

offsets = t2 - t1;
abs_offsets = abs(offsets);
[~, sort_idx] = sort(abs_offsets, 'descend');
top_N = 20; 
labeled_indices = sort_idx(1:top_N);


for k = 1:top_N
    idx = labeled_indices(k);
    
    this_color = net_colors(idx, :);

    text_color = this_color;
    if sum(text_color) > 2.2 
        text_color = text_color * 0.7; 
    end

    text(t1(idx) - 0.15, t2(idx) + 0.15, Label{idx}, ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', text_color, ... 
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom');
        
    plot(t1(idx), t2(idx), 'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', text_color);
end

xlabel('T-values (Without Controlling Diabetes)', 'FontSize', 12);
ylabel('T-values (With Controlling Diabetes)', 'FontSize', 12);
title(['Phenotypic Consistency (r = ' num2str(corr(t1(:), t2(:)), '%.3f') ')'], 'FontSize', 13);
grid on; 
axis square;


%% Envrionmental & Lifstyle Association

formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];




% 7 networks
UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),UKB_Lifestyle_0_0),brain_data_7n);
% UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = UKB_Lifestyle_0_0.Properties.VariableNames(2:end)';
results_lifestyle_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_lifestyle_7n(isnan(results_lifestyle_7n.t),:) = [];
% writetable(results_lifestyle_7n,'modal_lifestyle_7n_matrix.txt','Delimiter' ,' ')

res = results_lifestyle_7n(mafdr(results_lifestyle_7n.p,'BHFDR',true)<0.05,:);

XList = unique(results_lifestyle_7n.xname(mafdr(results_lifestyle_7n.p,'BHFDR',true)<0.001));

V = 'x1458';
Xt = results_lifestyle_7n_matrix.t(contains(results_lifestyle_7n_matrix.xname,V));
Xp = results_lifestyle_7n_matrix.p(contains(results_lifestyle_7n_matrix.xname,V));
X0 = Xt.*(Xp<(0.05/28));
X01 = X0(1:28);
X02 = X0(30:57)
X1 = icatb_vec2mat(X01(8:28));X1(find(eye(7))) = X01(1:7);
plot_matrix_Yeo7n(X1,8.5)
X2 = icatb_vec2mat(X02(8:28));X2(find(eye(7))) = X01(1:7);
plot_matrix_Yeo7n(X2,8.5)

for i = 4
    V = XList{i};
    Xt = results_lifestyle_7n_matrix.t(contains(results_lifestyle_7n_matrix.xname,V));
    Xp = results_lifestyle_7n_matrix.p(contains(results_lifestyle_7n_matrix.xname,V));
    X0 = Xt.*(Xp<(0.05/28));
    X01 = X0(1:28);
    X02 = X0(30:57);
    X1 = icatb_vec2mat(X01(8:28));X1(find(eye(7))) = X01(1:7);
    plot_matrix_Yeo7n(X1,10.5)
    X2 = icatb_vec2mat(X02(8:28));X2(find(eye(7))) = X02(1:7);
    plot_matrix_Yeo7n(X2,10.5)
    
end


% 7*7 network-modules
UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_Lifestyle_0_0),brain_data_7n_matrix);
% UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
XList = UKB_Lifestyle_0_0.Properties.VariableNames(2:end)';
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
results_lifestyle_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
% writetable(results_lifestyle_7n_matrix,'modal_lifestyle_7n_matrix.txt','Delimiter' ,' ')

% 360 regions
formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),UKB_Lifestyle_0_0),brain_data_360p);
% UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)';
results_lifestyle_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
% writetable(results_lifestyle_360p,'modal_lifestyle_7n.txt','Delimiter' ,' ')

i = 11;
V = XList{i};
Xt = results_lifestyle_360p.t(contains(results_lifestyle_360p.xname,V));
Xp = results_lifestyle_360p.p(contains(results_lifestyle_360p.xname,V));
X0 = Xt.*(Xp<(0.05/360));
X01 = X0(1:360);
X02 = X0(361:end);
[a, cb] = plot_cortical(parcel_to_surface(X01([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-8.5,8.5],'layout','grid');

[a, cb] = plot_cortical(parcel_to_surface(X02([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-8.5,8.5],'layout','grid');

% Interaction
% formula_interaction = ['AgeAttend_2_0 * BMI_2_0 + Sex_0_0 + HeadMotion_2_0 + ',...
%                        'SNR_2_0 + Vol_WB_TIV_2_0 + (1|Centre_0_0)'];

% UKB_BrainXY.SA2Y = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2Y(UKB_BrainXY.AgeAttend_2_0 < 65) = 0;
% UKB_BrainXY.SA2Y((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 7)) = 1;
% 
% UKB_BrainXY.SA2NA = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 <= 6)) = 0;
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 6)) = 1;
% 
% XList = {'SA2Y','SA2NA'}';
% results_pheno_360p_SA = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
% 
% V = 'SA2Y';
% X = results_pheno_360p_SA.t(strcmp(results_pheno_360p_SA.xname,V)).*...
%     (results_pheno_360p_SA.p(strcmp(results_pheno_360p_SA.xname,V))<(0.05/360));
% 
% [a, cb] = plot_cortical(parcel_to_surface(X([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-6,6]);


% UKB_BrainXY.SA2Y = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2Y(UKB_BrainXY.AgeAttend_2_0 < 65) = 0;
% UKB_BrainXY.SA2Y((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 7)) = 1;
% 
% UKB_BrainXY.SA2NA = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 <= 6)) = 0;
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 6)) = 1;
% 
% XList = {'SA2Y','SA2NA'}';
% results_pheno_7n_matrix_SA = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
% UKB_BrainXY.SA2Y = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2Y(UKB_BrainXY.AgeAttend_2_0 < 65) = 0;
% UKB_BrainXY.SA2Y((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 7)) = 1;
% 
% UKB_BrainXY.SA2NA = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 <= 6)) = 0;
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 6)) = 1;
% 
% XList = {'SA2Y','SA2NA'}';
% results_pheno_7n_SA = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);


%% Visualize Heatmap
UKB_OutcomesList = readtable("/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Outcomes_List.csv");
UKB_LifestyleList = readtable("/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Lifestyle_List.csv");
UKB_LifestyleList = [UKB_LifestyleList(1,:);UKB_LifestyleList];
UKB_LifestyleList.FieldID(1) = 189;
UKB_LifestyleList.Description{1} = 'Townsend deprivation index';
UKB_LifestyleList.Description{22} = 'Avoided activities or situations for previous stress in past month';

res_clean = results_pheno_7n;
refInfo = UKB_OutcomesList;

res_clean = results_lifestyle_7n;
refInfo = UKB_LifestyleList;



res_clean = [results_pheno_7n;results_lifestyle_7n];
refInfo = [UKB_LifestyleList;UKB_OutcomesList];

fdr_p = mafdr(res_clean.p, 'BHFDR', true);
res_clean.is_significant = fdr_p < 0.05;

t_table = unstack(res_clean(:, {'xname', 'yname', 't'}), 't', 'yname');
sig_table = unstack(res_clean(:, {'xname', 'yname', 'is_significant'}), 'is_significant', 'yname');

% Extract unique environment factors (X-axis raw keys: e.g., 'x738_0_0')
x_labels_raw = t_table.xname; 
y_labels = t_table.Properties.VariableNames(2:end); 

% Convert tables to numeric matrices
t_matrix = table2array(t_table(:, 2:end));
sig_matrix = table2array(sig_table(:, 2:end));
sig_matrix(isnan(sig_matrix)) = 0; % Clear NaN padding from unstack

x_labels_new = cell(size(x_labels_raw));
for i = 1:length(x_labels_raw)
    % Use regex to extract the main FieldID number between 'x' and the first '_'
    % For example, 'x738_0_0' will yield tokens{1} = '738'
    tokens = regexp(x_labels_raw{i}, '^x(\d+)_', 'tokens');
    
    if ~isempty(tokens)
        field_id_num = str2double(tokens{1}{1});
        
        % Look up the matching row in refInfo (UKB_OutcomesList)
        match_idx = find(refInfo.FieldID == field_id_num, 1);
        
        if ~isempty(match_idx)
            % If matched, extract Description (and optionally append Instance info if desired)
            x_labels_new{i} = refInfo.Abbreviation{match_idx};
        else
            % Fallback to raw label if no match found in refInfo
            x_labels_new{i} = x_labels_raw{i};
        end
    else
        x_labels_new{i} = x_labels_raw{i};
    end
end
y_labels_new = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN', ...
                'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN'};
t_matrix_final = t_matrix';
sig_matrix_final = sig_matrix';

figure('Color', 'w', 'Position', [100, 100, 1050, 650]); % Wider layout for Descriptions
h_im = imagesc(t_matrix_final, [-9.2,9.2]);
hold on;
% Apply user-defined Blue-Red colormap
colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)));
axis image; 
[num_rows, num_cols] = size(t_matrix_final);

% Overlay custom black grid lines
for r = 0:num_rows
    plot([0.5, num_cols+0.5], [r+0.5, r+0.5], 'k-', 'LineWidth', 0.8);
end
for c = 0:num_cols
    plot([c+0.5, c+0.5], [0.5, num_rows+0.5], 'k-', 'LineWidth', 0.8);
end

% Overlay significance markers (Asterisks)
for r = 1:num_rows
    for c = 1:num_cols
        if sig_matrix_final(r, c) == 1
            text(c, r, '*', 'Color', 'k', 'FontSize', 16, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
        end
    end
end

% Set the beautifully converted Descriptions on X-axis and Networks on Y-axis
set(gca, 'XTick', 1:num_cols, 'XTickLabel', x_labels_new, 'TickLabelInterpreter', 'none');
set(gca, 'YTick', 1:num_rows, 'YTickLabel', y_labels_new, 'TickLabelInterpreter', 'none');
% Rotate X-axis text by 45 or 90 degrees (45 recommended since descriptions are longer text)
xtickangle(90);
box on;
set(gca, 'Layer', 'top', 'LineWidth', 1);
% Configure the symmetric colorbar
c_bar = colorbar;
c_bar.Label.String = 't-statistic';
c_bar.Label.FontSize = 11;
max_abs_t = max(abs(t_matrix_final(:)));
caxis([-max_abs_t, max_abs_t]); 

title('Association Map (FDR Correction)', ...
      'FontSize', 13, 'FontWeight', 'bold');
  
  
  
  
res_clean_all = [results_pheno_7n; results_lifestyle_7n];
refInfo = [UKB_LifestyleList; UKB_OutcomesList];
res_clean = res_clean_all(contains(res_clean_all.yname, 'Red'), :);

fdr_p = mafdr(res_clean.p, 'BHFDR', true);
res_clean.is_significant = fdr_p < 0.05;
res_clean.signed_logp = -log10(res_clean.p) .* sign(res_clean.t);

x_labels_raw = res_clean.xname;
num_rows_total = height(res_clean);
categories = cell(num_rows_total, 1);
abbreviations = cell(num_rows_total, 1);

for i = 1:num_rows_total
    tokens = regexp(x_labels_raw{i}, '^x(\d+)_', 'tokens');
    if ~isempty(tokens)
        field_id_num = str2double(tokens{1}{1});
        match_idx = find(refInfo.FieldID == field_id_num, 1);
        if ~isempty(match_idx)
            categories{i} = char(refInfo.Category_sub{match_idx});
            abbreviations{i} = char(refInfo.Abbreviation{match_idx});
        else
            categories{i} = 'Unknown';
            abbreviations{i} = char(x_labels_raw{i});
        end
    else
        categories{i} = 'Unknown';
        abbreviations{i} = char(x_labels_raw{i});
    end
end

res_clean.Category = categories;
res_clean.Abbreviation = abbreviations;
res_clean = sortrows(res_clean, {'Category', 'Abbreviation'});

[unique_phenos, ~] = unique(res_clean.Abbreviation, 'stable');
num_phenos = length(unique_phenos);

[~, unique_pheno_idx] = unique(res_clean.Abbreviation, 'stable');
pheno_categories = res_clean.Category(unique_pheno_idx);

res_clean.X_coord = zeros(height(res_clean), 1);
for i = 1:height(res_clean)
    current_name = res_clean.Abbreviation{i};
    idx = find(strcmp(unique_phenos, current_name), 1);
    if ~isempty(idx)
        res_clean.X_coord(i) = idx;
    end
end

figure('Color', 'w', 'Position', [100, 100, 1200, 750]);
hold on;

networks = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN'};
Yeo7Color = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78] ./ 255;

plot([0.5, num_phenos + 0.5], [0, 0], 'k-', 'LineWidth', 1.2);

h_scatter = zeros(length(networks), 1);
for n = 1:length(networks)
    net_mask = contains(res_clean.yname, ['NC', num2str(n)]);
    net_data = res_clean(net_mask, :);
    if ~isempty(net_data)
        h_scatter(n) = scatter(net_data.X_coord, net_data.signed_logp, 65, ...
            'MarkerFaceColor', Yeo7Color(n, :), 'MarkerEdgeColor', [0.15, 0.15, 0.15], ...
            'LineWidth', 0.5, 'MarkerFaceAlpha', 0.85);
    end
end

sig_all_data = res_clean(res_clean.is_significant == 1, :);
if ~isempty(sig_all_data)
    [sig_unique_phenos, ~] = unique(sig_all_data.Abbreviation, 'stable');
    num_sig_phenos = length(sig_unique_phenos);
    sig_data_list = cell(num_sig_phenos, 1);
    
    for p = 1:num_sig_phenos
        p_mask = strcmp(sig_all_data.Abbreviation, sig_unique_phenos{p});
        p_rows = sig_all_data(p_mask, :);
        [~, max_idx] = max(abs(p_rows.signed_logp));
        sig_data_list{p} = p_rows(max_idx, :);
    end
    
    sig_data = vertcat(sig_data_list{:});
    txt_x = sig_data.X_coord;
    txt_y = sig_data.signed_logp;
    num_sig = length(txt_x);
    
    for i = 1:num_sig
        txt_y(i) = txt_y(i) + sign(txt_y(i)) * 0.4;
    end
    
    iterations = 150; alpha = 0.25; min_dist_x = 4.5; min_dist_y = 1.1;
    for iter = 1:iterations
        for i = 1:num_sig
            for j = 1:num_sig
                if i == j, continue; end
                dx = txt_x(i) - txt_x(j); dy = txt_y(i) - txt_y(j);
                if abs(dx) < min_dist_x && abs(dy) < min_dist_y
                    if dx == 0, dx = randn() * 0.1; end
                    if dy == 0, dy = randn() * 0.1; end
                    txt_x(i) = txt_x(i) + alpha * (sign(dx) * min_dist_x - dx) * 0.3;
                    txt_y(i) = txt_y(i) + alpha * (sign(dy) * min_dist_y - dy) * 0.6;
                end
            end
            txt_x(i) = max(1, min(num_phenos, txt_x(i)));
        end
    end
    
    for i = 1:num_sig
        plot([sig_data.X_coord(i), txt_x(i)], [sig_data.signed_logp(i), txt_y(i)], '-', 'Color', [0.65, 0.65, 0.65], 'LineWidth', 0.6);
        text(txt_x(i), txt_y(i), sprintf('%s', sig_data.Abbreviation{i}), 'FontSize', 8, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'BackgroundColor', [1 1 1 0.85], 'EdgeColor', [0.75 0.75 0.75], 'Margin', 1.5);
    end
end

y_limits = ylim;
y_range = y_limits(2) - y_limits(1);
new_y_min = y_limits(1) - y_range * 0.25; 
ylim([new_y_min, y_limits(2)]);

[unique_cats, cat_start_idx] = unique(pheno_categories, 'stable');
cat_end_idx = [cat_start_idx(2:end) - 1; num_phenos];
cat_colors = [0.8941 0.1020 0.1098; 0.2157 0.4941 0.7216; 0.3020 0.6863 0.2902; 0.5961 0.3059 0.6392; 1.0000 0.4980 0.0000; 1.0000 0.8200 0.1500; 0.6510 0.3373 0.1569];
num_cats = length(unique_cats);
if num_cats > size(cat_colors, 1), cat_colors = repmat(cat_colors, ceil(num_cats / size(cat_colors, 1)), 1); end

bar_y_bottom = new_y_min + y_range * 0.18;
bar_y_height = y_range * 0.04;

for k = 1:num_cats
    x_start = cat_start_idx(k) - 0.5;
    x_width = cat_end_idx(k) - cat_start_idx(k) + 1;
    
    rectangle('Position', [x_start, bar_y_bottom, x_width, bar_y_height], 'FaceColor', cat_colors(k, :), 'EdgeColor', [0.2, 0.2, 0.2], 'LineWidth', 0.8);
    
    x_center = x_start + x_width / 2;
    text(x_center, bar_y_bottom - y_range * 0.02, unique_cats{k}, 'FontSize', 9.5, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Rotation', 90, 'Interpreter', 'none');
     
    if k > 1
        plot([x_start, x_start], [bar_y_bottom + bar_y_height, y_limits(2)], '--', 'Color', [0.82, 0.82, 0.82], 'LineWidth', 0.8);
    end
end

set(gca, 'XTick', [], 'XTickLabel', {}, 'XColor', 'none'); 
xlim([0.5, num_phenos + 0.5]);
ylabel('Signed -log_{10}(p-value)', 'FontSize', 12, 'FontWeight', 'bold');
title('Signed Manhattan Plot across Networks & Categories', 'FontSize', 14, 'FontWeight', 'bold');

lgd = legend(h_scatter(h_scatter > 0), networks(h_scatter > 0), 'Location', 'northeastoutside');
lgd.Title.String = 'Networks'; lgd.Box = 'on';

grid on; set(gca, 'GridColor', [0.93, 0.93, 0.93], 'GridAlpha', 0.6); box on;
set(gca, 'Layer', 'top', 'LineWidth', 1.2, 'TickDir', 'out');




res_clean_all = [results_pheno_7n; results_lifestyle_7n];
refInfo = [UKB_LifestyleList; UKB_OutcomesList];

% =========================================================================
% 1. Data Filtering & Significance Correction
% =========================================================================
fdr_p = mafdr(res_clean_all.p, 'BHFDR', true);
res_clean_all.is_significant = fdr_p < 0.05;

x_labels_raw = res_clean_all.xname;
num_rows_total = height(res_clean_all);
abbreviations = cell(num_rows_total, 1);

for i = 1:num_rows_total
    tokens = regexp(x_labels_raw{i}, '^x(\d+)_', 'tokens');
    if ~isempty(tokens)
        field_id_num = str2double(tokens{1}{1});
        match_idx = find(refInfo.FieldID == field_id_num, 1);
        if ~isempty(match_idx)
            abbreviations{i} = char(refInfo.Abbreviation{match_idx});
        else
            abbreviations{i} = char(x_labels_raw{i});
        end
    else
        abbreviations{i} = char(x_labels_raw{i});
    end
end
res_clean_all.Abbreviation = abbreviations;

% Split into Red and Syn tables (Significant only)
red_table = res_clean_all(contains(res_clean_all.yname, 'Red') & res_clean_all.is_significant == 1, :);
syn_table = res_clean_all(contains(res_clean_all.yname, 'Syn') & res_clean_all.is_significant == 1, :);

left_phenos = unique(red_table.Abbreviation, 'stable');
right_phenos = unique(syn_table.Abbreviation, 'stable');

num_left = length(left_phenos);
num_right = length(right_phenos);
num_center = 7;

% =========================================================================
% 2. Coordinate Assignment (X, Y Framework)
% =========================================================================
figure('Color', 'w', 'Position', [50, 50, 1300, 800]);
hold on; axis off;

% Center Column (7 Networks) Positions
x_center = 0.5;
y_center_nodes = linspace(0.1, 0.9, num_center)';

% Left Column (Red Phenotypes) Positions
x_left = 0.1;
y_left_nodes = linspace(0.05, 0.95, num_left)';

% Right Column (Syn Phenotypes) Positions
x_right = 0.9;
y_right_nodes = linspace(0.05, 0.95, num_right)';

% Maps to instantly retrieve Y positions
left_map = containers.Map(left_phenos, y_left_nodes);
right_map = containers.Map(right_phenos, y_right_nodes);

% =========================================================================
% 3. Geometry Styles & Colors
% =========================================================================
networks = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN'};
Yeo7Color = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78] ./ 255;

t_max = max(abs(res_clean_all.t(res_clean_all.is_significant == 1)));

% Bezier curve helper for elegant, smooth stream lines
draw_curve = @(x1, y1, x2, y2, col, w) plot(...
    kron((1-linspace(0,1,50)).^2, x1) + kron(2*(1-linspace(0,1,50)).*linspace(0,1,50), (x1+x2)/2) + kron(linspace(0,1,50).^2, x2), ...
    kron((1-linspace(0,1,50)).^2, y1) + kron(2*(1-linspace(0,1,50)).*linspace(0,1,50), y1) + kron(linspace(0,1,50).^2, y2), ...
    '-', 'Color', [col, 0.45], 'LineWidth', w);

% =========================================================================
% 4. Draw Edges (Left to Center & Center to Right)
% =========================================================================
% Left Edges (Red)
for i = 1:height(red_table)
    p_name = red_table.Abbreviation{i};
    y_p = left_map(p_name);
    
    for n = 1:num_center
        if contains(red_table.yname{i}, ['NC', num2str(n)])
            y_n = y_center_nodes(n);
            line_w = 0.5 + 3.5 * (abs(red_table.t(i)) / t_max); % Scale width by t-stat
            draw_curve(x_left, y_p, x_center, y_n, Yeo7Color(n,:), line_w);
        end
    end
end

% Right Edges (Syn)
for i = 1:height(syn_table)
    p_name = syn_table.Abbreviation{i};
    y_p = right_map(p_name);
    
    for n = 1:num_center
        if contains(syn_table.yname{i}, ['NC', num2str(n)])
            y_n = y_center_nodes(n);
            line_w = 0.5 + 3.5 * (abs(syn_table.t(i)) / t_max);
            draw_curve(x_center, y_n, x_right, y_p, Yeo7Color(n,:), line_w);
        end
    end
end

% =========================================================================
% 5. Draw Nodes & Text Labels
% =========================================================================
% Left Labels (Right-aligned to face the center)
for i = 1:num_left
    plot(x_left, y_left_nodes(i), 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
    text(x_left - 0.01, y_left_nodes(i), left_phenos{i}, 'FontSize', 8.5, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end

% Right Labels (Left-aligned to face the center)
for i = 1:num_right
    plot(x_right, y_right_nodes(i), 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
    text(x_right + 0.01, y_right_nodes(i), right_phenos{i}, 'FontSize', 8.5, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end

% Center Network Badges
for n = 1:num_center
    scatter(x_center, y_center_nodes(n), 350, 'MarkerFaceColor', Yeo7Color(n,:), ...
            'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);
    text(x_center, y_center_nodes(n), networks{n}, 'FontSize', 9, 'FontWeight', 'bold', ...
         'Color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% =========================================================================
% 6. Outer Titles
% =========================================================================
text(x_left, 0.98, 'Red7nNC Associated Phenotypes', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
text(x_right, 0.98, 'Syn7nNC Associated Phenotypes', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
text(x_center, 0.98, 'Brain Networks', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

xlim([x_left - 0.2, x_right + 0.2]);
ylim([0, 1]);




res_clean_all = [results_pheno_7n; results_lifestyle_7n];
refInfo = [UKB_LifestyleList; UKB_OutcomesList];

% =========================================================================
% 1. Data Filtering & Categorization
% =========================================================================
fdr_p = mafdr(res_clean_all.p, 'BHFDR', true);
res_clean_all.is_significant = fdr_p < 0.001;

x_labels_raw = res_clean_all.xname;
num_rows_total = height(res_clean_all);
categories = cell(num_rows_total, 1);
abbreviations = cell(num_rows_total, 1);

for i = 1:num_rows_total
    tokens = regexp(x_labels_raw{i}, '^x(\d+)_', 'tokens');
    if ~isempty(tokens)
        field_id_num = str2double(tokens{1}{1});
        match_idx = find(refInfo.FieldID == field_id_num, 1);
        if ~isempty(match_idx)
            categories{i} = char(refInfo.Category_sub{match_idx});
            abbreviations{i} = char(refInfo.Abbreviation{match_idx});
        else
            categories{i} = 'Unknown';
            abbreviations{i} = char(x_labels_raw{i});
        end
    else
        categories{i} = 'Unknown';
        abbreviations{i} = char(x_labels_raw{i});
    end
end
res_clean_all.Category = categories;
res_clean_all.Abbreviation = abbreviations;

% Split and keep significant associations
red_table = res_clean_all(contains(res_clean_all.yname, 'Red') & res_clean_all.is_significant == 1, :);
syn_table = res_clean_all(contains(res_clean_all.yname, 'Syn') & res_clean_all.is_significant == 1, :);

% OPTIMIZED ORDER: Sort tables by Category first, then by Abbreviation
red_table = sortrows(red_table, {'Category', 'Abbreviation'});
syn_table = sortrows(syn_table, {'Category', 'Abbreviation'});

% Extract unique phenotypes while preserving the new categorized order
left_phenos = unique(red_table.Abbreviation, 'stable');
right_phenos = unique(syn_table.Abbreviation, 'stable');

num_left = length(left_phenos);
num_right = length(right_phenos);
num_center = 7;

% =========================================================================
% 2. Coordinate Assignment
% =========================================================================
figure('Color', 'w', 'Position', [50, 50, 1300, 850]);
hold on; axis off;

x_left = 0.1;   x_center = 0.5;   x_right = 0.9;
y_center_nodes = linspace(0.15, 0.85, num_center)';
y_left_nodes = linspace(0.05, 0.95, num_left)';
y_right_nodes = linspace(0.05, 0.95, num_right)';

left_map = containers.Map(left_phenos, y_left_nodes);
right_map = containers.Map(right_phenos, y_right_nodes);

% =========================================================================
% 3. Geometry Styles & Smooth Bezier Curves
% =========================================================================
networks = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN'};
Yeo7Color = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78] ./ 255;
t_max = max(abs(res_clean_all.t(res_clean_all.is_significant == 1)));

draw_curve = @(x1, y1, x2, y2, col, w) plot(...
    kron((1-linspace(0,1,50)).^2, x1) + kron(2*(1-linspace(0,1,50)).*linspace(0,1,50), (x1+x2)/2) + kron(linspace(0,1,50).^2, x2), ...
    kron((1-linspace(0,1,50)).^2, y1) + kron(2*(1-linspace(0,1,50)).*linspace(0,1,50), y1) + kron(linspace(0,1,50).^2, y2), ...
    '-', 'Color', [col, 0.45], 'LineWidth', w);

% =========================================================================
% 4. Plot Edges (Left to Center & Center to Right)
% =========================================================================
for i = 1:height(red_table)
    p_name = red_table.Abbreviation{i};
    y_p = left_map(p_name);
    for n = 1:num_center
        if contains(red_table.yname{i}, ['NC', num2str(n)])
            line_w = 0.5 + 3.5 * (abs(red_table.t(i)) / t_max);
            draw_curve(x_left, y_p, x_center, y_center_nodes(n), Yeo7Color(n,:), line_w);
        end
    end
end

for i = 1:height(syn_table)
    p_name = syn_table.Abbreviation{i};
    y_p = right_map(p_name);
    for n = 1:num_center
        if contains(syn_table.yname{i}, ['NC', num2str(n)])
            line_w = 0.5 + 3.5 * (abs(syn_table.t(i)) / t_max);
            draw_curve(x_center, y_center_nodes(n), x_right, y_p, Yeo7Color(n,:), line_w);
        end
    end
end

% =========================================================================
% 5. Render Nodes & Text
% =========================================================================
for i = 1:num_left
    plot(x_left, y_left_nodes(i), 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
    text(x_left - 0.01, y_left_nodes(i), left_phenos{i}, 'FontSize', 8.5, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end

for i = 1:num_right
    plot(x_right, y_right_nodes(i), 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
    text(x_right + 0.01, y_right_nodes(i), right_phenos{i}, 'FontSize', 8.5, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end

for n = 1:num_center
    scatter(x_center, y_center_nodes(n), 400, 'MarkerFaceColor', Yeo7Color(n,:), 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 1.2);
    text(x_center, y_center_nodes(n), networks{n}, 'FontSize', 9.5, 'FontWeight', 'bold', 'Color', 'w', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

text(x_left, 0.98, 'Red7nNC Associated Phenotypes', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
text(x_right, 0.98, 'Syn7nNC Associated Phenotypes', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
text(x_center, 0.98, 'Brain Networks', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

xlim([x_left - 0.2, x_right + 0.2]); ylim([0, 1]);










res_clean_all = [results_pheno_7n; results_lifestyle_7n];
refInfo = [UKB_LifestyleList; UKB_OutcomesList];

fdr_p = mafdr(res_clean_all.p, 'BHFDR', true);
res_clean_all.is_significant = fdr_p < 0.05;

x_labels_raw = res_clean_all.xname;
num_rows_total = height(res_clean_all);
categories = cell(num_rows_total, 1);
abbreviations = cell(num_rows_total, 1);

for i = 1:num_rows_total
    tokens = regexp(x_labels_raw{i}, '^x(\d+)_', 'tokens');
    if ~isempty(tokens)
        field_id_num = str2double(tokens{1}{1});
        match_idx = find(refInfo.FieldID == field_id_num, 1);
        if ~isempty(match_idx)
            categories{i} = char(refInfo.Category_sub{match_idx});
            abbreviations{i} = char(refInfo.Abbreviation{match_idx});
        else
            categories{i} = 'Unknown';
            abbreviations{i} = char(x_labels_raw{i});
        end
    else
        categories{i} = 'Unknown';
        abbreviations{i} = char(x_labels_raw{i});
    end
end
res_clean_all.Category = categories;
res_clean_all.Abbreviation = abbreviations;

red_table = res_clean_all(contains(res_clean_all.yname, 'Red') & res_clean_all.is_significant == 1, :);
syn_table = res_clean_all(contains(res_clean_all.yname, 'Syn') & res_clean_all.is_significant == 1, :);

red_table = sortrows(red_table, {'Category', 'Abbreviation'});
syn_table = sortrows(syn_table, {'Category', 'Abbreviation'});

left_phenos = unique(red_table.Abbreviation, 'stable');
right_phenos = unique(syn_table.Abbreviation, 'stable');

num_left = length(left_phenos);
num_right = length(right_phenos);
num_center = 7;

% Compute node degrees (number of significant networks per phenotype)
left_counts = groupsummary(red_table, 'Abbreviation');
right_counts = groupsummary(syn_table, 'Abbreviation');

left_deg_map = containers.Map(left_counts.Abbreviation, left_counts.GroupCount);
right_deg_map = containers.Map(right_counts.Abbreviation, right_counts.GroupCount);

figure('Color', 'w', 'Position', [50, 50, 1300, 850]);
hold on; axis off;

x_left = 0.1;   x_center = 0.5;   x_right = 0.9;
y_center_nodes = linspace(0.15, 0.85, num_center)';
y_left_nodes = linspace(0.05, 0.95, num_left)';
y_right_nodes = linspace(0.05, 0.95, num_right)';

left_map = containers.Map(left_phenos, y_left_nodes);
right_map = containers.Map(right_phenos, y_right_nodes);

networks = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN'};
Yeo7Color = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78] ./ 255;
t_max = max(abs(res_clean_all.t(res_clean_all.is_significant == 1)));

draw_curve = @(x1, y1, x2, y2, col, w) plot(...
    kron((1-linspace(0,1,50)).^2, x1) + kron(2*(1-linspace(0,1,50)).*linspace(0,1,50), (x1+x2)/2) + kron(linspace(0,1,50).^2, x2), ...
    kron((1-linspace(0,1,50)).^2, y1) + kron(2*(1-linspace(0,1,50)).*linspace(0,1,50), y1) + kron(linspace(0,1,50).^2, y2), ...
    '-', 'Color', [col, 0.45], 'LineWidth', w);

for i = 1:height(red_table)
    p_name = red_table.Abbreviation{i};
    y_p = left_map(p_name);
    for n = 1:num_center
        if contains(red_table.yname{i}, ['NC', num2str(n)])
            line_w = 0.5 + 3.5 * (abs(red_table.t(i)) / t_max);
            draw_curve(x_left, y_p, x_center, y_center_nodes(n), Yeo7Color(n,:), line_w);
        end
    end
end

for i = 1:height(syn_table)
    p_name = syn_table.Abbreviation{i};
    y_p = right_map(p_name);
    for n = 1:num_center
        if contains(syn_table.yname{i}, ['NC', num2str(n)])
            line_w = 0.5 + 3.5 * (abs(syn_table.t(i)) / t_max);
            draw_curve(x_center, y_center_nodes(n), x_right, y_p, Yeo7Color(n,:), line_w);
        end
    end
end

% Render Left Side Nodes (Size ranges dynamically from 15 to 120 based on degree)
for i = 1:num_left
    deg = left_deg_map(left_phenos{i});
    sz = 15 + (deg - 1) * 15; 
    scatter(x_left, y_left_nodes(i), sz, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.5);
    text(x_left - 0.015, y_left_nodes(i), left_phenos{i}, 'FontSize', 8.5, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end

% Render Right Side Nodes (Size ranges dynamically from 15 to 120 based on degree)
for i = 1:num_right
    deg = right_deg_map(right_phenos{i});
    sz = 15 + (deg - 1) * 15;
    scatter(x_right, y_right_nodes(i), sz, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.5);
    text(x_right + 0.015, y_right_nodes(i), right_phenos{i}, 'FontSize', 8.5, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
end

for n = 1:num_center
    scatter(x_center, y_center_nodes(n), 450, 'MarkerFaceColor', Yeo7Color(n,:), 'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 1.2);
    text(x_center, y_center_nodes(n), networks{n}, 'FontSize', 9.5, 'FontWeight', 'bold', 'Color', 'w', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

text(x_left, 0.98, 'Red7nNC Associated Phenotypes', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
text(x_right, 0.98, 'Syn7nNC Associated Phenotypes', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
text(x_center, 0.98, 'Brain Networks', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

xlim([x_left - 0.2, x_right + 0.2]); ylim([0, 1]);




%% Biomarker Association

UKB_Biomarker_0_0 = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_Biomarker_0_0.csv');
% 7 networks
UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_Biomarker_0_0),brain_data_7n);
UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; 

formula = ['AgeAttend_2_0+AgeAttend_2_0^2+AgeGap+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = UKB_Biomarker_0_0.Properties.VariableNames(2:end)'; 
results_biomarker_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

results_biomarker_7n.xname = strrep(results_biomarker_7n.xname,'_0','');
results_biomarker_7n.xname = strrep(results_biomarker_7n.xname,'_count','#')
results_biomarker_7n.xname = strrep(results_biomarker_7n.xname,'_pct','%')

results_biomarker_7n(isnan(results_biomarker_7n.t),:) = [];

p_thres = max(results_biomarker_7n.p(mafdr(results_biomarker_7n.p,'BHFDR',true)<0.05));

Yeo7Color = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78] ./ 255;

resRed = results_biomarker_7n(contains(results_biomarker_7n.yname,'Red7n'),:);
iv_Code = repmat([1:7],1,size(resRed,1)/7)';
plotVolcano(resRed.t, resRed.p, 'T', p_thres , iv_Code, Yeo7Color,resRed.xname,20,[700,550]);

resSyn = results_biomarker_7n(contains(results_biomarker_7n.yname,'Syn7n'),:);
iv_Code = repmat([1:7],1,size(resSyn,1)/7)';
plotVolcano(resSyn.t, resSyn.p, 'T', p_thres , iv_Code, Yeo7Color,resSyn.xname,20,[700,550]);






% 7x7 network-modules
UKB_BrainXY = innerjoin(UKB_Biomarker,brain_data_7n_matrix);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; 

formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = UKB_Biomarker.Properties.VariableNames(95:end)'; 
results_biomarker_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);

% writetable(results_biomarker_7n,'results_biomarker_7n.txt','Delimiter' ,',')
% writetable(results_biomarker_7n_matrix,'results_biomarker_7n_matrix.txt','Delimiter' ,',')

varNames = { ...
    'results_demo_360p', 'results_demo_7n', 'results_demo_7n_matrix', ...
    'results_lifestyle_360p', 'results_lifestyle_7n', 'results_lifestyle_7n_matrix', ...
    'results_pheno_360p', 'results_pheno_7n', 'results_pheno_7n_matrix', ...
    'results_pheno_adult_360p', 'results_pheno_adult_7n', 'results_pheno_adult_7n_matrix', ...
    'results_pheno_older_360p', 'results_pheno_older_7n', 'results_pheno_older_7n_matrix' ...
};

writePath = '/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/';
for i = 1:length(varNames)
    varName = varNames{i};
    if evalin('base', sprintf('exist(''%s'', ''var'')', varName))
        T = eval(varName);
        fileName = [writePath,sprintf('%s.txt', varName)];
        writetable(T, fileName, 'Delimiter', ' ');
        fprintf('Saved: %s\n', fileName);
    else
        warning('Variable %s does not exist in workspace, skipped.', varName);
    end
end

figure;
x = 'x709';lim = 8;
results_mat = results_lifestyle_7n_matrix;
tMat = icatb_vec2mat(results_mat.t(strcmp(results_mat.xname,x)&(contains(results_mat.yname,'RedYeo7n_Btw'))));
tMat(find(eye(7))) = results_mat.t(strcmp(results_mat.xname,x)&(contains(results_mat.yname,'RedYeo7n_In')));
pMat = icatb_vec2mat(results_mat.p(strcmp(results_mat.xname,x)&(contains(results_mat.yname,'RedYeo7n_Btw'))));
pMat(find(eye(7))) = results_mat.p(strcmp(results_mat.xname,x)&(contains(results_mat.yname,'RedYeo7n_In')));
imagesc(tMat.*(pMat<0.05/28),[-1*lim,lim]);set(gca,'DataAspectRatio',[1 1 1]);colorbar;colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)))
set(gca, 'XTick', 1:7, 'XTickLabel', Yeo7Names);set(gca, 'YTick', 1:7, 'YTickLabel', Yeo7Names);xtickangle(90)


tMat = icatb_vec2mat(results_mat.t(strcmp(results_mat.xname,x)&(contains(results_mat.yname,'SynYeo7n_Btw'))));
tMat(find(eye(7))) = results_mat.t(strcmp(results_mat.xname,x)&(contains(results_mat.yname,'SynYeo7n_In')));
pMat = icatb_vec2mat(results_mat.p(strcmp(results_mat.xname,x)&(contains(results_mat.yname,'SynYeo7n_Btw'))));
pMat(find(eye(7))) = results_mat.p(strcmp(results_mat.xname,x)&(contains(results_mat.yname,'SynYeo7n_In')));
imagesc(tMat.*(pMat<0.05/28),[-1*lim,lim]);set(gca,'DataAspectRatio',[1 1 1]);colorbar;colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)))
set(gca, 'XTick', 1:7, 'XTickLabel', Yeo7Names);set(gca, 'YTick', 1:7, 'YTickLabel', Yeo7Names);xtickangle(90)


%% Mediation Anlysis
UKB_BrainXY = innerjoin(innerjoin(innerjoin(innerjoin(UKB_Basic, UKB_BrainMRICov),UKB_Lifestyle_0_0),UKB_Outcomes_2_0),brain_data_7n);
YList = {'x20016_2_0','x4526_2_0','x2178_2_0'};
XList = unique(results_lifestyle_7n.xname(mafdr(results_lifestyle_7n.p,'BHFDR',true)<0.001));
MList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
Centre = [AsDummy(UKB_BrainXY.Centre_0_0)];
Centre(:,[sum(Centre)<1000]) = [];
COVList = {'AgeAttend_2_0','Sex_0_0','HeadMotion_2_0',...
    'BMI_2_0','Race_1','Race_2','Race_3','Race_4','Vol_WB_TIV_2_0'};

% === Initialize Variables ===
num_X = length(XList);
num_Y = length(YList);
num_M = length(MList);
total_runs = num_X * num_Y * num_M;

% Create a cell array to store pre-constructed individual tables safely within parfor
results_cell = cell(total_runs, 1);

% Extract the covariate matrix beforehand
COV0 = [UKB_BrainXY{:, COVList}, Centre];

% Define the standardized table variable names and types for initialization
varNames = {'X_Name', 'M_Name', 'Y_Name', ...
            'Beta_a', 't_a', 'p_a', ...
            'Beta_b', 't_b', 'p_b', ...
            'Beta_c_prime', 't_c_prime', 'p_c_prime', ...
            'Beta_ab', 't_ab', 'p_ab', ...
            'Mediation_Proportion', 'Status'};

fprintf('Starting %d Mediation Analyses via parfor...\n', total_runs);
% === 1. Initialize Variables and Numeric Matrices ===
num_X = length(XList);
num_Y = length(YList);
num_M = length(MList);
total_runs = num_X * num_Y * num_M;

% Use pure numeric matrices for stable parfor allocation (Avoids memory crash)
% 13 columns for all stats: beta_a, t_a, p_a, beta_b, t_b, p_b, c_prime, t_c_prime, p_c_prime, ab, t_ab, p_ab, prop
stats_matrix = zeros(total_runs, 13); 
% Store corresponding loop indices to reconstruct names later
idx_mapping = zeros(total_runs, 3); 

% Extract the covariate matrix beforehand to minimize worker communication overhead
COV0 = [UKB_BrainXY{:, COVList}, Centre];

fprintf('Starting %d Mediation Analyses safely (Memory-Optimized Pure Matrix)...\n', total_runs);
if ~isempty(gcp('nocreate'))
    pctRunOnAll maxNumCompThreads(1);
end

% === 1. Initialize Variables and Numeric Matrices ===
num_X = length(XList);
num_Y = length(YList);
num_M = length(MList);
total_runs = num_X * num_Y * num_M;

% 17¸öÍ³¼ÆÖ¸±ê£º[beta, t, p] x 5¸öpath (a,b,cp,c,ab) + 1¸öProportion + 1¸öType±êÇ©
stats_matrix = zeros(total_runs, 17); 
idx_mapping = zeros(total_runs, 3); 

% ÌáÇ°ÌáÈ¡Ð­±äÁ¿
COV0 = [UKB_BrainXY{:, COVList}, Centre];

fprintf('Starting %d Safe Mediation Analyses via ToolBox with parfor...\n', total_runs);

parfor idx = 1:total_runs
    % Deconstruct the flattened 1D index into m, i, j coordinates
    [m, i, j] = ind2sub([num_M, num_X, num_Y], idx);
    idx_mapping(idx, :) = [i, m, j]; % Track [X, M, Y] index
    
    % === Data extraction ===
    Xvar = UKB_BrainXY{:, XList{i}};
    Yvar = UKB_BrainXY{:, YList{j}};
    Mvar = UKB_BrainXY{:, MList{m}};
    
    % Remove rows with NaN
    [X2, Y2, M2, COV2] = nanrm(Xvar, Yvar, Mvar, COV0);
    
    % Remove near-constant covariates
    if ~isempty(COV2)
        covVar = var(COV2, 0, 1);
        constCols = covVar < 1e-8;
        if any(constCols)
            COV2(:, constCols) = [];
        end
    end
    
    % Z-score normalization
    X2 = zscore(X2); Y2 = zscore(Y2); M2 = zscore(M2);
    if ~isempty(COV2)
        COV2 = zscore(COV2);
    end
    
    % Run mediation safely inside try-catch
    try
        [paths, stats] = mediation(X2, Y2, M2, ...
            'covs', COV2, ...
            'verbose', 'no');  
        
        % Extract metrics from toolbox result structure
        beta_a  = stats.paths(1);  t_a  = stats.t(1);  p_a  = stats.p(1);
        beta_b  = stats.paths(2);  t_b  = stats.t(2);  p_b  = stats.p(2);
        beta_cp = stats.paths(3);  t_cp = stats.t(3);  p_cp = stats.p(3);
        beta_c  = stats.paths(4);  t_c  = stats.t(4);  p_c  = stats.p(4);
        beta_ab = stats.paths(5);  t_ab = stats.t(5);  p_ab = stats.p(5);
        
        % Calculate Mediation Proportion safely
        total_effect = beta_ab + beta_cp;
        prop_mediated = beta_ab / (total_effect + 1e-12);
        
        % Classify Mediation Type (0: Non-sig, 1: Complementary, 2: Suppression, 3: Full)
        med_type = 0;
        if p_ab < 0.05
            if p_cp >= 0.05
                med_type = 3; % Full Mediation
            elseif sign(beta_ab) == sign(beta_cp)
                med_type = 1; % Complementary Mediation
            else
                med_type = 2; % Competitive/Suppression
            end
        end
        
        % Pack everything into the numeric matrix row
        stats_matrix(idx, :) = [beta_a, t_a, p_a, ...
                                 beta_b, t_b, p_b, ...
                                 beta_cp, t_cp, p_cp, ...
                                 beta_c, t_c, p_c, ...
                                 beta_ab, t_ab, p_ab, ...
                                 prop_mediated, med_type];
    catch
        % On failure, fill row with NaNs
        stats_matrix(idx, :) = NaN(1, 17);
    end
end

fprintf('\nParallel computing finished! Assembling final elegant wide table...\n');

% Filter out combinations that failed completely (where all stats are NaN)
valid_idx = ~isnan(stats_matrix(:, 1));

X_Names_Final = XList(idx_mapping(valid_idx, 1))';
M_Names_Final = MList(idx_mapping(valid_idx, 2));
Y_Names_Final = YList(idx_mapping(valid_idx, 3))';
valid_stats = stats_matrix(valid_idx, :);

% Dynamic string mapping for mediation types
type_labels = {'Not Significant', 'Complementary', 'Suppression', 'Full'};
raw_types = valid_stats(:, 17);
raw_types(isnan(raw_types) | raw_types < 0 | raw_types > 3) = 0; % Fail-safe for indices
Mediation_Type = type_labels(raw_types + 1)';

% Construct Table
base_info_table = table(X_Names_Final, M_Names_Final, Y_Names_Final, Mediation_Type, ...
    'VariableNames', {'X_Name', 'M_Name', 'Y_Name', 'Mediation_Type'});

core_stats_table = array2table(valid_stats(:, 1:16), 'VariableNames', {...
    'Beta_path_a', 'T_path_a', 'P_path_a', ...
    'Beta_path_b', 'T_path_b', 'P_path_b', ...
    'Beta_path_cp', 'T_path_cp', 'P_path_cp', ...
    'Beta_path_c', 'T_path_c', 'P_path_c', ...
    'Beta_Mediation_ab', 'T_Mediation_ab', 'P_Mediation_ab', ...
    'Mediation_Proportion'});

% Stitch everything together
mediation_results_table = [base_info_table, core_stats_table];

% Print preview of the final dataset
disp(head(mediation_results_table));


fprintf('\nParallel computing finished! Assembling final table now...\n');

% Map the numeric indices back to their actual string names
X_Names_Final = XList(idx_mapping(:, 1))';
M_Names_Final = MList(idx_mapping(:, 2));
Y_Names_Final = YList(idx_mapping(:, 3))';

% Convert the core statistics to a table
stats_table = array2table(stats_matrix, 'VariableNames', {...
    'Beta_a', 't_a', 'p_a', ...
    'Beta_b', 't_b', 'p_b', ...
    'Beta_c_prime', 't_c_prime', 'p_c_prime', ...
    'Beta_ab', 't_ab', 'p_ab', ...
    'Mediation_Proportion'});

% Stitch names and statistics together
mediation_results_table = [table(X_Names_Final, M_Names_Final, Y_Names_Final, ...
    'VariableNames', {'X_Name', 'M_Name', 'Y_Name'}), stats_table];

% Display the first few rows of the final beautiful table
disp(head(mediation_results_table));

YList = {'x20016_2_0','x4526_2_0'};
XList = unique(results_lifestyle_7n.xname(mafdr(results_lifestyle_7n.p,'BHFDR',true)<0.001));
MList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
CList = {'Sex_0_0','HeadMotion_2_0','Euler_Num_WB_2_0','BMI_2_0','Vol_WB_TIV_2_0','Race_1','Race_2','Race_3','Race_4'};
Centre_mat = [AsDummy(UKB_BrainXY.Centre_0_0)];
Centre_mat(:, [sum(Centre_mat) < 1000]) = [];


resultMediationTable = batchMediationAnalysis(UKB_BrainXY, XList, YList, MList, CList,'Centre',Centre_mat);
resMediationEff = resultMediationTable(mafdr(resultMediationTable.P_Mediation_ab,'BHFDR',true)<0.05,:);



res_clean = resultMediationTable;
refInfo = UKB_LifestyleList;max_abs_t = 5.2;

fdr_p = mafdr(res_clean.P_Mediation_ab, 'BHFDR', true);
res_clean.is_significant = fdr_p < 0.05;
res_clean = res_clean(contains(res_clean.Y_Name,'x4526_2_0'),:)

t_table = unstack(res_clean(:, {'X_Name', 'M_Name', 'T_Mediation_ab'}), 'T_Mediation_ab', 'M_Name');
sig_table = unstack(res_clean(:, {'X_Name', 'M_Name', 'is_significant'}), 'is_significant', 'M_Name');

% Extract unique environment factors (X-axis raw keys: e.g., 'x738_0_0')
x_labels_raw = t_table.X_Name; 
y_labels = t_table.Properties.VariableNames(2:end); 

% Convert tables to numeric matrices
t_matrix = table2array(t_table(:, 2:end));
sig_matrix = table2array(sig_table(:, 2:end));
sig_matrix(isnan(sig_matrix)) = 0; % Clear NaN padding from unstack

x_labels_new = cell(size(x_labels_raw));
for i = 1:length(x_labels_raw)
    % Use regex to extract the main FieldID number between 'x' and the first '_'
    % For example, 'x738_0_0' will yield tokens{1} = '738'
    tokens = regexp(x_labels_raw{i}, '^x(\d+)_', 'tokens');
    
    if ~isempty(tokens)
        field_id_num = str2double(tokens{1}{1});
        
        % Look up the matching row in refInfo (UKB_OutcomesList)
        match_idx = find(refInfo.FieldID == field_id_num, 1);
        
        if ~isempty(match_idx)
            % If matched, extract Description (and optionally append Instance info if desired)
            x_labels_new{i} = refInfo.Abbreviation{match_idx};
        else
            % Fallback to raw label if no match found in refInfo
            x_labels_new{i} = x_labels_raw{i};
        end
    else
        x_labels_new{i} = x_labels_raw{i};
    end
end
y_labels_new = {'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN', ...
                'VN', 'SMN', 'DAN', 'VAN', 'LN', 'FPN', 'DMN'};
t_matrix_final = t_matrix';
sig_matrix_final = sig_matrix';

figure('Color', 'w', 'Position', [100, 100, 1050, 650]); % Wider layout for Descriptions
h_im = imagesc(t_matrix_final);
hold on;
% Apply user-defined Blue-Red colormap
colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)));
axis image; 
[num_rows, num_cols] = size(t_matrix_final);

% Overlay custom black grid lines
for r = 0:num_rows
    plot([0.5, num_cols+0.5], [r+0.5, r+0.5], 'k-', 'LineWidth', 0.8);
end
for c = 0:num_cols
    plot([c+0.5, c+0.5], [0.5, num_rows+0.5], 'k-', 'LineWidth', 0.8);
end

% Overlay significance markers (Asterisks)
for r = 1:num_rows
    for c = 1:num_cols
        if sig_matrix_final(r, c) == 1
            text(c, r, '*', 'Color', 'k', 'FontSize', 16, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
        end
    end
end

% Set the beautifully converted Descriptions on X-axis and Networks on Y-axis
set(gca, 'XTick', 1:num_cols, 'XTickLabel', x_labels_new, 'TickLabelInterpreter', 'none');
set(gca, 'YTick', 1:num_rows, 'YTickLabel', y_labels_new, 'TickLabelInterpreter', 'none');
% Rotate X-axis text by 45 or 90 degrees (45 recommended since descriptions are longer text)
xtickangle(90);
box on;
set(gca, 'Layer', 'top', 'LineWidth', 1);
% Configure the symmetric colorbar
c_bar = colorbar;
c_bar.Label.String = 't-statistic';
c_bar.Label.FontSize = 11;
caxis([-max_abs_t, max_abs_t]); 

title('Association Map (FDR Correction)', ...
      'FontSize', 13, 'FontWeight', 'bold');


% writetable(T,'/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/mediation_results_7nx.txt')

%% Interaction Effect
UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_BrainMRICov,UKB_Lifestyle_0_0),UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
formula = 'AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)';

X2List = unique(results_lifestyle_7n.xname(mafdr(results_lifestyle_7n.p,'BHFDR',true)<0.001));
% X2List = [unique(res_clean_all.xname(mafdr(res_clean_all.p, 'BHFDR', true)<0.05))];
X1List = {'AgeAttend_2_0'};
YList = brain_data_7n.Properties.VariableNames(2:end)';
results_interact_7n = LMM_withFormula_Interact(UKB_BrainXY,formula,X1List,X2List,YList);

CList = {'Sex_0_0','BMI_2_0','HeadMotion_2_0','Race_1','Race_2','Race_3','Race_4','Vol_WB_TIV_2_0','Centre_0_0'};
plot_interEff(UKB_BrainXY, 'Red7n_NC2', 'AgeAttend_2_0', 'x904_2_0', CList, Yeo7Color(2,:))


allResults = table();
for i = 1:14
    residualModel = fitlme(UKB_BrainXY, ...
        [YList{i}, ' ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + Vol_WB_TIV_2_0 + (1|Centre_0_0)']);
    UKB_BrainXY.Residual = residuals(residualModel);
    interactionModel = fitlm(UKB_BrainXY, ...
        'x20016_2_0 ~ Residual + AgeAttend_2_0 + AgeAttend_2_0 * Residual');
    coefTable = interactionModel.Coefficients;
    coefTable.YVariable = repmat(YList{i}, height(coefTable), 1);
    estimate = coefTable.Estimate;
    se = coefTable.SE;
    tStat = coefTable.tStat;
    pValue = coefTable.pValue;
    rowNames = {'(Intercept)';'AgeAttend_2_0';'Residual';'AgeAttend_2_0:Residual'}; 
    YVariable = repmat(YList(i), length(rowNames), 1);
    tempTable = table(YVariable, rowNames,estimate, se, tStat, pValue,  'VariableNames', {'YVariable', 'RowName','Estimate', 'SE', 'tStat', 'pValue'});
    allResults = [allResults; tempTable];
end
disp(allResults);


X = UKB_BrainXY.Residual; 
Z = UKB_BrainXY.AgeAttend_2_0;     
Y = UKB_BrainXY.x20016; 

Z_mean = mean(Z);
Z_std = std(Z);

Z_low = Z_mean - Z_std;
Z_mid = Z_mean;
Z_high = Z_mean + Z_std;

X_range = linspace(min(X), max(X), 100); 
Y_low = interactionModel.Coefficients.Estimate(1) + ... % Intercept
    interactionModel.Coefficients.Estimate(2) * Z_low + ...
    interactionModel.Coefficients.Estimate(3) * X_range + ...
    interactionModel.Coefficients.Estimate(4) * (Z_low .* X_range);
Y_mid = interactionModel.Coefficients.Estimate(1) + ... % Intercept
    interactionModel.Coefficients.Estimate(2) * Z_mid + ...
    interactionModel.Coefficients.Estimate(3) * X_range + ...
    interactionModel.Coefficients.Estimate(4) * (Z_mid .* X_range);
Y_high = interactionModel.Coefficients.Estimate(1) + ... % Intercept
    interactionModel.Coefficients.Estimate(2) * Z_high + ...
    interactionModel.Coefficients.Estimate(3) * X_range + ...
    interactionModel.Coefficients.Estimate(4) * (Z_high .* X_range);


figure;
plot(X_range, Y_low, 'b-', 'LineWidth', 1.5); hold on;
plot(X_range, Y_mid, 'g-', 'LineWidth', 1.5);
plot(X_range, Y_high, 'r-', 'LineWidth', 1.5);
xlabel('Residual');
ylabel('Fluid Intelligence Score');
legend(sprintf('Age = %.2f (Low)', Z_low), ...
sprintf('Age = %.2f (Mid)', Z_mid), ...
sprintf('Age = %.2f (High)', Z_high), ...
'Location', 'best');
title('Interaction Effect');
grid on;

figure;
for i = 1:7
    subplot(1, 7, i); 
    hold on;
    
    residualModel = fitlme(UKB_BrainXY, ...
        [YList{i}, ' ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + Vol_WB_TIV_2_0 + (1|Centre_0_0)']);
    UKB_BrainXY.Residual = residuals(residualModel);
    
    interactionModel = fitlm(UKB_BrainXY, ...
        'Residual ~ x20016 + AgeAttend_2_0 + AgeAttend_2_0 * x20016');
    
    Y = UKB_BrainXY.Residual; 
    X = UKB_BrainXY.AgeAttend_2_0; 
    Z = UKB_BrainXY.x20016; 
    
    Z_mean = mean(Z);
    Z_std = std(Z);
    
    Z_low = Z_mean - Z_std;
    Z_mid = Z_mean;
    Z_high = Z_mean + Z_std;
   
    X_range = linspace(min(X), max(X), 100); 
    
    Y_low = interactionModel.Coefficients.Estimate(1) + ...
        interactionModel.Coefficients.Estimate(2) * Z_low + ...
        interactionModel.Coefficients.Estimate(3) * X_range + ...
        interactionModel.Coefficients.Estimate(4) * (Z_low .* X_range);
    
    Y_mid = interactionModel.Coefficients.Estimate(1) + ...
        interactionModel.Coefficients.Estimate(2) * Z_mid + ...
        interactionModel.Coefficients.Estimate(3) * X_range + ...
        interactionModel.Coefficients.Estimate(4) * (Z_mid .* X_range);
    
    Y_high = interactionModel.Coefficients.Estimate(1) + ...
        interactionModel.Coefficients.Estimate(2) * Z_high + ...
        interactionModel.Coefficients.Estimate(3) * X_range + ...
        interactionModel.Coefficients.Estimate(4) * (Z_high .* X_range);
    
    SE_low = sqrt(interactionModel.Coefficients.SE(1)^2 + interactionModel.Coefficients.SE(2)^2 + interactionModel.Coefficients.SE(3)^2 + interactionModel.Coefficients.SE(4)^2); % ¼ò»¯°æ±¾
    SE_mid = SE_low;
    SE_high = SE_low; 
    
    CI_low_low = Y_low - 1.96 * SE_low;
    CI_low_high = Y_low + 1.96 * SE_low;
    
    CI_mid_low = Y_mid - 1.96 * SE_mid;
    CI_mid_high = Y_mid + 1.96 * SE_mid;
    
    CI_high_low = Y_high - 1.96 * SE_high;
    CI_high_high = Y_high + 1.96 * SE_high;
    
    col_low = Yeo7Colormap(i, :); 
    col_mid = Yeo7Colormap(i, :); 
    col_high = Yeo7Colormap(i, :); 
    
    fill([X_range, fliplr(X_range)], [CI_low_low, fliplr(CI_low_high)], [0.7, 0.7, 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    fill([X_range, fliplr(X_range)], [CI_mid_low, fliplr(CI_mid_high)], [0.7, 0.7, 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    fill([X_range, fliplr(X_range)], [CI_high_low, fliplr(CI_high_high)],[0.7, 0.7, 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
    plot(X_range, Y_low, '-.', 'Color', col_low, 'LineWidth', 2.5);
    plot(X_range, Y_mid, '--', 'Color', col_mid, 'LineWidth', 2.5);
    plot(X_range, Y_high, '-', 'Color', col_high, 'LineWidth', 2.5);
    
    xlabel('Residual');
    ylabel('Fluid Intelligence');
    title(Yeo7Names2{i});

    grid on;
end




%% GWAS

Xpheno = [array2table(brain_data_7n.eid),brain_data_7n];
Xpheno.Properties.VariableNames{2} = 'IID';Xpheno.Properties.VariableNames{1} = 'FID';
writetable(Xpheno,'/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/UKB_redsyn_pheno.txt','Delimiter' ,' ')

cd('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/')
phenofile = '/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/UKB_redsyn_pheno.txt';

covarfile = '/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/UKB_GeneCov_MRI.txt';
covarname = ['Sex_0_0,AgeAttend_2_0,BMI_2_0,HeadMotion_2_0,Vol_WB_TIV_2_0,Centre_2_0,PC1-PC20'];

n = Xpheno.Properties.VariableNames(3:end)';
parfor chr = 1:22
    bfile_qc = ['/public/home/ISTBI_data/UKB/gene/v3/QCed/ukb_imp_chr',num2str(chr),'_v3'];
    outfile = ['/public/home/zhangjie/ZJLab/UKBiobank_Project/project/',...
        'Detecting_Compensation_Effects_of_Aging_Using_PID/GWAS_result/gwas_chr',num2str(chr),'_7n'];
    plink_command = ...
        ['/public/home/zhangjie/DataAnalysis/wxr_toolbox/toobox_other/plink2_linux_x86_64_20191030/plink2 ',...
        '--bfile ',bfile_qc,' ',...
        '--linear  hide-covar --debug ',...
        '--pheno ',phenofile,' ',...
        '--pheno-name Red7n_NC1,Red7n_NC2,Red7n_NC3,Red7n_NC4,Red7n_NC5,Red7n_NC6,Red7n_NC7,Syn7n_NC1,Syn7n_NC2,Syn7n_NC3,Syn7n_NC4,Syn7n_NC5,Syn7n_NC6,Syn7n_NC7 '...
        '--covar ',covarfile,' ',...
        '--covar-name ',covarname,' ',...
        '--covar-variance-standardize ',...
        '--out ',outfile];
    unix(plink_command);
    disp(['GWAS of in chr',num2str(chr)])
end

PList = {'Red7n_NC1','Red7n_NC2','Red7n_NC3','Red7n_NC4','Red7n_NC5','Red7n_NC6','Red7n_NC7','Syn7n_NC1','Syn7n_NC2','Syn7n_NC3','Syn7n_NC4','Syn7n_NC5','Syn7n_NC6','Syn7n_NC7'}
cd('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/GWAS_result/')

for i = 1:14
    
    unix(['head -n 1 gwas_chr1_7n.',PList{i},'.glm.linear > gwas.',PList{i},'.glm.linear']);
    unix(['tail -n +2 gwas_chr1_7n.',PList{i},'.glm.linear >> gwas.',PList{i},'.glm.linear']);
    
    for chr = 2:22
        unix(['tail -n +2 gwas_chr',num2str(chr),'_7n.',PList{i},'.glm.linear >> gwas.',PList{i},'.glm.linear']);
    end
    
    unix(['sed -i ','''','1s/.*/CHR POS ID REF ALT PROVISIONAL_REF A1 OMITTED A1_FREQ TEST OBS_CT BETA SE T_STAT P ERRCODE/',...
        '''',' gwas.',PList{i},'.glm.linear']);
    
    unix(['gzip gwas.',PList{i},'.glm.linear'])
end

gwas_data = readtable('');
gwas_data = sortrows(gwas_data, 'P', 'ascend');

unique_phenos = unique(gwas_data.PHENO);
lead_snps = table(); 
window_bp = 500000;  
for i = 1:length(unique_phenos)
    sub_df = gwas_data(strcmp(gwas_data.PHENO, unique_phenos{i}), :);
    while ~isempty(sub_df)
        lead = sub_df(1, :);
        lead_snps = [lead_snps; lead];
        in_window = (sub_df.CHR == lead.CHR) & (abs(sub_df.pos - lead.pos) <= window_bp);
        sub_df(in_window, :) = []; 
    end
end
disp(lead_snps(:, {'SNP', 'CHR', 'pos', 'PHENO', 'Beta', 'P'}));

%% APOE4

APOE4_genetype = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/gene/APOE4_genetype.txt')
APOE4_genetype = removevars(APOE4_genetype, 'FID');
APOE4_genetype = removevars(APOE4_genetype, {'PAT','MAT','SEX','PHENOTYPE'});

% 7 networks
UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_Lifestyle_0_0,APOE4_genetype),UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

XList = {'allele1','APOE4_dosage'};

formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];

results_APOE4_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_APOE4_adult_7n = LMM_with_EffectSizes_Full(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0<65,:),formula,XList,YList);
results_APOE4_older_7n = LMM_with_EffectSizes_Full(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0>=65,:),formula,XList,YList);


%% Clinical Diagnosis

UKB_ICD10_TimeGap_2_0 = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_ICD10_TimeGap_2_0.csv', 'Delimiter', ',');
UKB_ICD10_TimeGap_2_0.Q99 = str2double(UKB_ICD10_TimeGap_2_0.Q99);
ICD10_Coding = UKB_ICD10_TimeGap_2_0.Properties.VariableNames(2:end);

UKB_ICD10_Past_2_0 = UKB_ICD10_TimeGap_2_0;
UKB_ICD10_Future_2_0 = UKB_ICD10_TimeGap_2_0;

XPast = UKB_ICD10_TimeGap_2_0{:,2:end}<0;
XFuture = UKB_ICD10_TimeGap_2_0{:,2:end}>0;
XHealthy = UKB_ICD10_TimeGap_2_0{:,2:end}>0;
XPastCont = zeros(size(UKB_ICD10_TimeGap_2_0{:,2:end}));
XPastCont(XFuture) = NaN;
XPastCont(XPast) = 1;

XFutureCont = UKB_ICD10_TimeGap_2_0{:,2:end};
XFutureCont(XPast) = NaN;

UKB_ICD10_Past_2_0{:,2:end} = XPastCont;
UKB_ICD10_Future_2_0{:,2:end} = XFutureCont;

UKB_BrainMRICov = readtable('/home1/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BrainMRICov.txt');

% 360 regions * age
UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_ICD10_Past_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; 

formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];





XList = ICD10_Coding(nansum(UKB_BrainXY{:,95:1225},1)>1000)';
results_ICD10Past_7n = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_ICD10Past_7n(isnan(results_ICD10Past_7n.t),:) = [];


UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_ICD10_TimeGap_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = ICD10_Coding(nansum(~isnan(UKB_BrainXY{:,95:1225}),1)>1000)';
XList = brain_data_7n.Properties.VariableNames(2:end)'; 
CList = {'AgeAttend_2_0','Sex_0_0','HeadMotion_2_0','BMI_2_0','Race_1','Race_2','Race_3','Race_4','Vol_WB_TIV_2_0','Centre_0_0'};
results_ICD10Cox_7n = CoxRegression(UKB_BrainXY, XList,  YList, CList);



p_thres = max(results_ICD10Past_7n.p(mafdr(results_ICD10Past_7n.p,'BHFDR',true)<0.05));


Yeo7Color = [120 18 134; 70 130 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78] ./ 255;

resRed = results_ICD10Past_7n(contains(results_ICD10Past_7n.yname,'Red7n'),:);
iv_Code = repmat([1:7],1,size(resRed,1)/7)';
plotVolcano(resRed.t, resRed.p, 'T', p_thres, iv_Code, Yeo7Color,resRed.xname,20);

resSyn = results_ICD10Past_7n(contains(results_ICD10Past_7n.yname,'Syn7n'),:);
iv_Code = repmat([1:7],1,size(resSyn,1)/7)';
plotVolcano(resSyn.t, resSyn.p, 'T', p_thres, iv_Code, Yeo7Color,resSyn.xname,20);



p_thres = max(results_ICD10Cox_7n.P_Value(mafdr(results_ICD10Cox_7n.P_Value,'BHFDR',true)<0.05));



resRed = results_ICD10Cox_7n(contains(results_ICD10Cox_7n.Independent_Var,'Red7n'),:);
iv_Code = repmat([1:7],1,size(resRed,1)/7)';
plotVolcano(resRed.HR, resRed.P_Value, 'HR', p_thres,  iv_Code, Yeo7Color,resRed.Survival_Var,20);

resSyn = results_ICD10Cox_7n(contains(results_ICD10Cox_7n.Independent_Var,'Syn7n'),:);
iv_Code = repmat([1:7],1,size(resSyn,1)/7)';
plotVolcano(resSyn.HR, resSyn.P_Value, 'HR', p_thres, iv_Code, Yeo7Color,resSyn.Survival_Var,20);


XList = unique(results_ICD10Past_7n.xname(mafdr(results_ICD10Past_7n.p,'BHFDR',true)<0.0001));

UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_ICD10_Past_2_0),brain_data_7n_matrix);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; 
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
results_ICD10Past_7n_matrix = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_ICD10Past_7n_matrix(isnan(results_ICD10Past_7n_matrix.t),:) = [];


UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_ICD10_Past_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; 
formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
results_ICD10Past_360p = LMM_with_EffectSizes_Full(UKB_BrainXY,formula,XList,YList);
results_ICD10Past_360p(isnan(results_ICD10Past_360p.t),:) = [];


V = 'F32';
Xt = results_ICD10Past_7n_matrix.t(contains(results_ICD10Past_7n_matrix.xname,V));
Xp = results_ICD10Past_7n_matrix.p(contains(results_ICD10Past_7n_matrix.xname,V));
X0 = Xt.*(Xp<(0.05/28));
X01 = X0(1:28);
X02 = X0(30:57)
X1 = icatb_vec2mat(X01(8:28));X1(find(eye(7))) = X01(1:7);
plot_matrix_Yeo7n(X1,10)
X2 = icatb_vec2mat(X02(8:28));X2(find(eye(7))) = X02(1:7);
plot_matrix_Yeo7n(X2,10)

for i = 1:6
    V = XList{i};
    Xt = results_ICD10Past_7n_matrix.t(contains(results_ICD10Past_7n_matrix.xname,V));
    Xp = results_ICD10Past_7n_matrix.p(contains(results_ICD10Past_7n_matrix.xname,V));
    X0 = Xt.*(Xp<(0.05/28));
    X01 = X0(1:28);
    X02 = X0(30:57);
    X1 = icatb_vec2mat(X01(8:28));X1(find(eye(7))) = X01(1:7);
    plot_matrix_Yeo7n(X1,10);title([V,'-syn'])
    X2 = icatb_vec2mat(X02(8:28));X2(find(eye(7))) = X02(1:7);
    plot_matrix_Yeo7n(X2,10);title([V,'-red'])
    
end

i = 4;
V = XList{i};
Xt = results_ICD10Past_360p.t(contains(results_ICD10Past_360p.xname,V));
Xp = results_ICD10Past_360p.p(contains(results_ICD10Past_360p.xname,V));
X0 = Xt.*(Xp<(0.05/360));
X01 = X0(1:360);
X02 = X0(361:end);
[a, cb] = plot_cortical(parcel_to_surface(X01([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-8.5,8.5],'layout','grid');

[a, cb] = plot_cortical(parcel_to_surface(X02([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-8.5,8.5],'layout','grid');




%% Interaction Effect of Lifestyle

UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];

X1List = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
% X1List = UKB_Lifestyle_0_0.Properties.VariableNames(95:end)';
YList = brain_data_7n.Properties.VariableNames(2:end)';
X2List_raw = {'AgeAttend_2_0'};
X2List = repmat(X2List_raw, length(X1List), 1);

results_interact_7n = LMM_withFormula_Interact(UKB_BrainXY,formula,X1List,X2List,YList);



UKB_BrainXY = innerjoin(innerjoin(UKB_BrainMRICov,UKB_ICD10_Past_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
X1List = ICD10_Coding(nansum(UKB_BrainXY{:,95:1225},1)>1000)';
YList = brain_data_7n.Properties.VariableNames(2:end)';
X2List_raw = {'AgeAttend_2_0'};
X2List = repmat(X2List_raw, length(X1List), 1);

results_interact_7n = LMM_withFormula_Interact(UKB_BrainXY,formula,X1List,X2List,YList);


%% Structural Mediation

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

UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_BrainMRICov,UKB_Freesurfer_APARC(:,[1,740,741,744])),UKB_Glasser360_thickness_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
UKB_BrainXY(UKB_BrainXY.Euler_Num_WB_2_0<-217,:) = [];
XList = repmat({'AgeAttend_2_0'},360,1);
YList = brain_data_360p.Properties.VariableNames(2:end)'; 
MList = UKB_Glasser360_thickness_2_0.Properties.VariableNames(2:end)' ;
CList = {'Sex_0_0','HeadMotion_2_0','Euler_Num_WB_2_0','BMI_2_0','Vol_WB_TIV_2_0','Race_1','Race_2','Race_3','Race_4'};

Centre_mat = [AsDummy(UKB_BrainXY.Centre_0_0)];
Centre_mat(:, [sum(Centre_mat) < 1000]) = [];

results_Mediation_360p_Redundancy_thk = batchMediationAnalysis_1to1(UKB_BrainXY, XList, YList(1:360), MList([181:360,1:180]), CList,'Centre',Centre_mat);
results_Mediation_360p_Synergy_thk = batchMediationAnalysis_1to1(UKB_BrainXY, XList, YList(360+[1:360]), MList([181:360,1:180]), CList,'Centre',Centre_mat);
X = results_Mediation_360p_Synergy_thk.T_Value(5:5:end);[a, cb] = plot_cortical(parcel_to_surface(X([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-5.5,5.5]);

cov_matrix = [UKB_BrainXY{:, [10,20,24,67:70,77,83,95:97]},Centre_mat];
X_matrix   = UKB_BrainXY{:, 98:457};
X_aligned  = X_matrix(:, [181:360, 1:180]); 

M_red_mat  = UKB_BrainXY{:, 458:817}; 
M_syn_mat  = UKB_BrainXY{:, 818:1177}; 

r_red_thk = zeros(360, 1);  p_red_thk = zeros(360, 1);
r_syn_thk = zeros(360, 1);  p_syn_thk = zeros(360, 1);

fprintf('Starting pinpoint pair-wise partial correlation for 360 brain nodes...\n');

for i = 1:360
    x_i = X_aligned(:, i);
    m_red_i = M_red_mat(:, i);
    m_syn_i = M_syn_mat(:, i);
    
    [r_r, p_r] = partialcorr(m_red_i, x_i, cov_matrix, 'Rows', 'pairwise');
    r_red_thk(i) = r_r;
    p_red_thk(i) = p_r;
    
    [r_s, p_s] = partialcorr(m_syn_i, x_i, cov_matrix, 'Rows', 'pairwise');
    r_syn_thk(i) = r_s;
    p_syn_thk(i) = p_s;
end


UKB_BrainXY = innerjoin(innerjoin(innerjoin(UKB_BrainMRICov,UKB_Freesurfer_APARC(:,[1,740,741,744])),UKB_Glasser360_area_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
UKB_BrainXY(UKB_BrainXY.Euler_Num_WB_2_0<-217,:) = [];
XList = repmat({'AgeAttend_2_0'},360,1);
YList = brain_data_360p.Properties.VariableNames(2:end)'; 
MList = UKB_Glasser360_thk_2_0.Properties.VariableNames(2:end)' ;
CList = {'Sex_0_0','HeadMotion_2_0','Euler_Num_WB_2_0','BMI_2_0','Vol_WB_TIV_2_0','Race_1','Race_2','Race_3','Race_4'};

Centre_mat = [AsDummy(UKB_BrainXY.Centre_0_0)];
Centre_mat(:, [sum(Centre_mat) < 1000]) = [];

results_Mediation_360p_Redundancy_area = batchMediationAnalysis_1to1(UKB_BrainXY, XList, YList(1:360), MList([181:360,1:180]), CList,'Centre',Centre_mat);
results_Mediation_360p_Synergy_area = batchMediationAnalysis_1to1(UKB_BrainXY, XList, YList(360+[1:360]), MList([181:360,1:180]), CList,'Centre',Centre_mat);
X = results_Mediation_360p_Redundancy_area.T_Value(5:5:end);[a, cb] = plot_cortical(parcel_to_surface(X([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-5.5,5.5]);

cov_matrix = [UKB_BrainXY{:, [10,20,24,67:70,77,83,95:97]},Centre_mat];
X_matrix   = UKB_BrainXY{:, 98:457};
X_aligned  = X_matrix(:, [181:360, 1:180]); 

M_red_mat  = UKB_BrainXY{:, 458:817}; 
M_syn_mat  = UKB_BrainXY{:, 818:1177}; 

r_red_area = zeros(360, 1);  p_red_area = zeros(360, 1);
r_syn_area = zeros(360, 1);  p_syn_area = zeros(360, 1);

fprintf('Starting pinpoint pair-wise partial correlation for 360 brain nodes...\n');

for i = 1:360
    x_i = X_aligned(:, i);
    m_red_i = M_red_mat(:, i);
    m_syn_i = M_syn_mat(:, i);
    
    [r_r, p_r] = partialcorr(m_red_i, x_i, cov_matrix, 'Rows', 'pairwise');
    r_red_area(i) = r_r;
    p_red_area(i) = p_r;
    
    [r_s, p_s] = partialcorr(m_syn_i, x_i, cov_matrix, 'Rows', 'pairwise');
    r_syn_area(i) = r_s;
    p_syn_area(i) = p_s;
end

X1 = r_red_thk.*(p_red_thk<0.05/360);
X2 = r_syn_thk.*(p_syn_thk<0.05/360);
X3 = r_red_area.*(p_red_area<0.05/360);
X4 = r_syn_area.*(p_syn_area<0.05/360);
X = [r_red_thk,r_syn_thk,r_red_area,r_syn_area];
X = [p_red_thk,p_syn_thk,p_red_area,p_syn_area];

[a, cb] = plot_cortical(parcel_to_surface(X3([181:360,1:180]),'glasser_360_fsa5'),'color_range',[-.02,.02]);

% resultTable = batchMediationAnalysis(dataTable, XList, YList, MList, COVList, varargin);
% resultTable = batchMediationAnalysis_1to1(dataTable, XList, YList, MList, COVList, varargin)

%% Forest Plot

%% Repeatiablty



load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/rotate_parcellation-master/perm_centroid_info_HCPMMP1.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP')
X = results_age_360p.t([1:360]);
null_X = X(perm_id);
X7n = netMean1D(X', Yeo7MMP);
X7n_null = netMean1D(null_X', Yeo7MMP);
X7n > X7n_null;
X7n_perm_p = nanmean(X7n > X7n_null);




load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/UKB_PID_Data.mat', 'Red360p_NC', 'Syn360p_NC')
Xsbj = netMean1D(Red360p_NC,Yeo7MMP);
for i = 1:1000
    Xsbj_null(:,:,i) = netMean1D(Red360p_NC,Yeo7MMP(perm_id(:,i)));
end

for j = 1:size(Red360p_NC,1)
    perm_p_sbj(j,:) = nanmean(Xsbj(j,:) > squeeze(Xsbj_null(j,:,:))');
    
end
hist(perm_p_sbj)
colormap(Yeo7Color)

Xsbj = netMean1D(Syn360p_NC,Yeo7MMP);
for i = 1:1000
    Xsbj_null(:,:,i) = netMean1D(Syn360p_NC,Yeo7MMP(perm_id(:,i)));
end

for j = 1:size(Syn360p_NC,1)
    perm_p_sbj(j,:) = nanmean(Xsbj(j,:) > squeeze(Xsbj_null(j,:,:))');
    
end
hist(perm_p_sbj)


for i = 1:7
    [h, p, stat] = kstest2(redundancy(Yeo7MMP==i), redundancy(Yeo7MMP!=i));
end


load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/red_syn/ukb_3_0_glasser_red_syn.mat', 'redundancy', 'synergy')



redundancy_all = nanmean(redundancy,3);
synergy_all = nanmean(synergy,3);

SimilarityMatrix = corr(redundancy_all, 'Type', 'Pearson');
data_clusters = kmeans(SimilarityMatrix, 7);
[aligned_xlabel, dice_matrix, mapping_table] = AlignClusters(data_clusters,Yeo7MMP)

SimilarityMatrix = corr(synergy_all, 'Type', 'Pearson');
data_clusters = kmeans(SimilarityMatrix, 7);
[aligned_xlabel, dice_matrix, mapping_table] = AlignClusters(data_clusters,Yeo7MMP)


load('/public/home/zhangjie/ZJLab/UKBiobank_Project/code/Redundancy_Synergy_Information/UKB_PhiID_All/data/ukb_2_0_dk_syn_red.mat')
DK_info_ukb = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/code/Redundancy_Synergy_Information/UKB_PhiID_All/data/DK_label_ukb.csv');

DK_Yeo7 = nan(84,1);
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'VIS')) = 1;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'SOM')) = 2;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'DAN')) = 3;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'VAN')) = 4;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'LIM')) = 5;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'FPN')) = 6;
DK_Yeo7(strcmp(DK_info_ukb.yeo_7, 'DMN')) = 7;
DK_Yeo7(69:84) = 8;

[RedYeo7n_aparc] = netMeans2D(redundancy,DK_Yeo7(1:68));
[SynYeo7n_aparc] = netMeans2D(synergy,DK_Yeo7(1:68));
RedYeo7n_In = wxr_mat2dia3d(RedYeo7n_aparc);
RedYeo7n_Btw = wxr_mat2vec3d(RedYeo7n_aparc);
SynYeo7n_In = wxr_mat2dia3d(SynYeo7n_aparc);
SynYeo7n_Btw = wxr_mat2vec3d(SynYeo7n_aparc);
UKB_RedSyn_aparc_2_0 = [array2table(ukb_2_0_BOLD_eID),array2table(RedYeo7n_In),array2table(RedYeo7n_Btw),...
    array2table(SynYeo7n_In),array2table(SynYeo7n_Btw)];
UKB_RedSyn_aparc_2_0.Properties.VariableNames{1} = 'eid';

ukb_RedSynYeo7n_2_0 = removevars(ukb_RedSynYeo7n_2_0, {'SynYeo7n_Total','RedYeo7n_Total'});

[b2f,f2b] = idFinderNum(ukb_RedSynYeo7n_2_0.eid,UKB_RedSyn_aparc_2_0.eid);
[r_ind p] = corr(UKB_RedSyn_aparc_2_0{f2b,[30:57,2:29]}, ukb_RedSynYeo7n_2_0{b2f,[2:29,30:57]})


resVec = diag(r_ind);resVec = resVec([1:28]);

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


z = find((aparc_68_fsa5.* glasser_360_fsa5)~=0);

parfor i = 1:length(f2b)
    X_aparc = nanmean(redundancy(:,:,f2b(i)));
    X_glasser = brain_data_360p{b2f(i),2:end}(1:360);
    aparc_68_fsa5 = parcel_to_surface(X_aparc,'aparc_fsa5');
    glasser_360_fsa5 = parcel_to_surface(X_glasser,'glasser_360_fsa5');
    
    [r_red_fs5(i,1),p_red_fs5(i,1)] = corr(aparc_68_fsa5(z)',glasser_360_fsa5(z)','Rows','pairwise','Type','spearman')
    disp(['sub-',num2str(i),': Spatial Correlation'])
end

parfor i = 1:length(f2b)
    X_aparc = nanmean(synergy(:,:,f2b(i)));
    X_glasser = brain_data_360p{b2f(i),2:end}(361:720);
    aparc_68_fsa5 = parcel_to_surface(X_aparc,'aparc_fsa5');
    glasser_360_fsa5 = parcel_to_surface(X_glasser,'glasser_360_fsa5');
    
    [r_syn_fs5(i,1),p_syn_fs5(i,1)] = corr(aparc_68_fsa5(z)',glasser_360_fsa5(z)','Rows','pairwise','Type','spearman')
    disp(['sub-',num2str(i),': Spatial Correlation'])
end

scatter(X,Xgradients)




X_aparc = nanmean(nanmean(redundancy(:,:,f2b)),3);
X_glasser = nanmean(brain_data_360p{b2f,2:end}(:,1:360),1);

aparc_68_fsa5 = parcel_to_surface(X_aparc,'aparc_fsa5');
glasser_360_fsa5 = parcel_to_surface(X_glasser,'glasser_360_fsa5');
aparc_68_fsa5(aparc_68_fsa5==0)=NaN;
glasser_360_fsa5(glasser_360_fsa5==0)=NaN;
% redundancy_aparc_68_fsa5 = parcel_to_surface(1:68,'aparc_fsa5');
% glasser_360_fsa5 = parcel_to_surface(1:360,'glasser_360_fsa5');
r_perm_red = corr(aparc_68_fsa5',glasser_360_fsa5','Rows','pairwise','Type','spearman');
p_perm_red = perm_sphere_p(aparc_68_fsa5',glasser_360_fsa5',perm_id(:,1:1000),'spearman');

X_aparc = nanmean(nanmean(synergy(:,:,f2b)),3);
X_glasser = nanmean(brain_data_360p{b2f,2:end}(:,361:720),1);

aparc_68_fsa5 = parcel_to_surface(X_aparc,'aparc_fsa5');
glasser_360_fsa5 = parcel_to_surface(X_glasser,'glasser_360_fsa5');
aparc_68_fsa5(aparc_68_fsa5==0)=NaN;
glasser_360_fsa5(glasser_360_fsa5==0)=NaN;
% redundancy_aparc_68_fsa5 = parcel_to_surface(1:68,'aparc_fsa5');
% glasser_360_fsa5 = parcel_to_surface(1:360,'glasser_360_fsa5');
r_perm_syn = corr(aparc_68_fsa5',glasser_360_fsa5','Rows','pairwise','Type','spearman');
p_perm_syn = perm_sphere_p(aparc_68_fsa5',glasser_360_fsa5',perm_id(:,1:1000),'spearman');

scatter(aparc_68_fsa5',glasser_360_fsa5')

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
z = find((aparc_68_fsa5.* glasser_360_fsa5)~=0);
scatter(aparc_68_fsa5(z), glasser_360_fsa5(z), 36, Yeo7_fsa5(z), 'filled');
colormap(Yeo7Color);
caxis([0.5 7.5]); 
xlabel('DK-68 (fsa5)');
ylabel('HCP-MMP (fsa5)');
grid on;
hcb = colorbar;
set(hcb, 'Ticks', 1:7, 'TickLabels', {'VN','SMN','DAN','VAN','LN','FPN','DMN'});
hold on;
x_data = aparc_68_fsa5(z);
y_data = glasser_360_fsa5(z);
p = polyfit(x_data, y_data, 1);
x_range = [min(x_data), max(x_data)];
y_fit = polyval(p, x_range);
plot(x_range, y_fit, 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, 'LineStyle', '--');
hold off;

load('/public/home/zhangjie/ZJLab/ZJLab_Toolset/rotate_parcellation-master/perm_centroid_info_fs5.mat');

[r p] = corr(aparc_68_fsa5,glasser_360_fsa5,'Rows','pairwise','Type','spearman');
for i = 1:1000
    [r_null(i,1) p_null(i,1)] = corr(aparc_68_fsa5,glasser_360_fsa5(perm_id(:,i)),'Rows','pairwise','Type','spearman');
end
