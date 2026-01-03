%% Data loading
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/ukb_synred_workspace_20250105.mat', 'ukb_RedSynYeo7n_2_0')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'Red360p_NC', 'Syn360p_NC')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/UKB_glasser_RedSyn_Yeo7n.mat', 'ukb_RedYeo7n_2_0', 'ukb_SynYeo7n_2_0')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/ukb_synred_workspace_20250105.mat', 'Yeo7MMP')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/ukb_synred_workspace_20250105.mat', 'modal_pheno')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/ukb_synred_workspace_20250105.mat', 'modal_pheno_7n')
load('/public/home/zhangjie/DataAnalysis/wxr_toolbox/toolbox_matlab/cbrewer2-master/cbrewer2/colorbrewer.mat', 'colorbrewer')
UKB_Outcomes_2_0 = readtable('/public/home/zhangjie/UKB_Outcomes_2_0.csv')
UKB_Lifestyle_0_0 = readtable('/public/home/zhangjie/UKB_Lifestyle_MRI_0_0.csv');
Yeo7Names = {'Visual';'Somatomotor';'DorsalAttention';'VentralAttention';'Limbic';'Frontoparietal';'Default'};

UKB_Lifestyle_0_0.x1160_L = double(UKB_Lifestyle_0_0.x1160 >= 9);
UKB_Lifestyle_0_0.x1160_S = double(UKB_Lifestyle_0_0.x1160 < 6);
UKB_Lifestyle_0_0.x1160_S(isnan(UKB_Lifestyle_0_0.x1160)) = NaN;
UKB_Lifestyle_0_0.x1160_L(isnan(UKB_Lifestyle_0_0.x1160)) = NaN;
UKB_Lifestyle_0_0 = removevars(UKB_Lifestyle_0_0, 'x1160');
x709 = UKB_Lifestyle_0_0.x709;
UKB_Lifestyle_0_0.x709 = double(UKB_Lifestyle_0_0.x709 == 1);
UKB_Lifestyle_0_0.x709(isnan(x709)) = NaN;

brain_data_360p = [ukb_RedSynYeo7n_2_0(:,1), array2table(Red360p_NC),array2table(Syn360p_NC)];
Red7n_NC = netMean1D(brain_data_360p{:,2:361},Yeo7MMP);
Syn7n_NC = netMean1D(brain_data_360p{:,362:721},Yeo7MMP);
brain_data_7n = [brain_data_360p(:,1),array2table(Red7n_NC),array2table(Syn7n_NC)];
brain_data_7n_matrix = ukb_RedSynYeo7n_2_0;

clear Red360p_NC Red7n_NC results_7n results_lifestyle_7nx results_pheno_adult results_pheno_older ...
    Syn360p_NC Syn7n_NC T T_ab ukb_SynYeo7n_2_0 ukb_RedYeo7n_2_0 ukb_RedSynYeo7n_2_0 modal_pheno_7n modal_pheno

%% Demographic Variables Asociation

% 360 regions * age
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; 

formula = ['Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_360p = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_360p = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_360p = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_360p = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);

% 7*7 network-modules
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n_matrix);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; 

formula = ['Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_7n_matrix = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_7n_matrix = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_7n_matrix = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_7n_matrix = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);

% 7 networks
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; 

formula = ['Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'AgeAttend_2_0'};
results_age_7n = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'Sex_0_0'};
results_sex_7n = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = {'BMI_2_0'};
results_bmi_7n = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+(1|Centre_0_0)'];
XList = {'Vol_WB_TIV_2_0'};
results_tiv_7n = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);


results_demo_360p = [results_age_360p;results_sex_360p;results_bmi_360p;results_tiv_360p];
results_demo_7n_matrix = [results_age_7n_matrix;results_sex_7n_matrix;results_bmi_7n_matrix;results_tiv_7n_matrix];
results_demo_7n = [results_age_7n;results_sex_7n;results_bmi_7n;results_tiv_7n];


%% Phenotype Association
formula = ['Sex_0_0+AgeAttend_2_0+BMI_2_0+AgeAttend_2_0^2+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];

% 360 regions
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_360p);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_360p.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = {'x20016','x4526','x2178'}';
results_pheno_360p = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
results_pheno_adult_360p = LMM_withFormula2(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0<65,:),formula,XList,YList);
results_pheno_older_360p = LMM_withFormula2(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0>=65,:),formula,XList,YList);

% 7*7 network-modules
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n_matrix);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
results_pheno_7n_matrix = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
results_pheno_adult_7n_matrix = LMM_withFormula2(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0<65,:),formula,XList,YList);
results_pheno_older_7n_matrix = LMM_withFormula2(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0>=65,:),formula,XList,YList);

% 7 networks
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';

XList = UKB_Outcomes_2_0.Properties.VariableNames(2:end)';
results_pheno_7n = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
results_pheno_adult_7n = LMM_withFormula2(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0<65,:),formula,XList,YList);
results_pheno_older_7n = LMM_withFormula2(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0>=65,:),formula,XList,YList);


%% Envrionmental & Lifstyle Association

formula = ['Sex_0_0+AgeAttend_2_0+AgeAttend_2_0^2+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];

% 360 regions
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_360p);
% UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
XList = {'x1160_L','x1160_S','x1558','x1568','x1588','x20489','x20497','x709','x738'};
YList = brain_data_360p.Properties.VariableNames(2:end)';
results_lifestyle_360p = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
% writetable(results_lifestyle_360p,'modal_lifestyle_7n.txt','Delimiter' ,' ')

% 7*7 network-modules
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n_matrix);
% UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = UKB_Lifestyle_0_0.Properties.VariableNames(95:end)';
results_lifestyle_7n_matrix = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
% writetable(results_lifestyle_7n_matrix,'modal_lifestyle_7n_matrix.txt','Delimiter' ,' ')

% 7 networks
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n);
% UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
XList = UKB_Lifestyle_0_0.Properties.VariableNames(95:end)';
results_lifestyle_7n = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
results_lifestyle_7n(isnan(results_lifestyle_7n.t),:) = [];
% writetable(results_lifestyle_7n,'modal_lifestyle_7n_matrix.txt','Delimiter' ,' ')


% UKB_BrainXY.SA2Y = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2Y(UKB_BrainXY.AgeAttend_2_0 < 65) = 0;
% UKB_BrainXY.SA2Y((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 7)) = 1;
% 
% UKB_BrainXY.SA2NA = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 <= 6)) = 0;
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 6)) = 1;
% 
% XList = {'SA2Y','SA2NA'}';
% results_pheno_360p_SA = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
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
% results_pheno_7n_matrix_SA = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
% UKB_BrainXY.SA2Y = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2Y(UKB_BrainXY.AgeAttend_2_0 < 65) = 0;
% UKB_BrainXY.SA2Y((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 7)) = 1;
% 
% UKB_BrainXY.SA2NA = nan(size(UKB_BrainXY,1),1);
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 <= 6)) = 0;
% UKB_BrainXY.SA2NA((UKB_BrainXY.AgeAttend_2_0 >= 65)&(UKB_BrainXY.x20016 > 6)) = 1;
% 
% XList = {'SA2Y','SA2NA'}';
% results_pheno_7n_SA = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);

%% Biomarker Association

UKB_Biomarker = readtable('/public/home/zhangjie/ZJLab/UKBiobank_Project/data/table/UKB_BiomarkerMRI_0_0.csv','Delimiter' ,',','TreatAsEmpty','NA');

% 7 networks
UKB_BrainXY = innerjoin(UKB_Biomarker,brain_data_7n);
% UKB_BrainXY.AgeGap = UKB_BrainXY.AgeAttend_2_0 - UKB_BrainXY.AgeAttend_0_0;
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n.Properties.VariableNames(2:end)'; 

formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = UKB_Biomarker.Properties.VariableNames(95:end)'; 
results_biomarker_7n = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);

% 7x7 network-modules
UKB_BrainXY = innerjoin(UKB_Biomarker,brain_data_7n_matrix);
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = brain_data_7n_matrix.Properties.VariableNames(2:end)'; 

formula = ['AgeAttend_2_0+AgeAttend_2_0^2+Sex_0_0+BMI_2_0+HeadMotion_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = UKB_Biomarker.Properties.VariableNames(95:end)'; 
results_biomarker_7n_matrix = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);

writetable(results_biomarker_7n,'results_biomarker_7n.txt','Delimiter' ,',')
writetable(results_biomarker_7n_matrix,'results_biomarker_7n_matrix.txt','Delimiter' ,',')

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
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n);
YList = {'x20016','x4526','x2178'};
XList = {'x1160_L','x1160_S','x1558','x1568','x1588','x20489','x20497','x709','x738'};
MList = brain_data_7n.Properties.VariableNames(2:end)'; % Y.Properties.RedSyniableNames(2:end)';
Centre = [AsDummy(UKB_BrainXY.Centre_0_0)];
Centre(:,[sum(Centre)<1000]) = [];
COVList = {'AgeAttend_2_0','Sex_0_0','HeadMotion_2_0',...
    'BMI_2_0','Race_1','Race_2','Race_3','Race_4','Vol_WB_TIV_2_0'};
COV0 = [UKB_BrainXY{:,COVList},Centre];
for m = 1:length(MList)
    for i = 1:length(XList)
        parfor j = 1:length(YList)
            % clear MediationRes
            % === Data extraction ===
            Xvar = UKB_BrainXY{:, XList{i}};Yvar = UKB_BrainXY{:, YList{j}};Mvar = UKB_BrainXY{:, MList{m}};
            
            % Remove rows with NaN in any variable
            [X2, Y2, M2, COV2] = nanrm(Xvar, Yvar, Mvar, COV0);
            
            % === Step 1. Remove near-constant covariates ===
            if ~isempty(COV2)
                covVar = var(COV2, 0, 1);
                constCols = covVar < 1e-8;
                if any(constCols)
                    fprintf('Removing %d near-constant covariates\n', sum(constCols));
                    COV2(:, constCols) = [];
                end
            end
            
            % === Step 2. Z-score normalization for all numeric variables ===
            % (important for numerical stability and comparability)
            X2 = zscore(X2);Y2 = zscore(Y2);M2 = zscore(M2);
            if ~isempty(COV2)
                COV2 = zscore(COV2);
            end
            
            % === Step 3. Check condition number of design matrix (optional) ===
            if ~isempty(COV2)
                condVal = cond([ones(size(COV2,1),1), COV2]);
                if condVal > 1e12
                    warning('High condition number in covariates: %g (possible collinearity)', condVal);
                end
            end
            
            % === Step 4. Run mediation with stable parameters ===
            try
                [paths, stats] = mediation(X2, Y2, M2, ...
                    'covs', COV2, ...% 'boot', 'bootsamples', 1000, ...
                    'verbose', 'no');  % suppress unnecessary text output
            catch ME
                warning('Mediation failed at i=%d, j=%d, m=%d: %s', i, j, m, ME.message);
                paths = []; stats = [];
            end
            MediateResult = [stats.paths;stats.mean;stats.ste;stats.t;stats.df;stats.p];
            MediateResults(:,:,i,j,m) = MediateResult;
            disp(['Mediate Analysis: X=',XList{i},' ','Y=',YList{j},' M=',MList{m}]);
        end
    end
end
MediateResults(:,:,:,3,2)
x = squeeze(MediateResults(:,5,:,3,2));
x = squeeze(MediateResults(6,5,:,3,:));


Yeo7Names2 = {'VN','SMN','DAN','VAN','LN','FPN','DMN'};
MList = {'red-VN','red-SMN','red-DAN','red-VAN','red-LN','red-FPN','red-DMN',...
    'syn-VN','syn-SMN','syn-DAN','syn-VAN','syn-LN','syn-FPN','syn-DMN'}

StatNames = {'path','mean','ste','t','df','p'};
PathList = {'a','b','c','c''','ab'};
numRows = length(XList)*length(YList)*length(MList)*length(PathList);

X_col = cell(numRows,1);
Y_col = cell(numRows,1);
M_col = cell(numRows,1);
Path_col = cell(numRows,1);
StatCols = nan(numRows,6);  

idx = 1;
for xi = 1:length(XList)
    for yi = 1:length(YList)
        for mi = 1:length(MList)
            for pi = 1:length(PathList)
                X_col{idx} = XList{xi};
                Y_col{idx} = YList{yi};
                M_col{idx} = MList{mi};
                Path_col{idx} = PathList{pi};
                StatCols(idx,:) = MediateResults(:,pi,xi,yi,mi)'; 
                idx = idx + 1;
            end
        end
    end
end

T = table(X_col, M_col, Y_col, Path_col, StatCols(:,1), StatCols(:,2), StatCols(:,3), StatCols(:,4), StatCols(:,5), StatCols(:,6), ...
    'VariableNames', {'X','M','Y','Path', StatNames{:}});

head(T)

writetable(T,'/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/mediation_results_7nx.txt')

%% Interaction Effect
UKB_BrainXY = innerjoin(innerjoin(UKB_Lifestyle_0_0,UKB_Outcomes_2_0),brain_data_7n);
YList = brain_data_7n.Properties.VariableNames(2:end)';
allResults = table();

for i = 1:14
    residualModel = fitlme(UKB_BrainXY, ...
        [YList{i}, ' ~ Sex_0_0 + BMI_2_0 + HeadMotion_2_0 + Race_1 + Race_2 + Race_3 + Race_4 + Vol_WB_TIV_2_0 + (1|Centre_0_0)']);
    UKB_BrainXY.Residual = residuals(residualModel);
    interactionModel = fitlm(UKB_BrainXY, ...
        'x20016 ~ Residual + AgeAttend_2_0 + AgeAttend_2_0 * Residual');
    coefTable = interactionModel.Coefficients;
    coefTable.YVariable = repmat(YList{i}, height(coefTable), 1);
    estimate = coefTable.Estimate;
    se = coefTable.SE;
    tStat = coefTable.tStat;
    pValue = coefTable.pValue;
    rowNames = {'(Intercept)';'AgeAttend_2_0';'Residual';'AgeAttend_2_0:Residual'};  % »ñÈ¡ÐÐÃû (¼´ÏµÊýµÄÃû³Æ)
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

results_APOE4_7n = LMM_withFormula2(UKB_BrainXY,formula,XList,YList);
results_APOE4_adult_7n = LMM_withFormula2(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0<65,:),formula,XList,YList);
results_APOE4_older_7n = LMM_withFormula2(UKB_BrainXY(UKB_BrainXY.AgeAttend_2_0>=65,:),formula,XList,YList);



