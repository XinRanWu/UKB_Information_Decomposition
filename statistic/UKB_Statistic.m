%% load data
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/data/table');
UKB_UsedData_Brain = readtable('UKB_UsedData_Brain.csv');
UKB_BrainImagingVars = readtable('UKB_BrainImagingVars.csv');
UKB_AllVars = readtable('UKB_AllVars_NumCol.csv');
UKB_MentalDiagnosis = readtable('UKB_MentalDiagnosis.csv');


BrainTotalList = {'FA_WB_Total_2_0','MD_WB_Total_2_0','MO_WB_Total_2_0','L1_WB_Total_2_0','L2_WB_Total_2_0',...
    'L3_WB_Total_2_0','ICVF_WB_Total_2_0','ISOVF_WB_Total_2_0','Area_WB_Total_2_0','Thickness_WB_Total_2_0',...
    'WMHyperInt_WB_Total_2_0','Vol_WB_SCG_2_0','Vol_WB_CSF_2_0','Vol_WB_Vent3_2_0','Vol_WB_Vent4_2_0',...
    'Vol_WB_Vent5_2_0','Vol_WB_Cort_2_0','Vol_WB_CWM_2_0'}';

MentalList = UKB_MentalDiagnosis.Properties.VariableNames(2:end)';
BrainList = UKB_BrainImagingVars.Properties.VariableNames(2:end)';
BrainAsegList = BrainList((contains(BrainList,'Vol'))&(contains(BrainList,'_2_0')));
BrainAsegList = BrainAsegList([31:38,47:54]);

%compute 
formula = ['Sex_0_0+AgeAttend_2_0+HeadMotion_2_0+BMI_2_0+',...
    'Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)'];
XList = UKB_AllVars.Properties.VariableNames(3:end)';
XList = XList(contains(XList,'_2_'));
XList = XList(~(contains(XList,'_2_0_1')|contains(XList,'Intensity_')|...
    contains(XList,'VolRatio_')|contains(XList,'HBF_')|contains(XList,'Vol_WB_TIV_2_0')));
XList = XList(~(contains(XList,'_0_0')|contains(XList,'_1_0')|contains(XList,'_3_0')));
XListAll = [{'TDI_0_0';'Income_2_0';'EduAge_0_0';'Handedness_0_0'};...
    XList;BrainTotalList;BrainAsegList;MentalList(sum(UKB_UsedData_Brain{:,1885:1900},1)>200)];

clear BrainAsegList BrainList BrainTotalList MentalList XList
clear UKB_AllVars UKB_BrainImagingVars UKB_BrainMRICov UKB_MentalDiagnosis UKB_rfMRICov


dt = datetime('now', 'TimeZone', 'Asia/Shanghai');
parpool('local', 20);
%% Part1 left_right
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/left_right');
csvFiles = dir('*.csv');
for i = 1:length(csvFiles)
    cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/left_right');
    disp(i);
    csvFileName = csvFiles(i).name; 
    name = csvFileName(1:end-4);
    tablePath = ['/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/score_table_allsubject/left_right','/',name];
    mkdir(tablePath);
    ukb_RedSyn_table = readtable(csvFileName);
    get_ukb_statistic_result(ukb_RedSyn_table,UKB_UsedData_Brain,formula,XListAll,tablePath)
end

%% Part2 extra_intra
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/extra_intra');
csvFiles = dir('*.csv');
for i = 1:length(csvFiles)
    cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/extra_intra');
    disp(i);
    csvFileName = csvFiles(i).name; 
    name = csvFileName(1:end-4);
    tablePath = ['/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/score_table_allsubject/extra_intra','/',name];
    mkdir(tablePath);
    ukb_RedSyn_table = readtable(csvFileName);
    get_ukb_statistic_result(ukb_RedSyn_table,UKB_UsedData_Brain,formula,XListAll,tablePath)
end


%% Part3 homotopy
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/homotopy');
csvFiles = dir('*.csv');
for i = 1:length(csvFiles)
    cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/homotopy');
    disp(i);
    csvFileName = csvFiles(i).name; 
    name = csvFileName(1:end-4);
    tablePath = ['/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/score_table_allsubject/homotopy','/',name];
    mkdir(tablePath);
    ukb_RedSyn_table = readtable(csvFileName);
    get_ukb_statistic_result(ukb_RedSyn_table,UKB_UsedData_Brain,formula,XListAll,tablePath)
end

