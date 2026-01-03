cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/data/table');
UKB_UsedData_Brain = readtable('UKB_UsedData_Brain.csv');

%% Part1 left_right
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/left_right');

csvFiles = dir('*.csv');
left_right = cell(size(csvFiles,1),2);
for i = 1:length(csvFiles)
    disp(i);
    csvFileName = csvFiles(i).name; 
    tableData = readtable(csvFileName, 'Delimiter', ',');
    Data = innerjoin(UKB_UsedData_Brain,tableData(:,[1,8:end]));
    left_right{i,1} = csvFileName;
    column_name = tableData.Properties.VariableNames(:,8:end);
    value_cell = get_age_p_tvalue_lmm(Data,column_name);
    left_right{i,2} = value_cell;
    clear csvFileName tableData column_name pvalue_cell
end

%% Part2 extra_intra
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/extra_intra');

csvFiles = dir('*.csv');
extra_intra = cell(size(csvFiles,1),2);
for i = 1:length(csvFiles)
    disp(i);
    csvFileName = csvFiles(i).name; 
    tableData = readtable(csvFileName, 'Delimiter', ',');
    Data = innerjoin(UKB_UsedData_Brain,tableData(:,[1,8:end]));
    extra_intra{i,1} = csvFileName;
    column_name = tableData.Properties.VariableNames(:,8:end);
    value_cell = get_age_p_tvalue_lmm(Data,column_name);
    extra_intra{i,2} = value_cell;
    clear csvFileName tableData column_name pvalue_cell
end


%% Part3 homotopy
cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/table_allsubject/homotopy');

csvFiles = dir('*.csv');
homotopy = cell(size(csvFiles,1),2);
for i = 1:length(csvFiles)
    disp(i);
    csvFileName = csvFiles(i).name; 
    tableData = readtable(csvFileName, 'Delimiter', ',');
    Data = innerjoin(UKB_UsedData_Brain,tableData(:,[1,8:end]));
    homotopy{i,1} = csvFileName;
    column_name = tableData.Properties.VariableNames(:,8:end);
    value_cell = get_age_p_tvalue_lmm(Data,column_name);
    homotopy{i,2} = value_cell;
    clear csvFileName tableData column_name pvalue_cell
end

cd('/public/home/yangsy/Brain/Brain_project/UKB_02_PhiID_glasser/result/age_effect');
save('lmm_age_p_tvalue_allsubject_without_agesquare.mat','extra_intra','left_right','homotopy','-v7.3');


