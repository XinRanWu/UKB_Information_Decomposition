%% Part1 left_right
cd('E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age\score_table_lower\left_right');
% Traverse all subdirectories
subDirs = dir('*.*'); 
subDirs(1) = []; % Remove '.' from the current directory
subDirs(1) = []; % Remove '..'from the upper-level directory

for i = 1:length(subDirs)
    if subDirs(i).isdir
        disp(i);
        % If it is a directory, iterate recursively
        cd(subDirs(i).name); % Go to subdirectory
        csvFiles = dir('*.csv'); % Only the CSV file in the current directory is found
        for j = 1:length(csvFiles)
            csvFileName = csvFiles(j).name; 
            name = csvFileName(1:end-4);
            tablePath = ['E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age\score_table_lower\left_right','\',subDirs(i).name,'\',csvFileName];
            table = readtable(csvFileName);
            Pval = table.PVal;
            Pfdr = multicmp(Pval,'fdr');
            table.Pfdr = Pfdr;
            writetable(table,tablePath);
        end
        cd('..'); % Returns the upper-level directory
    end
end


%% Part2 extra_intra
cd('E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age\score_table_lower\extra_intra');
% Traverse all subdirectories
subDirs = dir('*.*'); 
subDirs(1) = []; % Remove '.' from the current directory
subDirs(1) = []; % Remove '..'from the upper-level directory

for i = 1:length(subDirs)
    if subDirs(i).isdir
        disp(i);
        % If it is a directory, iterate recursively
        cd(subDirs(i).name); % Go to subdirectory
        csvFiles = dir('*.csv'); % Only the CSV file in the current directory is found
        for j = 1:length(csvFiles)
            csvFileName = csvFiles(j).name; 
            name = csvFileName(1:end-4);
            tablePath = ['E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age\score_table_lower\extra_intra','\',subDirs(i).name,'\',csvFileName];
            table = readtable(csvFileName);
            Pval = table.PVal;
            Pfdr = multicmp(Pval,'fdr');
            table.Pfdr = Pfdr;
            writetable(table,tablePath);
        end
        cd('..'); % Returns the upper-level directory
    end
end


%% Part3 homotopy
cd('E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age\score_table_lower\homotopy');
% Traverse all subdirectories
subDirs = dir('*.*'); 
subDirs(1) = []; % Remove '.' from the current directory
subDirs(1) = []; % Remove '..'from the upper-level directory

for i = 1:length(subDirs)
    if subDirs(i).isdir
        disp(i);
        % If it is a directory, iterate recursively
        cd(subDirs(i).name); % Go to subdirectory
        csvFiles = dir('*.csv'); % Only the CSV file in the current directory is found
        for j = 1:length(csvFiles)
            csvFileName = csvFiles(j).name; 
            name = csvFileName(1:end-4);
            tablePath = ['E:\Master\Brain_project\UKB_02_PhiID_glasser\result\score_table_split_age\score_table_lower\homotopy','\',subDirs(i).name,'\',csvFileName];
            table = readtable(csvFileName);
            Pval = table.PVal;
            Pfdr = multicmp(Pval,'fdr');
            table.Pfdr = Pfdr;
            writetable(table,tablePath);
        end
        cd('..'); % Returns the upper-level directory
    end
end



