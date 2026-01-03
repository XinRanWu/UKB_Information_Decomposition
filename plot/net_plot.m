%% Part1 
clear;
cd('C:\Users\AS\Desktop\PID\data');
load("lmm_age_p_tvalue_allsubject_without_agesquare.mat")

syn_line = extra_intra{3,2};
syn_line_t = [syn_line{:,3}];
syn_line_p = [syn_line{:,2}];

syn_line_p = multicmp(syn_line_p,'fdr');
syn_line_t(syn_line_p > 0.05) = nan;

syn_net = array2network(syn_line_t,1);

red_line = extra_intra{2,2};
red_line_t = [red_line{:,3}];
red_line_p = [red_line{:,2}];

red_line_p = multicmp(red_line_p,'fdr');
red_line_t(red_line_p > 0.05) = nan;

red_net = array2network(red_line_t,1);

yeo7name = {'VIS','SMN','DAN','SAN','limbic','FPN','DMN'};

figure(1)
subplot(1,2,1)
h=heatmap(red_net);
title('red')
% 设置颜色条为'turbo'
colormap(turbo);
% 获取数据的最大值和最小值
max_val = max(red_net(:));
min_val = min(red_net(:));
% 计算颜色条的范围，使得0位于中央
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
% 显示颜色条
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = yeo7name;
h.YDisplayLabels = yeo7name;

subplot(1,2,2)
h=heatmap(syn_net);
title('syn')
% 设置颜色条为'turbo'
colormap(turbo);
% 获取数据的最大值和最小值
max_val = max(syn_net(:));
min_val = min(syn_net(:));
% 计算颜色条的范围，使得0位于中央
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
% 显示颜色条
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = yeo7name;
h.YDisplayLabels = yeo7name;


%% Part 2
cd('C:\Users\AS\Desktop\PID\data');
load("lmm_age_p_tvalue_allsubject_without_agesquare.mat")
idx = [7;8;6;2;5;4;3;1];
datacell = left_right(idx,:);
id_red = [1;3;5;7];
id_syn = [2;4;6;8];

red_p = zeros(8,4);
red_t = zeros(8,4);
syn_p = zeros(8,4);
syn_t = zeros(8,4);
for i = 1:size(datacell,1)
    disp(i);
    data = datacell{i,2};
    data = data(:,2:3);
    data = cell2mat(data);
    red_p(i,:) = data(id_red,1)';
    red_t(i,:) = data(id_red,2)';
    syn_p(i,:) = data(id_syn,1)';
    syn_t(i,:) = data(id_syn,2)';
end

clear idx id_red id_syn data i

red_p(red_p > 0.05) = nan;
id_nan = isnan(red_p);
red_t(id_nan) = nan;

syn_p(syn_p > 0.05) = nan;
id_nan = isnan(syn_p);
syn_t(id_nan) = nan;

columnname = {'TOTAL','VIS','SMN','DAN','SAN','limbic','FPN','DMN'}';
rowname = {'R','L','R+L','R-L'};

figure(1)
subplot(1,2,1)
h=heatmap(red_t);
title('Redundancy t value')
% 设置颜色条为'turbo'
colormap(turbo);
% 获取数据的最大值和最小值
max_val = max(red_t(:));
min_val = min(red_t(:));
% 计算颜色条的范围，使得0位于中央
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
% 显示颜色条
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = rowname;
h.YDisplayLabels = columnname;

subplot(1,2,2)
h=heatmap(syn_t);
title('Synergy t value')
% 设置颜色条为'turbo'
colormap(turbo);
% 获取数据的最大值和最小值
max_val = max(syn_t(:));
min_val = min(syn_t(:));
% 计算颜色条的范围，使得0位于中央
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
% 显示颜色条
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = rowname;
h.YDisplayLabels = columnname;

%% Part 3
cd('C:\Users\AS\Desktop\PID\data');
load("lmm_age_p_tvalue_allsubject_without_agesquare.mat")

extra_net = extra_intra{1,2};
extra_red = cell2mat(extra_net(1:7,2:3));
extra_syn = cell2mat(extra_net(8:14,2:3));
extra_p = zeros(7,2);
extra_t = zeros(7,2);
extra_p(1:7,1) = extra_red(1:7,1); 
extra_p(1:7,2) = extra_syn(1:7,1);
extra_t(1:7,1) = extra_red(1:7,2);
extra_t(1:7,2) = extra_syn(1:7,2);

extra_p(extra_p > 0.05) = nan;
id_nan = isnan(extra_p);
extra_t(id_nan) = nan;

rowname = {'Redundancy','Synergy'};
columnname = {'VIS','SMN','DAN','SAN','limbic','FPN','DMN'}';

figure(1)
h=heatmap(extra_t);
title('Extra net t value')
% 设置颜色条为'turbo'
colormap(turbo);
% 获取数据的最大值和最小值
max_val = max(extra_t(:));
min_val = min(extra_t(:));
% 计算颜色条的范围，使得0位于中央
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
% 显示颜色条
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = rowname;
h.YDisplayLabels = columnname;

%% Part 4
cd('C:\Users\AS\Desktop\PID\data');
load("lmm_age_p_tvalue_allsubject_without_agesquare.mat")
total = extra_intra{4,2};
total_p = zeros(2,2);
total_t = zeros(2,2);

total_p(1,1) = total{1,2};
total_p(2,1) = total{3,2};
total_p(1,2) = total{2,2};
total_p(2,2) = total{4,2};

total_t(1,1) = total{1,3};
total_t(2,1) = total{3,3};
total_t(1,2) = total{2,3};
total_t(2,2) = total{4,3};

total_p(total_p > 0.05) = nan;
id_nan = isnan(total_p);
total_t(id_nan) = nan;

rowname = {'Redundancy','Synergy'};
columnname = {'Intra net','Extra net'}';

figure(1)
h=heatmap(total_t);
title('Total net t value')
% 设置颜色条为'turbo'
colormap(turbo);
% 获取数据的最大值和最小值
max_val = max(total_t(:));
min_val = min(total_t(:));
% 计算颜色条的范围，使得0位于中央
if max_val > abs(min_val)
    clim([-max_val, max_val]);
else
    clim([min_val, -min_val]);
end
% 显示颜色条
colorbar;
h.MissingDataColor = 'w';
h.XDisplayLabels = rowname;
h.YDisplayLabels = columnname;

