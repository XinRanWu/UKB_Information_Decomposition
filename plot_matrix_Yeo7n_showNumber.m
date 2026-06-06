% resVec = [diag(X);icatb_mat2vec(X)'];
function plot_matrix_Yeo7n_showNumber(resVec,clim)
figure;
load('/public/home/zhangjie/DataAnalysis/wxr_toolbox/toolbox_matlab/cbrewer2-master/cbrewer2/colorbrewer.mat', 'colorbrewer');
Yeo7Names = {'Visual';'Somatomotor';'DorsalAttention';'VentralAttention';'Limbic';'Frontoparietal';'Default'};
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
h.ColorLimits = [-1*clim, clim];
end