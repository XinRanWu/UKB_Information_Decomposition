% resVec = [diag(X);icatb_mat2vec(X)'];
function plot_matrix_Yeo7n(resVec,clim)
figure;
load('/public/home/zhangjie/DataAnalysis/wxr_toolbox/toolbox_matlab/cbrewer2-master/cbrewer2/colorbrewer.mat', 'colorbrewer');
Yeo7Names = {'Visual';'Somatomotor';'DorsalAttention';'VentralAttention';'Limbic';'Frontoparietal';'Default'};
imagesc(resVec,[-1*clim,clim]);set(gca,'DataAspectRatio',[1 1 1]);colorbar;colormap(mkcmap(colorbrewer.div.RdBu{1,11}(end:-1:1,:)))
set(gca, 'XTick', 1:7, 'XTickLabel', Yeo7Names);set(gca, 'YTick', 1:7, 'YTickLabel', Yeo7Names);xtickangle(90)
end