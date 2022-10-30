%% single subject maps
%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paths
base_dir = '/Volumes/Data/zoocon/Rev_Hippocampgoal/';
out_dir = [base_dir, 'Stats/'];
data_path = [base_dir, 'Data/'];
scripts_path = [base_dir, 'Scripts/'];
color_path = [base_dir, 'Colormaps/']; % for python colormaps
plots_path = [base_dir, 'Figures/'];

% env
addpath(scripts_path);
addpath(color_path);
% load data 
% load('/Volumes/Data/zoocon/Rev_Hippocampgoal/Stats/HIPP_MERGE_BL/HIPP_MERGE_BL_stats.mat')
load([out_dir, 'HIPP_MERGE_BL/HIPP_MERGE_BL_stats.mat'])

%%%%%%%%%%%%%%%%%%%%%%%% Vis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask = all_stats(1).stats.pos_clust.pixels;
nsubs = size(all_stats(1).stats.data1,1);
% preallocate
out_mat = zeros([nsubs,1]);
for isub = 1:nsubs
    out_xlabels{isub}= ['sub',num2str(isub, '%02d')];
    tmp1 = squeeze(all_stats(1).stats.data1(isub,:,:) - all_stats(1).stats.data2(isub,:,:));
    
    out_mat(isub) = mean(tmp1(mask)); % subset
end

%% sign test
[p,h,stats]= signrank(out_mat)

%% Figure
% Single subject maps
f1 = figure;
figure;bar(categorical(out_xlabels), out_mat);
ax=gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ax.YLabel.String = 'Average Similarity';
title('Single Subject Similarity at Significant Timepoints')
saveas(f1, [plots_path, 'single_subjects_hippocmpal_cluster_mass.eps'], 'epsc');
close

% 
% ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
% ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
% set(ax,'TickDir','out','box','off')
% colormap(viridis)
% c = colorbar;
% set(c, 'FontSize',15)
% title('analytic zmap')
% % saveas(f1, [plots_path, 'Fig4_group_level_con_div_zmap.eps'], 'epsc');
% close 