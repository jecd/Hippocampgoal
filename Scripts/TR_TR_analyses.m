%% TR by TR analysis figures 
% code to reproduce TR * TR analysis and figures

%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% paths
base_dir = '/Volumes/Data/zoocon/Hippocampgoal/';
out_dir = [base_dir, 'Stats/'];
data_path = [base_dir, 'Data/'];
scripts_path = [base_dir, 'Scripts/'];
color_path = [base_dir, 'Colormaps/']; % for python colormaps
plots_path = [base_dir, 'Figures/'];

% env
addpath(scripts_path);
addpath(color_path);

% load data 
load([data_path, '/TR_TR_RDMs.mat'])

%% do stats
% this function will take in agregated subject data and do comparisons based on user supplied contrasts
% loops over ROIs

% add interaction contrast
ROI_names = fieldnames(RDMs);
contrast_names = fieldnames(RDMs.(ROI_names{1}).mean_data);
% loop through ROIs and add the contrast
for iROI = 1:length(ROI_names)
    cur_roi = ROI_names{iROI}; % name for index
    
    % names for index 
    cName1a = 'diff_start_same_end_scon_cue'; % converging same con
    cName1b = 'same_start_diff_end_scon_cue'; % diverging same con
    cName2a = 'diff_start_same_end_dcon_cue'; % converging diff con
    cName2b = 'same_start_diff_end_dcon_cue'; % diverging diff con
    
    data1a = RDMs.(cur_roi).mean_data.(cName1a).indiv;
    data1b = RDMs.(cur_roi).mean_data.(cName1b).indiv;
    data2a = RDMs.(cur_roi).mean_data.(cName2a).indiv;
    data2b = RDMs.(cur_roi).mean_data.(cName2b).indiv;
    
    % converging - diverging scon
    z = data1a - data1b;
    out_struct = [];
    out_struct.data = squeeze(mean(z,1)); % mean data 
    out_struct.indiv = z; % indiv
    out_struct.dimord = 'nsub*tr*tr';
    RDMs.HIPP_MERGE_BL.mean_data.con_div_scon = out_struct;
    
    % converging - diverging dcon
    z = data2a - data2b;
    out_struct = [];
    out_struct.data = squeeze(mean(z,1)); % mean data 
    out_struct.indiv = z; % indiv
    out_struct.dimord = 'nsub*tr*tr';
    RDMs.HIPP_MERGE_BL.mean_data.con_div_dcon = out_struct;
end

contrasts = {
            {'diff_start_same_end_scon_cue', 'same_start_diff_end_scon_cue'}, % converging > diverging
            {'diff_start_same_end_scon_cue', 'diff_start_same_end_dcon_cue'}, % converging scon > converging dcon
            {'same_start_diff_end_scon_cue', 'same_start_diff_end_dcon_cue'}, % diverging scon > diverging dcon
            {'diff_start_same_end_dcon_cue', 'same_start_diff_end_dcon_cue'},  % converging dcon > diverging dcon
            {'con_div_scon', 'con_div_dcon'}
            }; 

tc_stats_xcor_FIR(RDMs, contrasts, out_dir); % ~ 7min with 10k perms

%%%%%%%%%%%%%%%%%%%%%%%% Main Paper Figs %%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4 

cur_ROI = 'HIPP_MERGE_BL';
load([out_dir, cur_ROI,'/', cur_ROI, '_stats.mat'])
interest_contrast = {{'diff_start_same_end_scon_cue', 'same_start_diff_end_scon_cue'}};
data_struct = all_stats.stats; %% FIX ME 

% apply window and figure out masking
manual_lag = 10; % in TRs % This makes P1 activation the zeropoint on the figure
tmp_mask = zeros(size(data_struct.zmap, 1), size(data_struct.zmap,2));
tmp_mask(manual_lag:end, manual_lag:end) = 1;
tmp_mask = logical(tmp_mask);
out_mat_dims = sqrt(length(find(tmp_mask)));

% loop through grabbing the data and replacing
tmp_data1_mat = zeros([size(data_struct.data1,1), out_mat_dims ,out_mat_dims]);
tmp_data2_mat = zeros([size(data_struct.data2,1), out_mat_dims, out_mat_dims]);
for isub = 1:size(data_struct.data1,1)
    % c1
    tmp_data = squeeze(data_struct.data1(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data1_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
    % c2
    tmp_data = squeeze(data_struct.data2(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data2_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
end
data_struct.data1 = tmp_data1_mat; % subset back into data
data_struct.data2 = tmp_data2_mat;

%% make group level plots for main figure
ref_keep = ones(size(data_struct.zmap,1), size(data_struct.zmap,2));
ref_keep = tril(ref_keep);
ref_keep = ref_keep(manual_lag:end, manual_lag:end);

% con - div
f1 = figure;
d1_h = imagesc(reshape(data_struct.zmap(tmp_mask),  [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('analytic zmap')
saveas(f1, [plots_path, 'Fig4_group_level_con_div_zmap.eps'], 'epsc');
close 

% thresh
d = reshape(data_struct.zmapthresh(tmp_mask), [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]);
tmp = bwconncomp(d);
tmp_idx = cellfun(@(x) length(x)==1, tmp.PixelIdxList);
d(cell2mat(tmp.PixelIdxList(tmp_idx))) = 0;
f1=figure;d1_h=imagesc(d);
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
title('zmap thresh')
colormap(viridis)
caxis manual
caxis([0,3])
c = colorbar;
set(c, 'FontSize',15)
saveas(f1, [plots_path, 'Fig4_group_level_con_div_zmap_thresh.eps'], 'epsc');
close 

% con
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data1,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
caxis manual
caxis([0,0.035])
c = colorbar;
set(c, 'FontSize',15)
title('converging')
saveas(f1, [plots_path, 'Fig4_group_level_converging.eps'], 'epsc');
close 

% div 
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data2,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
caxis manual
caxis([0,0.035])
c = colorbar;
set(c, 'FontSize',15)
title('diverging')
saveas(f1, [plots_path, 'Fig4_group_level_diverging.eps'], 'epsc');
close 

% localize zones
zz = reshape(data_struct.zmap(tmp_mask), [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]);
f1 = figure('Position', [10,10,1200,500]); hold on
plot(mean(zz(:,1:5),2), 'LineWidth', 2)
xline(1, '--', 'LineWidth',1)
xline(6, '--', 'LineWidth',1)
xline(11, '--', 'LineWidth',1)
xline(16, '--', 'LineWidth',1)
xline(21, '--', 'LineWidth',1)
xline(26, '--', 'LineWidth',1)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ylabel('Mean Pattern Similarity (z)', 'FontSize', 15)
title('P1 w/ All')
saveas(f1, [plots_path, 'Fig4_p1_w_all.eps'], 'epsc');
close 

%% reviwer figures
cur_ROI = 'HIPP_MERGE_BL';
load([out_dir, cur_ROI,'/', cur_ROI, '_stats.mat'])
% contrasts = {{'diff_start_same_end_scon_cue',% 'same_start_diff_end_scon_cue'}}; not use
data_struct = all_stats(5).stats;

% apply window and figure out masking
manual_lag = 10; % in TRs % This makes P1 activation the zeropoint on the figure
tmp_mask = zeros(size(data_struct.zmap, 1), size(data_struct.zmap,2));
tmp_mask(manual_lag:end, manual_lag:end) = 1;
tmp_mask = logical(tmp_mask);
out_mat_dims = sqrt(length(find(tmp_mask)));

% loop through grabbing the data and replacing
tmp_data1_mat = zeros([size(data_struct.data1,1), out_mat_dims ,out_mat_dims]);
tmp_data2_mat = zeros([size(data_struct.data2,1), out_mat_dims, out_mat_dims]);
for isub = 1:size(data_struct.data1,1)
    % c1
    tmp_data = squeeze(data_struct.data1(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data1_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
    % c2
    tmp_data = squeeze(data_struct.data2(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data2_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
end
data_struct.data1 = tmp_data1_mat; % subset back into data
data_struct.data2 = tmp_data2_mat;

% make group level plots for main figure
ref_keep = ones(size(data_struct.zmap,1), size(data_struct.zmap,2));
ref_keep = tril(ref_keep);
ref_keep = ref_keep(manual_lag:end, manual_lag:end);

% con - div
f1 = figure;
d1_h = imagesc(reshape(data_struct.zmap(tmp_mask),  [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('analytic zmap')
saveas(f1, [plots_path, 'Reviewer_fig_group_level_interaction_zmap.eps'], 'epsc');
close 

% thresh
d = reshape(data_struct.zmapthresh(tmp_mask), [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]);
tmp = bwconncomp(d);
tmp_idx = cellfun(@(x) length(x)==1, tmp.PixelIdxList);
d(cell2mat(tmp.PixelIdxList(tmp_idx))) = 0;
f1=figure;d1_h=imagesc(d);
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
title('zmap thresh')
colormap(viridis)
caxis manual
caxis([0,3])
c = colorbar;
set(c, 'FontSize',15)
saveas(f1, [plots_path, 'Reviewer_fig_group_level_interaction_zmap_thresh.eps'], 'epsc');
close 

% con > div same context
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data1,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
caxis manual
caxis([-0.03,0.02])
c = colorbar;
set(c, 'FontSize',15)
title('converging > diverging same context ')
saveas(f1, [plots_path, 'Reviewer_fig_group_level_c_d_scon.eps'], 'epsc');
close 

% con > div dcon
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data2,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
caxis manual
caxis([-0.03,0.02])
c = colorbar;
set(c, 'FontSize',15)
title('converging > diverging diff. context')
saveas(f1, [plots_path, 'Reviewer_fig_group_level_c_d_dcon.eps'], 'epsc');
% % saveas(f1, [plots_path, 'Fig4_group_level_diverging.eps'], 'epsc');
close 

%% Converging DIVERGING DCON
cur_ROI = 'HIPP_MERGE_BL';
load([out_dir, cur_ROI,'/', cur_ROI, '_stats.mat'])
% contrasts = {{'diff_start_same_end_scon_cue',% 'same_start_diff_end_scon_cue'}}; not use
data_struct = all_stats(4).stats;

% apply window and figure out masking
manual_lag = 10; % in TRs % This makes P1 activation the zeropoint on the figure
tmp_mask = zeros(size(data_struct.zmap, 1), size(data_struct.zmap,2));
tmp_mask(manual_lag:end, manual_lag:end) = 1;
tmp_mask = logical(tmp_mask);
out_mat_dims = sqrt(length(find(tmp_mask)));

% loop through grabbing the data and replacing
tmp_data1_mat = zeros([size(data_struct.data1,1), out_mat_dims ,out_mat_dims]);
tmp_data2_mat = zeros([size(data_struct.data2,1), out_mat_dims, out_mat_dims]);
for isub = 1:size(data_struct.data1,1)
    % c1
    tmp_data = squeeze(data_struct.data1(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data1_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
    % c2
    tmp_data = squeeze(data_struct.data2(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data2_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
end
data_struct.data1 = tmp_data1_mat; % subset back into data
data_struct.data2 = tmp_data2_mat;

% make group level plots for main figure
ref_keep = ones(size(data_struct.zmap,1), size(data_struct.zmap,2));
ref_keep = tril(ref_keep);
ref_keep = ref_keep(manual_lag:end, manual_lag:end);

% con - div
f1 = figure;
d1_h = imagesc(reshape(data_struct.zmap(tmp_mask),  [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('analytic zmap')
saveas(f1, [plots_path, 'Reviewer_fig_group_level_con_div_dcon_zmap.eps'], 'epsc');
close 

% thresh
d = reshape(data_struct.zmapthresh(tmp_mask), [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]);
tmp = bwconncomp(d);
tmp_idx = cellfun(@(x) length(x)==1, tmp.PixelIdxList);
d(cell2mat(tmp.PixelIdxList(tmp_idx))) = 0;
f1=figure;d1_h=imagesc(d);
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
title('zmap thresh')
colormap(viridis)
% caxis manual
% caxis([0,3])
c = colorbar;
set(c, 'FontSize',15)
saveas(f1, [plots_path, 'Reviewer_fig_group_level_con_div_dcon_zmap_thresh.eps'], 'epsc');
close 

% con diff context
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data1,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
caxis manual
caxis([0,0.035])
c = colorbar;
set(c, 'FontSize',15)
title('converging diff context ')
saveas(f1, [plots_path, 'Reviewer_fig_group_level_con_dcon.eps'], 'epsc');
close 

% div dcon
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data2,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
caxis manual
caxis([0,0.035])
c = colorbar;
set(c, 'FontSize',15)
title('diverging diff. context')
saveas(f1, [plots_path, 'Reviewer_fig_group_level_div_dcon.eps'], 'epsc');
% saveas(f1, [plots_path, 'Fig4_group_level_diverging.eps'], 'epsc');
close 


%% Figure S3
ref_keep = ones(size(data_struct.zmap,1), size(data_struct.zmap,2));
ref_keep = tril(ref_keep);
ref_keep = ref_keep(manual_lag:end, manual_lag:end);

% conv - div
f1=figure('Position', [10,10,1400,1400]);
z=zscore(data_struct.data1 - data_struct.data2, 0, [1]); % zscore difference
for i = 1:size(data_struct.data1,1)
    subplot(6,4,i)
    d1_h=imagesc(squeeze(z(i,:,:)));title(['sub ', num2str(i), ' con>div']);ax = gca;
    set(d1_h,'AlphaData',ref_keep)
    ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
    ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
    caxis manual
    caxis([-3,3])
    colorbar
end
saveas(f1, [plots_path, 'FigS3_indiv_react_con_div.eps'], 'epsc');
close 

%% Figure S4B
cur_ROI = 'V12_BL';
load([out_dir, cur_ROI,'/', cur_ROI, '_stats.mat'])
contrasts = {{'diff_start_same_end_scon_cue', 'same_start_diff_end_scon_cue'}};
data_struct = all_stats(1).stats;

% apply window and figure out masking
manual_lag = 10; % in TRs % This makes P1 activation the zeropoint on the figure
tmp_mask = zeros(size(data_struct.zmap, 1), size(data_struct.zmap,2));
tmp_mask(manual_lag:end, manual_lag:end) = 1;
tmp_mask = logical(tmp_mask);
out_mat_dims = sqrt(length(find(tmp_mask)));

% loop through grabbing the data and replacing
tmp_data1_mat = zeros([size(data_struct.data1,1), out_mat_dims ,out_mat_dims]);
tmp_data2_mat = zeros([size(data_struct.data2,1), out_mat_dims, out_mat_dims]);
for isub = 1:size(data_struct.data1,1)
    % c1
    tmp_data = squeeze(data_struct.data1(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data1_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
    % c2
    tmp_data = squeeze(data_struct.data2(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data2_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
end
data_struct.data1 = tmp_data1_mat; % subset back into data
data_struct.data2 = tmp_data2_mat;

% make group level plots for main figure
ref_keep = ones(size(data_struct.zmap,1), size(data_struct.zmap,2));
ref_keep = tril(ref_keep);
ref_keep = ref_keep(manual_lag:end, manual_lag:end);

% con - div
f1 = figure;
d1_h = imagesc(reshape(data_struct.zmap(tmp_mask),  [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('analytic zmap')
saveas(f1, [plots_path, 'FigS4B_group_level_con_div_zmap.eps'], 'epsc');
close 

% thresh
d = reshape(data_struct.zmapthresh(tmp_mask), [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]);
tmp = bwconncomp(d);
tmp_idx = cellfun(@(x) length(x)==1, tmp.PixelIdxList);
d(cell2mat(tmp.PixelIdxList(tmp_idx))) = 0;
f1=figure;d1_h=imagesc(d);
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
title('zmap thresh')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
saveas(f1, [plots_path, 'FigS4B_group_level_con_div_zmap_thresh.eps'], 'epsc');
close 

% con
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data1,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('converging')
saveas(f1, [plots_path, 'FigS4B_group_level_converging.eps'], 'epsc');
close 

% div 
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data2,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('diverging')
saveas(f1, [plots_path, 'FigS4B_group_level_diverging.eps'], 'epsc');
close 

%% Figure S4D
cur_ROI = 'BA4ap_BL';
load([out_dir, cur_ROI,'/', cur_ROI, '_stats.mat'])
contrasts = {{'diff_start_same_end_scon_cue', 'same_start_diff_end_scon_cue'}};
data_struct = all_stats(1).stats;

% apply window and figure out masking
manual_lag = 10; % in TRs % This makes P1 activation the zeropoint on the figure
tmp_mask = zeros(size(data_struct.zmap, 1), size(data_struct.zmap,2));
tmp_mask(manual_lag:end, manual_lag:end) = 1;
tmp_mask = logical(tmp_mask);
out_mat_dims = sqrt(length(find(tmp_mask)));

% loop through grabbing the data and replacing
tmp_data1_mat = zeros([size(data_struct.data1,1), out_mat_dims ,out_mat_dims]);
tmp_data2_mat = zeros([size(data_struct.data2,1), out_mat_dims, out_mat_dims]);
for isub = 1:size(data_struct.data1,1)
    % c1
    tmp_data = squeeze(data_struct.data1(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data1_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
    % c2
    tmp_data = squeeze(data_struct.data2(isub,:,:)); % subset
    tmp_data = reshape(tmp_data(tmp_mask), [out_mat_dims,out_mat_dims]); % window
    tmp_data2_mat(isub, : ,:) = tmp_data;
    tmp_data = [];
    
end
data_struct.data1 = tmp_data1_mat; % subset back into data
data_struct.data2 = tmp_data2_mat;

% make group level plots for main figure
ref_keep = ones(size(data_struct.zmap,1), size(data_struct.zmap,2));
ref_keep = tril(ref_keep);
ref_keep = ref_keep(manual_lag:end, manual_lag:end);

% con - div
f1 = figure;
d1_h = imagesc(reshape(data_struct.zmap(tmp_mask),  [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('analytic zmap')
saveas(f1, [plots_path, 'FigS4D_group_level_con_div_zmap.eps'], 'epsc');
close 

% thresh
d = reshape(data_struct.zmapthresh(tmp_mask), [size(tmp_mask,1) - 9,size(tmp_mask,2) - 9]);
tmp = bwconncomp(d);
tmp_idx = cellfun(@(x) length(x)==1, tmp.PixelIdxList);
d(cell2mat(tmp.PixelIdxList(tmp_idx))) = 0;
f1=figure;d1_h=imagesc(d);
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
title('zmap thresh')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
saveas(f1, [plots_path, 'FigS4D_group_level_con_div_zmap_thresh.eps'], 'epsc');
close 

% con
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data1,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('converging')
saveas(f1, [plots_path, 'FigS4D_group_level_converging.eps'], 'epsc');
close 

% div 
f1 = figure;d1_h=imagesc(squeeze(mean(data_struct.data2,1)));
set(d1_h,'AlphaData',ref_keep)
ax=gca;
ax.XAxis.TickValues = [1,6,11,16,21,26];ax.XAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.XAxis.FontSize = 15;
ax.YAxis.TickValues = [1,6,11,16,21,26];ax.YAxis.TickLabels = { 'p1', 'p2', 'p3', 'p4', 'p5','End'}; ax.YAxis.FontSize = 15;
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
title('diverging')
saveas(f1, [plots_path, 'FigS4D_group_level_diverging.eps'], 'epsc');
close 
