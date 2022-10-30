%% Colinearity Example Figures For Reviewer

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

%% 
% see this blog for whats in an SPM file - http://andysbrainblog.blogspot.com/2013/10/whats-in-spmmat-file.html
% load data
load('/Volumes/Data/zoocon/images_ffx_par/sub01/sequence_models_fir/Block1/sequence05/SPM.mat')

% probably should do this for all the sequences

% visualize design matrix scaled for display
f1=figure;imagesc(SPM.xX.nKX);
ax=gca;
ax.XAxis.FontSize = 12;ax.XAxis.Label.String = "Parameters";
ax.YAxis.FontSize = 12;ax.YAxis.Label.String = "Timepoints (TRs)";
set(ax,'TickDir','out','box','off')
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
% title('Single Subject Design Matrix')
saveas(f1, [plots_path, 'TR_TR_design_matrix.eps'], 'epsc');
close 


% visualize colinearity
colinear_matrix = triu(squareform(abs(1-pdist(SPM.xX.nKX',  'cosine'))),1);
f2 = figure;imagesc(colinear_matrix);
ax=gca;
ax.XAxis.FontSize = 12;ax.XAxis.Label.String = "Parameters";
ax.YAxis.FontSize = 12;ax.YAxis.Label.String = "Parameters";
colormap(viridis)
c = colorbar;
set(c, 'FontSize',15)
% title('Single Subject Colinearity Matrix')
saveas(f2, [plots_path, 'TR_TR_colinearity_matrix.eps'], 'epsc');
close

% Blue orthogonal vectors
% 1 means colinear
% 0 means orthog
% inbetweem

% get a VIF
%vif() computes variance inflation coefficients  
%VIFs are also the diagonal elements of the inverse of the correlation matrix [1], a convenient result that eliminates the need to set up the various regressions
%[1] Belsley, D. A., E. Kuh, and R. E. Welsch. Regression Diagnostics. Hoboken, NJ: John Wiley & Sons, 1980.
R0 = corrcoef(SPM.xX.nKX); % correlation matrix
V=diag(inv(R0))';

f3 = figure('Position',[10 10 1120 320]);plot(V);
ax = gca;
ax.XAxis.FontSize = 12;ax.XAxis.Label.String = "Parameters";
ax.YAxis.FontSize = 12;ax.YAxis.Label.String = "Variance Inflation Factor";
title('Single Subject Model Variance Inflation Factor')
saveas(f3, [plots_path, 'TR_TR_VIF.eps'], 'epsc');
close 

% 
% absolute_correlation = sum(colinear_matrix, 2);
% figure;plot(absolute_correlation(1:49))

%% loop through everything

% get a VIF
%vif() computes variance inflation coefficients  
%VIFs are also the diagonal elements of the inverse of the correlation matrix [1], a convenient result that eliminates the need to set up the various regressions
%[1] Belsley, D. A., E. Kuh, and R. E. Welsch. Regression Diagnostics. Hoboken, NJ: John Wiley & Sons, 1980.
R0 = corrcoef(SPM.xX.nKX); % correlation matrix
V=diag(inv(R0))';
